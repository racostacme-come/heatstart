import numpy as np
import pytest
from numpy.testing import assert_allclose
from scipy.linalg import expm

from heatstart import HeatGrid, evolve, exact_discrete


@pytest.fixture
def grid():
    return HeatGrid([0, 0.1, 0.3, 0.7, 1], [1, 2, 0.5, 3], [2, 1, 4, 2])


def dense(grid):
    n = grid.mass.size
    return np.column_stack([grid.stiffness_action(np.eye(n)[:, i]) for i in range(n)])


@pytest.mark.parametrize("method", ["be", "cn", "rannacher"])
def test_constant_and_conservation(grid, method):
    assert_allclose(evolve(grid, np.full(4, 7.0), 0.02, 50, method), 7, atol=2e-13)
    t = evolve(grid, [1, 0, 4, 2], 0.1, 20, method)
    assert_allclose(t @ grid.mass, grid.heat(t[0]), atol=2e-13)
    assert np.max(np.diff([grid.energy(row) for row in t])) < 1e-13


@pytest.mark.parametrize("method,theta", [("be", 1), ("cn", 0.5)])
def test_step_against_dense_solve(grid, method, theta):
    initial, dt = np.array([1, 3, 2, 0]), 0.037
    k, m = dense(grid), np.diag(grid.mass)
    expected = np.linalg.solve(m + theta * dt * k, (m - (1 - theta) * dt * k) @ initial)
    assert_allclose(evolve(grid, initial, dt, 1, method)[-1], expected, atol=2e-14)


def test_four_half_steps_then_cn(grid):
    initial, dt = [0, 1, 4, 1], 0.02
    startup = evolve(grid, initial, dt / 2, 4, "be")[-1]
    expected = evolve(grid, startup, dt, 3, "cn")[-1]
    assert_allclose(evolve(grid, initial, dt, 5, "rannacher")[-1], expected)


def test_exact_against_matrix_exponential(grid):
    initial = np.array([3, 1, 2, 7])
    expected = expm(-0.13 * dense(grid) / grid.mass[:, None]) @ initial
    assert_allclose(exact_discrete(grid, initial, 0.13), expected, atol=2e-13)
    assert_allclose(exact_discrete(grid, initial, 0), initial, atol=2e-14)


def test_face_resistance_and_flux(grid):
    dx = np.diff(grid.edges)
    resistance = dx[:-1] / (2 * grid.conductivity[:-1]) + dx[1:] / (2 * grid.conductivity[1:])
    assert_allclose(grid.conductance, 1 / resistance)
    k = dense(grid)
    assert_allclose(k, k.T)
    assert_allclose(k @ np.ones(4), 0, atol=1e-14)
    assert np.linalg.eigvalsh(k)[0] > -1e-13


def test_energy_identity(grid):
    dt = 0.12
    old, new = evolve(grid, [1, 3, 7, 2], dt, 1)
    mid = (old + new) / 2
    assert_allclose(grid.energy(new) - grid.energy(old), -dt * mid @ grid.stiffness_action(mid))


def test_stable_cn_can_be_negative():
    grid = HeatGrid(np.linspace(0, 1, 41), np.ones(40), np.ones(40))
    initial = np.zeros(40)
    initial[20] = 1
    dt = 20 * grid.explicit_bound
    cn = evolve(grid, initial, dt, 1)[-1]
    assert cn.min() < -0.5
    assert grid.energy(cn) < grid.energy(initial)
    for method in ["be", "rannacher"]:
        assert evolve(grid, initial, dt, 1, method)[-1].min() >= 0


@pytest.mark.parametrize("method,order", [("be", 1), ("cn", 2), ("rannacher", 2)])
def test_time_order(grid, method, order):
    initial = np.array([1, 2, 0, 3])
    exact = exact_discrete(grid, initial, 0.2)
    errors = [
        np.linalg.norm(evolve(grid, initial, 0.2 / n, n, method)[-1] - exact)
        for n in [160, 320]
    ]
    assert np.log2(errors[0] / errors[1]) > order - 0.12


def test_continuum_cell_average_convergence():
    errors = []
    for n in [20, 40, 80]:
        grid = HeatGrid(np.linspace(0, 1, n + 1), np.ones(n), np.ones(n))
        cosine = np.cos(np.pi * grid.centres) * np.sinc(1 / (2 * n))
        exact = 2 + np.exp(-(np.pi**2) * 0.1) * cosine
        errors.append(np.max(np.abs(exact_discrete(grid, 2 + cosine, 0.1) - exact)))
    assert np.all(np.log2(np.array(errors[:-1]) / errors[1:]) > 1.99)


@pytest.mark.parametrize(
    "edges,k,c",
    [
        ([0, 1], [1], [1]),
        ([0, 1, 1], [1, 1], [1, 1]),
        ([0, 2, 1], [1, 1], [1, 1]),
        ([0, 1, 2], [0, 1], [1, 1]),
        ([0, 1, 2], [1, 1], [-1, 1]),
        ([0, np.nan, 2], [1, 1], [1, 1]),
        ([0, 1, 2], [1], [1, 1]),
        ([0, 1, 2], [1, 1], [np.inf, 1]),
    ],
)
def test_bad_grid(edges, k, c):
    with pytest.raises(ValueError):
        HeatGrid(edges, k, c)


@pytest.mark.parametrize(
    "kwargs",
    [
        {"dt": 0},
        {"dt": -1},
        {"dt": np.nan},
        {"dt": np.inf},
        {"steps": -1},
        {"steps": 2.5},
        {"steps": True},
        {"startup_steps": -1},
        {"method": "unknown"},
        {"initial": [1, 2]},
        {"initial": [1, 2, 3, np.nan]},
    ],
)
def test_bad_evolution(grid, kwargs):
    args = dict(initial=[1, 2, 3, 4], dt=0.01, steps=3)
    args.update(kwargs)
    with pytest.raises(ValueError):
        evolve(grid, **args)


def test_readonly_and_no_input_mutation(grid):
    t = np.array([1.0, 2, 3, 4])
    copy = t.copy()
    evolve(grid, t, 0.01, 4)
    assert_allclose(t, copy)
    with pytest.raises(ValueError):
        grid.mass[0] = 0
    assert evolve(grid, t, 0.01, 0).shape == (1, 4)


@pytest.mark.parametrize("method", ["be", "cn", "rannacher"])
def test_analytical_two_cell_mode(method):
    grid = HeatGrid([0, 0.3, 1], [2, 4], [3, 1])
    initial = np.array([4.0, 1.0])
    lam = grid.conductance[0] * np.sum(1 / grid.mass)
    mean = grid.heat(initial) / grid.mass.sum()
    dt, steps = 0.02, 7
    z = lam * dt
    cn = (1 - z / 2) / (1 + z / 2)
    amplification = {
        "be": (1 + z) ** -steps,
        "cn": cn**steps,
        "rannacher": (1 + z / 2) ** -4 * cn ** (steps - 2),
    }[method]
    assert_allclose(
        evolve(grid, initial, dt, steps, method)[-1],
        mean + amplification * (initial - mean),
        atol=3e-14,
    )
    assert_allclose(
        exact_discrete(grid, initial, steps * dt),
        mean + np.exp(-lam * steps * dt) * (initial - mean),
        atol=3e-14,
    )


def test_cn_sufficient_positivity_and_be_large_step(grid):
    for method, dt in [("cn", 2 * grid.explicit_bound), ("be", 10.0)]:
        columns = np.column_stack([evolve(grid, e, dt, 1, method)[-1] for e in np.eye(4)])
        assert columns.min() >= -1e-14
        assert_allclose(columns @ np.ones(4), 1, atol=3e-13)


@pytest.mark.parametrize("time", [-1, np.inf, np.nan])
def test_invalid_exact_time(grid, time):
    with pytest.raises(ValueError):
        exact_discrete(grid, [1, 2, 3, 4], time)
