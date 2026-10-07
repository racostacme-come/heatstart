"""Cell-centred finite volumes and factored, banded theta integration.

Unit cross-sectional area, insulated ends, static positive material properties.
Temperatures represent cell averages approximated at cell centres.
"""

from dataclasses import dataclass
from numbers import Integral

import numpy as np
from scipy.linalg import cho_solve_banded, cholesky_banded, eigh_tridiagonal


def _vector(value, name, size=None, positive=False):
    arr = np.array(value, dtype=float, copy=True)
    if arr.ndim != 1 or (size is not None and arr.size != size):
        raise ValueError(f"{name} must be a one-dimensional array of length {size}")
    if not np.all(np.isfinite(arr)) or (positive and np.any(arr <= 0)):
        raise ValueError(f"{name} must be finite" + (" and positive" if positive else ""))
    return arr


@dataclass(frozen=True, init=False)
class HeatGrid:
    """Read-only mesh data; conductivity and volumetric capacity are cell arrays.

    Edges must be strictly increasing. At least two cells are required.
    In SI units: edges [m], conductivity [W/(m K)], capacity [J/(m^3 K)].
    """

    edges: np.ndarray
    conductivity: np.ndarray
    capacity: np.ndarray
    mass: np.ndarray
    conductance: np.ndarray
    diagonal: np.ndarray

    def __init__(self, edges, conductivity, capacity):
        edges = _vector(edges, "edges")
        if edges.size < 3:
            raise ValueError("at least two cells are required")
        dx = np.diff(edges)
        if not np.all(np.isfinite(dx)) or np.any(dx <= 0):
            raise ValueError("edges must increase with finite positive widths")
        k = _vector(conductivity, "conductivity", dx.size, positive=True)
        c = _vector(capacity, "capacity", dx.size, positive=True)
        with np.errstate(over="ignore", under="ignore", divide="ignore", invalid="ignore"):
            mass = c * dx
            g = 1 / (dx[:-1] / (2 * k[:-1]) + dx[1:] / (2 * k[1:]))
            diagonal = np.r_[g, 0.0] + np.r_[0.0, g]
        if any(not np.all(np.isfinite(v)) or np.any(v <= 0) for v in (mass, g, diagonal)):
            raise ValueError("derived coefficients overflow or underflow; rescale the problem")
        for name, value in (
            ("edges", edges),
            ("conductivity", k),
            ("capacity", c),
            ("mass", mass),
            ("conductance", g),
            ("diagonal", diagonal),
        ):
            value.setflags(write=False)
            object.__setattr__(self, name, value)

    @property
    def centres(self):
        return self.edges[:-1] + np.diff(self.edges) / 2

    @property
    def explicit_bound(self):
        """Sufficient forward-Euler positivity bound, min_i M_i/K_ii."""
        return float(np.min(self.mass / self.diagonal))

    def stiffness_action(self, temperature):
        """Return K T using paired face fluxes; both exterior fluxes vanish."""
        temperature = _vector(temperature, "temperature", self.mass.size)
        flux = self.conductance * (temperature[:-1] - temperature[1:])
        out = np.zeros_like(temperature)
        out[:-1] += flux
        out[1:] -= flux
        return out

    def heat(self, temperature):
        """Heat relative to the chosen temperature datum, per unit area."""
        return float(self.mass @ _vector(temperature, "temperature", self.mass.size))

    def energy(self, temperature):
        """Quadratic diffusion functional; this is not thermodynamic heat."""
        t = _vector(temperature, "temperature", self.mass.size)
        return float(0.5 * self.mass @ t**2)


class _ThetaStep:
    def __init__(self, grid, dt, theta):
        self.grid, self.dt, self.theta = grid, dt, theta
        with np.errstate(over="ignore", invalid="ignore"):
            band = np.zeros((2, grid.mass.size))
            band[0] = grid.mass + theta * dt * grid.diagonal
            band[1, :-1] = -theta * dt * grid.conductance
        if not np.all(np.isfinite(band)):
            raise ValueError("time-step matrix overflow; rescale the problem")
        self.factor = cholesky_banded(band, lower=True)

    def __call__(self, t):
        rhs = self.grid.mass * t - (1 - self.theta) * self.dt * self.grid.stiffness_action(t)
        out = cho_solve_banded((self.factor, True), rhs)
        if not np.all(np.isfinite(out)):
            raise FloatingPointError("nonfinite temperature; rescale the problem")
        return out


def evolve(grid, initial, dt, steps, method="cn", startup_steps=2):
    """Return temperatures at times 0, dt, ..., steps*dt (shape steps+1, cells).

    Methods: 'be' (backward Euler), 'cn' (Crank--Nicolson), 'rannacher'.
    Rannacher replaces each of the first startup_steps macro intervals by two
    backward-Euler half-steps, then uses CN. Startup is performed once per call.
    No universal positivity guarantee is claimed for subsequent CN intervals.
    """
    if not np.isscalar(dt) or not np.isfinite(dt) or dt <= 0:
        raise ValueError("dt must be finite and positive")
    for name, value in (("steps", steps), ("startup_steps", startup_steps)):
        if isinstance(value, bool) or not isinstance(value, Integral) or value < 0:
            raise ValueError(f"{name} must be a nonnegative integer")
    if method not in ("be", "cn", "rannacher"):
        raise ValueError("method must be be, cn, or rannacher")
    t = _vector(initial, "initial", grid.mass.size)
    out = np.empty((steps + 1, t.size))
    out[0] = t
    if steps == 0:
        return out
    full = _ThetaStep(grid, dt, 1.0 if method == "be" else 0.5)
    half = _ThetaStep(grid, dt / 2, 1.0) if method == "rannacher" else None
    for i in range(steps):
        t = half(half(t)) if half is not None and i < startup_steps else full(t)
        out[i + 1] = t
    return out


def exact_discrete(grid, initial, time):
    """Dense eigenvector reference exp(-M^-1 K t)T0; O(N^2) storage.

    This removes time integration error only, not spatial discretization error.
    """
    if not np.isscalar(time) or not np.isfinite(time) or time < 0:
        raise ValueError("time must be finite and nonnegative")
    t = _vector(initial, "initial", grid.mass.size)
    root = np.sqrt(grid.mass)
    values, vectors = eigh_tridiagonal(
        grid.diagonal / grid.mass, -grid.conductance / (root[:-1] * root[1:])
    )
    # The analytic generator is positive semidefinite; zero may round negative.
    values = np.maximum(values, 0)
    return (vectors @ (np.exp(-time * values) * (vectors.T @ (root * t)))) / root
