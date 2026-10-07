"""Deterministic, dimensionless experiments; no course data or fitted parameters."""

import json
from pathlib import Path

import numpy as np
from scipy.linalg import expm

from .solver import HeatGrid, evolve, exact_discrete

METHODS = ("cn", "be", "rannacher")


def weighted_error(grid, actual, reference):
    return float(np.sqrt(grid.mass @ (actual - reference) ** 2 / grid.mass.sum()))


def pulse_case():
    n = 128
    edges = np.linspace(0, 1, n + 1)
    left = (edges[:-1] + edges[1:]) / 2 < 0.5
    grid = HeatGrid(edges, np.where(left, 1.0, 4.0), np.where(left, 1.0, 2.0))
    initial = np.zeros(n)
    initial[n // 2 - 1] = 1
    return grid, initial, 20 * grid.explicit_bound


def run(output):
    """Write JSON, CSV and three plots, returning the complete metric dictionary."""
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.ticker import NullFormatter

    output = Path(output)
    output.mkdir(parents=True, exist_ok=True)
    plt.rcParams.update(
        {
            "font.size": 10,
            "axes.spines.top": False,
            "axes.spines.right": False,
            "savefig.dpi": 180,
        }
    )
    colors = {"cn": "#bd493c", "be": "#64748b", "rannacher": "#127d86"}
    grid, initial, dt = pulse_case()
    steps = 12
    histories = {m: evolve(grid, initial, dt, steps, m) for m in METHODS}
    exact = exact_discrete(grid, initial, dt)
    final = exact_discrete(grid, initial, dt * steps)
    metrics = {
        "schema_version": 1,
        "pulse": {
            "cells": 128,
            "dt": dt,
            "steps": steps,
            "dt_over_explicit_bound": dt / grid.explicit_bound,
            "methods": {},
        },
        "temporal": [],
        "spatial": [],
    }
    for method, history in histories.items():
        energies = np.array([grid.energy(t) for t in history])
        metrics["pulse"]["methods"][method] = {
            "first_min": float(history[1].min()),
            "all_min": float(history[1:].min()),
            "first_energy_ratio": energies[1] / energies[0],
            "max_energy_increase": float(np.max(np.diff(energies))),
            "relative_heat_drift": float(
                np.max(np.abs(history @ grid.mass - grid.heat(initial))) / grid.heat(initial)
            ),
            "first_time_error": weighted_error(grid, history[1], exact),
            "final_time_error": weighted_error(grid, history[-1], final),
        }
    # Independent dense exponential on a small heterogeneous, nonuniform mesh.
    small = HeatGrid([0, 0.07, 0.22, 0.5, 0.8, 1], [1, 3, 0.4, 2, 1], [2, 1, 3, 2, 1])
    k = np.column_stack([small.stiffness_action(e) for e in np.eye(5)])
    initial_small = np.array([1.0, 0, 4, 2, 1])
    reference = expm(-0.17 * k / small.mass[:, None]) @ initial_small
    metrics["matrix_exponential_max_error"] = float(
        np.max(np.abs(exact_discrete(small, initial_small, 0.17) - reference))
    )
    # Fixed spatial grid: isolate time discretization with an eigenvector oracle.
    time_final = 0.02
    time_reference = exact_discrete(grid, initial, time_final)
    for count in [20, 40, 80, 160, 320]:
        row = {"steps": count, "dt": time_final / count}
        for m in METHODS:
            row[m] = weighted_error(
                grid, evolve(grid, initial, row["dt"], count, m)[-1], time_reference
            )
        metrics["temporal"].append(row)
    metrics["temporal_last_orders"] = {
        m: float(np.log2(metrics["temporal"][-2][m] / metrics["temporal"][-1][m]))
        for m in METHODS
    }
    # Exact time propagation on uniform meshes: isolate continuum spatial error.
    for n in [16, 32, 64, 128, 256]:
        uniform = HeatGrid(np.linspace(0, 1, n + 1), np.ones(n), np.ones(n))
        cosine_avg = np.cos(np.pi * uniform.centres) * np.sinc(1 / (2 * n))
        analytic = 2 + np.exp(-(np.pi**2) * 0.1) * cosine_avg
        discrete = exact_discrete(uniform, 2 + cosine_avg, 0.1)
        metrics["spatial"].append(
            {"cells": n, "cell_average_l2_error": weighted_error(uniform, discrete, analytic)}
        )
    errors = [r["cell_average_l2_error"] for r in metrics["spatial"]]
    metrics["spatial_orders"] = np.log2(np.array(errors[:-1]) / errors[1:]).tolist()
    # Complete histories are recorded, not only the values selected for the paper.
    profiles = np.column_stack(
        [
            grid.centres,
            initial,
            exact,
            *(histories[m][1] for m in METHODS),
            final,
            *(histories[m][-1] for m in METHODS),
        ]
    )
    np.savetxt(
        output / "profiles.csv",
        profiles,
        delimiter=",",
        comments="",
        header="x,initial,exact_first,cn_first,be_first,rannacher_first,"
        "exact_final,cn_final,be_final,rannacher_final",
    )
    for m, history in histories.items():
        np.savetxt(
            output / f"history_{m}.csv",
            np.column_stack([np.arange(steps + 1) * dt, history]),
            delimiter=",",
            comments="",
            header="time," + ",".join(f"cell_{i}" for i in range(128)),
        )
    (output / "metrics.json").write_text(
        json.dumps(metrics, indent=2, allow_nan=False) + "\n", encoding="utf-8", newline="\n"
    )
    fig, axes = plt.subplots(1, 2, figsize=(9, 3.2), layout="constrained")
    for ax, step, truth in zip(axes, [1, steps], [exact, final], strict=True):
        ax.plot(grid.centres, truth, color="#172638", lw=2, label="exact semidiscrete")
        for m in METHODS:
            ax.plot(
                grid.centres,
                histories[m][step],
                color=colors[m],
                label=m,
                linestyle="--" if m == "cn" else "-",
            )
        ax.axhline(0, color="black", lw=0.5)
        ax.axvline(0.5, color="grey", lw=0.6, linestyle=":")
        ax.set(
            xlim=(0.35, 0.65),
            xlabel="Position x",
            ylabel="Temperature excess",
            title=f"After {step} macro step{'s' if step > 1 else ''}",
        )
    axes[1].legend(fontsize=8)
    fig.savefig(output / "pulse.png")
    plt.close(fig)
    fig, axes = plt.subplots(1, 2, figsize=(9, 3), layout="constrained")
    for m in METHODS:
        history = histories[m]
        axes[0].plot(
            np.arange(steps + 1),
            [grid.energy(t) / grid.energy(initial) for t in history],
            "o-",
            ms=3,
            color=colors[m],
            label=m,
        )
        axes[1].plot(
            np.arange(steps + 1), history.min(axis=1), "o-", ms=3, color=colors[m], label=m
        )
    axes[0].set(xlabel="Macro step", ylabel="Quadratic energy / initial", yscale="log")
    axes[1].set(xlabel="Macro step", ylabel="Minimum temperature excess")
    axes[0].legend()
    fig.savefig(output / "diagnostics.png")
    plt.close(fig)
    fig, axes = plt.subplots(1, 2, figsize=(9, 3.2), layout="constrained")
    for m in METHODS:
        axes[0].loglog(
            [r["dt"] for r in metrics["temporal"]],
            [r[m] for r in metrics["temporal"]],
            "o-",
            color=colors[m],
            label=m,
        )
    axes[0].set(xlabel="Time step (fixed 128 cells)", ylabel="Time error, mass RMS")
    axes[0].legend()
    axes[1].loglog([1 / r["cells"] for r in metrics["spatial"]], errors, "o-", color="#127d86")
    axes[1].set_xticks([1 / 256, 1 / 64, 1 / 16], ["1/256", "1/64", "1/16"])
    axes[1].xaxis.set_minor_formatter(NullFormatter())
    axes[1].set(
        xlabel="Cell width (exact time propagation)", ylabel="Continuum cell-average error"
    )
    fig.savefig(output / "convergence.png")
    plt.close(fig)
    return metrics
