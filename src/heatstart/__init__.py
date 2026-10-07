"""Conservative one-dimensional insulated heat diffusion."""

from .solver import HeatGrid, evolve, exact_discrete

__all__ = ["HeatGrid", "evolve", "exact_discrete"]
