"""Estimate Ka from Tmax and Kel for a one-compartment oral PK model."""

from .core import EstimationResult, estimate_ka, tmax_from_rates

__all__ = ["EstimationResult", "estimate_ka", "tmax_from_rates"]
__version__ = "1.0.0"

