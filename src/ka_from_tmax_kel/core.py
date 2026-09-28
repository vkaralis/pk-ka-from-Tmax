"""Numerical routines for estimating the absorption rate constant."""

from __future__ import annotations

from dataclasses import dataclass
from math import isfinite, log, log1p


@dataclass(frozen=True)
class EstimationResult:
    """Result returned by :func:`estimate_ka`."""

    ka: float
    tmax: float
    kel: float
    iterations: int
    residual: float


def _positive(value: float, name: str) -> float:
    value = float(value)
    if not isfinite(value) or value <= 0:
        raise ValueError(f"{name} must be a finite number greater than zero")
    return value


def tmax_from_rates(ka: float, kel: float) -> float:
    """Return Tmax for positive first-order absorption and elimination rates.

    ``ka`` and ``kel`` must use the same inverse-time unit. The returned Tmax
    uses the corresponding time unit. At ``ka == kel`` the continuous limit,
    ``1 / kel``, is returned.
    """

    ka = _positive(ka, "ka")
    kel = _positive(kel, "kel")
    relative_difference = (ka - kel) / kel
    if abs(relative_difference) < 1e-8:
        # log1p avoids cancellation close to Ka = Kel.
        return log1p(relative_difference) / (kel * relative_difference) if relative_difference else 1.0 / kel
    return log(ka / kel) / (ka - kel)


def estimate_ka(
    tmax: float,
    kel: float,
    *,
    rtol: float = 1e-10,
    max_iterations: int = 200,
) -> EstimationResult:
    """Estimate the unique positive Ka compatible with ``tmax`` and ``kel``.

    The equation is solved by bisection, without third-party dependencies.
    When ``tmax == 1 / kel``, the analytical limiting solution is ``Ka = Kel``.
    """

    tmax = _positive(tmax, "tmax")
    kel = _positive(kel, "kel")
    rtol = _positive(rtol, "rtol")
    if not isinstance(max_iterations, int) or max_iterations <= 0:
        raise ValueError("max_iterations must be a positive integer")

    limit_tmax = 1.0 / kel
    if abs(tmax - limit_tmax) <= rtol * max(tmax, limit_tmax):
        return EstimationResult(kel, tmax, kel, 0, limit_tmax - tmax)

    # Tmax decreases monotonically as Ka increases. Select the appropriate
    # side of Ka = Kel, then expand the finite bracket when necessary.
    if tmax < limit_tmax:
        low, high = kel, 2.0 * kel
        while tmax_from_rates(high, kel) > tmax:
            high *= 2.0
            if not isfinite(high):
                raise ArithmeticError("could not bracket Ka; check the input scale")
    else:
        high = kel
        low = kel * 0.5
        while tmax_from_rates(low, kel) < tmax:
            low *= 0.5
            if low == 0.0:
                raise ArithmeticError("could not bracket Ka; check the input scale")

    midpoint = (low + high) / 2.0
    for iteration in range(1, max_iterations + 1):
        midpoint = (low + high) / 2.0
        calculated = tmax_from_rates(midpoint, kel)
        residual = calculated - tmax
        if abs(residual) <= rtol * tmax or (high - low) <= rtol * midpoint:
            return EstimationResult(midpoint, tmax, kel, iteration, residual)
        if calculated > tmax:
            low = midpoint
        else:
            high = midpoint

    raise RuntimeError(f"Ka did not converge after {max_iterations} iterations")

