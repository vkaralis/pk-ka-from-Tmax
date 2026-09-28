"""Command-line interface."""

from __future__ import annotations

import argparse
from pathlib import Path

from .core import estimate_ka, tmax_from_rates


def _parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Estimate Ka from Tmax and Kel using the one-compartment oral PK equation."
    )
    parser.add_argument("tmax", type=float, help="observed Tmax (for example, hours)")
    parser.add_argument("kel", type=float, help="elimination rate constant (same inverse-time unit)")
    parser.add_argument("--plot", type=Path, metavar="FILE", help="save a diagnostic plot (requires matplotlib)")
    parser.add_argument("--precision", type=int, default=8, help="significant digits in output (default: 8)")
    return parser


def _save_plot(path: Path, tmax: float, kel: float, ka: float) -> None:
    try:
        import matplotlib.pyplot as plt
    except ImportError as exc:
        raise SystemExit("Plotting requires: pip install 'ka-from-tmax-kel[plot]'") from exc

    lower = min(ka, kel) * 0.2
    upper = max(ka, kel) * 3.0
    points = 300
    x = [lower * (upper / lower) ** (i / (points - 1)) for i in range(points)]
    y = [tmax_from_rates(value, kel) for value in x]
    fig, ax = plt.subplots()
    ax.plot(x, y, label="Predicted Tmax")
    ax.axhline(tmax, color="tab:red", linestyle="--", label="Observed Tmax")
    ax.scatter([ka], [tmax], color="black", zorder=3, label=f"Ka = {ka:.6g}")
    ax.set(xlabel="Ka", ylabel="Tmax", xscale="log", title="Ka estimation from Tmax and Kel")
    ax.grid(True, alpha=0.3)
    ax.legend()
    fig.tight_layout()
    fig.savefig(path, dpi=160)
    plt.close(fig)


def main(argv: list[str] | None = None) -> int:
    args = _parser().parse_args(argv)
    if args.precision < 1:
        raise SystemExit("--precision must be at least 1")
    try:
        result = estimate_ka(args.tmax, args.kel)
    except (ValueError, ArithmeticError, RuntimeError) as exc:
        raise SystemExit(f"Error: {exc}") from exc

    print(f"Ka = {result.ka:.{args.precision}g}")
    print(f"Residual Tmax = {result.residual:.3g}")
    print(f"Iterations = {result.iterations}")
    if args.plot:
        _save_plot(args.plot, args.tmax, args.kel, result.ka)
        print(f"Plot saved to {args.plot}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

