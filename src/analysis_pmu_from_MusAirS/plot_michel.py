#!/usr/bin/env python3
"""Plot independent pi-, K-, and KLong-produced muon Michel angular distributions."""
import argparse
from pathlib import Path
import sys

from plot_pmu import read_data

CURVES = (
    ("piMuon", "r", "-", r"$\mu^-(\pi^-)$"),
    ("KMuon", "g", "-", r"$\mu^-(K^-)$"),
    ("KLongMuon", "b", "-", r"$\mu^-(K_L^0)$"),
)


def polar_curve(rows, points=500):
    """Mirror cos(theta) bins; this is a radial display, not dN/dtheta."""
    import numpy as np
    from scipy.interpolate import CubicSpline

    cosine = np.array([row[1] for row in rows])
    counts = np.array([row[2] for row in rows])
    if np.any(np.abs(cosine) >= 1):
        raise ValueError("Michel bin centers must be strictly inside (-1, 1)")
    theta = np.arccos(cosine)
    angles = np.concatenate((theta, 2 * np.pi - theta))
    values = np.concatenate((counts, counts))
    order = np.argsort(angles)
    angles, values = angles[order], values[order]
    # Both nearest points have identical mirrored values. Close the curve with
    # that value and periodic derivatives; do not duplicate any bin-center angle.
    angles = np.concatenate(([0], angles, [2 * np.pi]))
    values = np.concatenate(([values[0]], values, [values[0]]))
    spline = CubicSpline(angles, values, bc_type="periodic")
    dense = np.linspace(0, 2 * np.pi, points)
    radial = spline(dense)
    if np.any(radial < 0):
        print("Warning: negative spline overshoot clipped to zero; increase sampling.", file=sys.stderr)
        radial = np.maximum(radial, 0)
    return dense, radial


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("data", type=Path, help="michelData.csv from generate_michel.C")
    parser.add_argument("-o", "--output", type=Path, help="default: beside CSV, michelComponents.pdf")
    parser.add_argument("--overwrite", action="store_true")
    parser.add_argument("--show", action="store_true")
    args = parser.parse_args()
    output = args.output or args.data.with_name("michelComponents.pdf")
    try:
        if output.resolve() == args.data.resolve():
            raise ValueError("Output must not replace input CSV")
        if output.exists() and not args.overwrite:
            raise ValueError(f"Refusing to overwrite {output}; use --overwrite")
        data = read_data(args.data, names=[item[0] for item in CURVES], normalized=False)
        import matplotlib
        if not args.show:
            matplotlib.use("Agg")
        import matplotlib.pyplot as plt

        figure, axes = plt.subplots(subplot_kw={"projection": "polar"}, figsize=(4, 4))
        maximum = 0
        for name, color, style, label in CURVES:
            theta, radial = polar_curve(data[name])
            maximum = max(maximum, float(radial.max()))
            axes.plot(theta, radial, color=color, linestyle=style, linewidth=1.5, label=label)
        axes.set_yticklabels([])
        axes.tick_params(labelsize=12)
        axes.set_ylim(0, 1.1 * maximum)
        axes.text(0.7, 0.3 * maximum, r"$\vec{p}(\mu)$", ha="center", va="center", fontsize=18)
        axes.annotate("", xy=(0, 0.7 * maximum), xytext=(0, 0),
                      arrowprops={"arrowstyle": "->", "linewidth": 5, "color": "red"})
        axes.legend(loc="upper right", bbox_to_anchor=(1.15, -0.1), ncol=3,
                    fancybox=True, shadow=True, fontsize=12, handlelength=0.8)
        output.parent.mkdir(parents=True, exist_ok=True)
        figure.savefig(output, bbox_inches="tight")
        if args.show:
            plt.show()
        plt.close(figure)
        print(f"Saved {output}")
    except (OSError, ValueError, ImportError) as exc:
        print(f"Michel plotting failed: {exc}", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
