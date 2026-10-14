"""Reproduce the second (1-sigma interval) figure of dev/2601reviseArticle/ratio5percent.ipynb."""

import argparse
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from scipy.interpolate import UnivariateSpline
from scipy.optimize import brentq
import uproot


def find_intersections(x, y, delta_chi2=1.0):
    x, y = np.asarray(x, dtype=float), np.asarray(y, dtype=float)
    if (len(x) < 4 or x.shape != y.shape or not np.all(np.isfinite(x))
            or not np.all(np.isfinite(y)) or np.any(np.diff(x) <= 0)):
        raise ValueError("Need >=4 finite points with strictly increasing ratios")
    if not np.isfinite(delta_chi2) or delta_chi2 <= 0:
        raise ValueError("Delta chi-square threshold must be positive")
    spline = UnivariateSpline(x, y, s=0)
    target = float(y.min() + delta_chi2)
    grid = np.linspace(x.min(), x.max(), 1000)
    differences = spline(grid) - target
    # Preserve the notebook's midpoint markers for visual reproduction.
    midpoints, precise = [], []
    for left, right, fl, fr in zip(grid[:-1], grid[1:], differences[:-1], differences[1:]):
        if fl * fr < 0:
            midpoints.append(float((left + right) / 2))
            precise.append(float(brentq(lambda ratio: float(spline(ratio) - target), left, right)))
        elif fl == 0:
            midpoints.append(float(left))
            precise.append(float(left))
    if differences[-1] == 0:
        midpoints.append(float(grid[-1]))
        precise.append(float(grid[-1]))
    if len(precise) != 2:
        raise ValueError(f"Expected two threshold crossings; found {len(precise)}. Widen the scan range.")
    return spline, target, midpoints, precise


def plot_scan(input_path, output_dir, delta_chi2=1.0):
    with uproot.open(input_path) as source:
        x = source["tree"]["ratio"].array(library="np")
        y = source["tree"]["result"].array(library="np")
    spline, target, markers, precise = find_intersections(x, y, delta_chi2)
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    fig, ax = plt.subplots(figsize=(8, 6))
    x_fill = np.linspace(markers[0], markers[1], 100)
    ax.fill_between(x_fill, spline(x_fill), target, color="red", alpha=0.2,
                    label=rf"$\Delta\chi^2={delta_chi2:g}$ Interval")
    ax.plot(x, y, label="ratio vs result")
    ax.axhline(target, color="r", linestyle="--", label=f"Target Y = {target}")
    ax.scatter(markers, spline(markers), color="b", label="Intersection Points")
    ax.tick_params(labelsize=20)
    ax.set_xlabel(r"K-$\pi$ Ratio", fontsize=22)
    ax.set_ylabel(r"$\Delta \chi^2$", fontsize=22)
    ax.set_ylim(0, 5)
    ax.set_xlim(0.04, 0.06)
    ax.grid()
    for extension in ("pdf", "png"):
        fig.savefig(output_dir / f"measure_ratio.{extension}", bbox_inches="tight", dpi=150)
    plt.close(fig)
    summary = {
        "input": str(Path(input_path).resolve()),
        "minimum_chi2": float(y.min()),
        "best_fit_ratio_grid": float(x[np.argmin(y)]),
        "delta_chi2_threshold": delta_chi2,
        "target_chi2": target,
        "notebook_midpoint_intersections": markers,
        "refined_intersections": precise,
    }
    (output_dir / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    return summary


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("input", type=Path, help="ROOT file with tree/ratio and tree/result")
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--delta-chi2", type=float, default=1.0)
    args = parser.parse_args()
    print(json.dumps(plot_scan(args.input, args.output_dir, args.delta_chi2), indent=2))
