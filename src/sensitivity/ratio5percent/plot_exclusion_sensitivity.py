"""Redraw michel_sensitivity_plot.ipynb with light-blue shading and violet text."""

import argparse
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import uproot

REPO_ROOT = Path(__file__).resolve().parents[3]
DEFAULT_INPUT = REPO_ROOT / "dev/2601reviseArticle/sensitivity_revise/eigen_sensitiviy_output.root"
OUTPUT_NAME = "241013kpiRatioNSensitivity"
SHADE_COLOR = "#f6f6ffff"
TEXT_COLOR = "#8a8af5ff"


def build_figure(ratios, event_counts):
    ratios, event_counts = np.asarray(ratios, dtype=float), np.asarray(event_counts, dtype=float)
    if (ratios.ndim != 1 or ratios.shape != event_counts.shape or len(ratios) < 2
            or not np.all(np.isfinite(ratios)) or not np.all(np.isfinite(event_counts))
            or np.any(ratios <= 0) or np.any(event_counts <= 0)):
        raise ValueError("Need matching arrays of >=2 finite, positive ratios and event counts")
    order = np.argsort(ratios)
    ratios, event_counts = ratios[order], event_counts[order]
    if np.any(np.diff(ratios) <= 0):
        raise ValueError("Ratios must be unique")
    x_pos, x_delta = 0.10, 0.01
    if ratios[0] > x_pos - x_delta or ratios[-1] < x_pos + x_delta:
        raise ValueError("Ratios must cover [0.09, 0.11] for the notebook annotation")

    fig, ax = plt.subplots(figsize=(6.4, 4.8))
    ax.plot(ratios, event_counts, "--", label=r"$N_\mathrm{evt}$ to $3\sigma$")
    ax.fill_between(ratios, 0, event_counts, color=SHADE_COLOR,
                    label="Kaon not significant")
    ax.set_yscale("log")
    ax.set_ylim(3e4, 2e7)
    ax.set_xlim(0.01, 0.20)
    ax.tick_params(labelsize=14)
    # Preserve the notebook label and fractional x coordinates verbatim.
    ax.set_xlabel(r"$K$-$\pi$ ratio (%)", fontsize=16)
    ax.set_ylabel(r"Events number ($N_\mathrm{evt}$)", fontsize=16)

    # Finalize the log-axis transform before measuring the screen-space angle.
    fig.canvas.draw()
    y_pos = np.interp(x_pos, ratios, event_counts)
    y_near = np.interp([x_pos - x_delta, x_pos + x_delta], ratios, event_counts)
    p1 = ax.transData.transform((x_pos - x_delta, y_near[0]))
    p2 = ax.transData.transform((x_pos + x_delta, y_near[1]))
    angle = np.degrees(np.arctan2(p2[1] - p1[1], p2[0] - p1[0]))
    ax.text(x_pos, y_pos * 0.5, "Kaon not significant", rotation=angle,
            rotation_mode="anchor", ha="center", va="center", fontsize=16,
            fontweight="bold", color=TEXT_COLOR)
    return fig


def plot_exclusion(input_path=DEFAULT_INPUT, output_dir=Path(__file__).resolve().parent / "results"):
    with uproot.open(input_path) as source:
        ntuple = source["sensitivityResult"]
        ratios = ntuple["ratio"].array(library="np")
        event_counts = ntuple["targetStatistics"].array(library="np")
    fig = build_figure(ratios, event_counts)
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    try:
        for extension in ("pdf", "png"):
            fig.savefig(output_dir / f"{OUTPUT_NAME}.{extension}", bbox_inches="tight", dpi=150)
    finally:
        plt.close(fig)
    return {"input": str(Path(input_path).resolve()), "points": len(ratios),
            "pdf": str((output_dir / f"{OUTPUT_NAME}.pdf").resolve()),
            "png": str((output_dir / f"{OUTPUT_NAME}.png").resolve()),
            "shade_color": SHADE_COLOR, "text_color": TEXT_COLOR}


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, default=DEFAULT_INPUT)
    parser.add_argument("--output-dir", type=Path, default=Path(__file__).resolve().parent / "results")
    args = parser.parse_args()
    print(json.dumps(plot_exclusion(args.input, args.output_dir), indent=2))
