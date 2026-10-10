#!/usr/bin/env python3
"""Compare direct and weighted marginal Michel sampling with the exact expectation."""
import argparse
from pathlib import Path
import sys

from plot_michel import CURVES
from plot_pmu import read_data


def event_count(csv_path):
    metadata = csv_path.with_name("michel_provenance.txt").read_text()
    for line in metadata.splitlines():
        if line.startswith("events_per_component="):
            value = int(line.split("=", 1)[1])
            if value > 0:
                return value
    raise ValueError(f"Missing events_per_component in provenance beside {csv_path}")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("direct", type=Path)
    parser.add_argument("reference", type=Path)
    parser.add_argument("polarization", type=Path)
    parser.add_argument("-o", "--output", type=Path, required=True)
    args = parser.parse_args()
    try:
        if args.output.exists():
            raise ValueError(f"Refusing to overwrite {args.output}")
        inputs = (args.direct, args.reference, args.polarization)
        if any(args.output.resolve() == path.resolve() for path in inputs):
            raise ValueError("Output must not replace an input")
        names = [item[0] for item in CURVES]
        datasets = [read_data(path, names=names, normalized=False) for path in inputs[:2]]
        events = [event_count(path) for path in inputs[:2]]
        polarization = read_data(args.polarization, names=[name + "Pmu" for name in names])
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        import numpy as np

        figure, axes = plt.subplots(2, 3, figsize=(12, 6), sharex=True,
                                    gridspec_kw={"height_ratios": [2, 1]})
        for column, (name, _, _, label) in enumerate(CURVES):
            source = polarization[name + "Pmu"]
            mean = sum(row[1] * row[2] for row in source) / sum(row[2] for row in source)
            rows = datasets[0][name]
            x = np.array([row[1] for row in rows])
            # These generators both use 100 equal bins on [-1,1].
            if len(rows) != 100 or not np.allclose(x, np.linspace(-0.99, 0.99, 100), rtol=0, atol=1e-12):
                raise ValueError("Comparison expects 100 equal cos(theta) bins on [-1,1]")
            expected = (1 + mean * x / 3) / 2
            for data, count, color, style, method in zip(
                    datasets, events, ("#0072B2", "#D55E00"), ("-", "--"),
                    ("Direct, unweighted", "Original, weighted")):
                other_x = np.array([row[1] for row in data[name]])
                if other_x.shape != x.shape or not np.allclose(other_x, x, rtol=0, atol=1e-12):
                    raise ValueError("Histogram binning must match")
                y = np.array([row[2] for row in data[name]]) / (count * 0.02)
                error = np.array([row[3] for row in data[name]]) / (count * 0.02)
                axes[0, column].plot(x, y, color=color, linestyle=style,
                                     linewidth=1.3, label=f"{method} ({count:.0e}/source)")
                residual = np.divide(y - expected, error, out=np.zeros_like(y), where=error > 0)
                axes[1, column].plot(x, residual, color=color, linestyle=style, linewidth=0.9)
            axes[0, column].plot(x, expected, color="black", linestyle=":", linewidth=1,
                                 label="Exact marginal expectation")
            axes[0, column].set_title(label)
            axes[0, column].legend(fontsize=8)
            axes[0, column].grid(alpha=0.25)
            axes[1, column].axhline(0, color="black", linewidth=0.7)
            axes[1, column].grid(alpha=0.25)
            axes[1, column].set_xlabel(r"$\cos\theta$")
        axes[0, 0].set_ylabel(r"$(1/N)\,dN/d\cos\theta$")
        axes[1, 0].set_ylabel("(MC - exact) / bin error")
        figure.suptitle("Independent Michel sources: marginal shape comparison")
        figure.tight_layout()
        args.output.parent.mkdir(parents=True, exist_ok=True)
        figure.savefig(args.output, bbox_inches="tight")
        plt.close(figure)
        print(f"Saved {args.output}")
    except (OSError, ValueError, ImportError) as error:
        print(f"Comparison failed: {error}", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
