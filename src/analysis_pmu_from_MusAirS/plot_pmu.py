#!/usr/bin/env python3
"""Plot the portable CSV produced by analyze_pmu.C; no PyROOT needed."""
import argparse
import csv
import math
from pathlib import Path
import sys

CURVES = (
    ("piMuonPmu", "r", "o", r"$\mu^-(\pi^-)$"),
    ("KMuonPmu", "g", ">", r"$\mu^-(K^-)$"),
    ("KLongMuonPmu", "b", "<", r"$\mu^-(K^0_L)$"),
)


def read_data(path):
    series = {name: [] for name, *_ in CURVES}
    with path.open(newline="", encoding="utf-8") as stream:
        reader = csv.DictReader(stream)
        required = {"histogram", "bin", "x", "y", "yerr"}
        if not required.issubset(reader.fieldnames or []):
            raise ValueError(f"CSV must contain columns: {', '.join(sorted(required))}")
        for line, row in enumerate(reader, 2):
            name = row["histogram"]
            if name not in series:
                continue
            try:
                index = int(row["bin"])
                x, y, error = (float(row[key]) for key in ("x", "y", "yerr"))
            except (ValueError, TypeError) as exc:
                raise ValueError(f"Invalid numeric data on CSV line {line}") from exc
            if not all(math.isfinite(value) for value in (x, y, error)):
                raise ValueError(f"Non-finite data on CSV line {line}")
            if not -1 <= x <= 1 or y < 0 or error < 0:
                raise ValueError(f"Invalid range on CSV line {line}")
            series[name].append((index, x, y, error))
    reference_x = None
    for name, rows in series.items():
        rows.sort()
        if not rows or [row[0] for row in rows] != list(range(1, len(rows) + 1)):
            raise ValueError(f"Missing or duplicate bins for {name}")
        x_values = [row[1] for row in rows]
        if any(left >= right for left, right in zip(x_values, x_values[1:])):
            raise ValueError(f"Bin centers are not increasing for {name}")
        if reference_x is not None and x_values != reference_x:
            raise ValueError("Histogram bin centers do not match")
        reference_x = x_values
        if not math.isclose(sum(row[2] for row in rows), 1, rel_tol=1e-8, abs_tol=1e-8):
            raise ValueError(f"{name} is empty or not normalized to 1")
    return series


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("data", type=Path, help="polarizationData.csv from analyze_pmu.C")
    parser.add_argument("-o", "--output", type=Path, help="default: beside CSV, combinedPmu_revise.pdf")
    parser.add_argument("--overwrite", action="store_true", help="allow replacing an existing figure")
    parser.add_argument("--show", action="store_true", help="also open an interactive figure")
    args = parser.parse_args()
    output = args.output or args.data.with_name("combinedPmu_revise.pdf")
    try:
        if output.resolve() == args.data.resolve():
            raise ValueError("Output must not replace the input CSV")
        if output.exists() and not args.overwrite:
            raise ValueError(f"Refusing to overwrite {output}; use --overwrite")
        series = read_data(args.data)
        import matplotlib
        if not args.show:
            matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        from matplotlib import font_manager

        # Use the original legend font only if it is installed; avoid noisy fallbacks.
        font = "Times New Roman" if any(
            item.name == "Times New Roman" for item in font_manager.fontManager.ttflist
        ) else "DejaVu Serif"
        figure, axes = plt.subplots(figsize=(5, 5))
        for name, color, marker, label in CURVES:
            # Zero bins cannot be shown on a log axis. Do not add artificial counts.
            positive = [row for row in series[name] if row[2] > 0]
            axes.errorbar(
                [row[1] for row in positive], [row[2] for row in positive],
                yerr=[row[3] for row in positive], fmt=marker, label=label,
                capsize=2, color=color,
            )
        axes.tick_params(labelsize=12)
        axes.set_xlabel(r"$\mathcal{P}_\mu$", fontsize=16)
        axes.set_ylabel("counts", fontsize=16)  # Match the historical figure.
        axes.set_ylim(1e-4, 1)
        axes.set_yscale("log")
        axes.legend(prop={"family": font, "size": 12})
        output.parent.mkdir(parents=True, exist_ok=True)
        figure.savefig(output, bbox_inches="tight")
        if args.show:
            plt.show()
        plt.close(figure)
        print(f"Saved {output}")
    except (OSError, ValueError, ImportError) as exc:
        print(f"PMU plotting failed: {exc}", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
