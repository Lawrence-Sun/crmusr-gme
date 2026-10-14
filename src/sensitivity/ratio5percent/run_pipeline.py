"""Generate Asimov distributions, scan K/pi, and reproduce the notebook figure."""

import argparse
import json
from pathlib import Path

from analysis import generate_distributions, scan_distributions
from plot_ratio import plot_scan

REPO_ROOT = Path(__file__).resolve().parents[3]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--templates", type=Path, default=REPO_ROOT / "dev/megaMichelDis.root")
    parser.add_argument("--pion-name", default="piMichelDis")
    parser.add_argument("--kaon-name", default="kaMichelDis")
    parser.add_argument("--output-dir", type=Path, default=Path(__file__).resolve().parent / "results")
    parser.add_argument("--events", type=float, default=10_000_000)
    parser.add_argument("--true-ratio", type=float, default=0.05)
    parser.add_argument("--scan-start", type=float, default=0.01)
    parser.add_argument("--scan-end", type=float, default=0.10)
    parser.add_argument("--points", type=int, default=1000)
    parser.add_argument("--ratio-axis", choices=("notebook", "physical"), default="notebook")
    args = parser.parse_args()
    if not args.templates.is_file():
        parser.error(f"Template file not found: {args.templates}; supply --templates")
    if not args.events > 0 or args.events == float("inf"):
        parser.error("--events must be finite and positive")
    distribution_path = args.output_dir / "michelDis.root"
    scan_path = args.output_dir / "measure_ratio_chi2test.root"
    if args.templates.resolve() in (distribution_path.resolve(), scan_path.resolve()):
        parser.error("Output paths must not overwrite the input templates")
    generate_distributions(args.templates, distribution_path, args.true_ratio,
                           args.scan_start, args.scan_end, args.points,
                           args.pion_name, args.kaon_name)
    scan_distributions(distribution_path, scan_path, args.events, args.ratio_axis)
    summary = plot_scan(scan_path, args.output_dir)
    print(json.dumps(summary, indent=2))
    print(f"Outputs: {args.output_dir.resolve()}")


if __name__ == "__main__":
    main()
