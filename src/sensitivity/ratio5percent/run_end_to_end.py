"""Raw MusAirS -> Michel templates -> ratio measurement and pion exclusion plots."""

import argparse
from datetime import datetime, timezone
import json
from pathlib import Path
import shutil
import subprocess

import numpy as np
import uproot

from analysis import generate_distributions, load_templates, scan_distributions, scan_exclusion
from plot_exclusion_sensitivity import plot_exclusion
from plot_ratio import plot_scan

HERE = Path(__file__).resolve().parent
UPSTREAM = HERE.parents[1] / "analysis_pmu_from_MusAirS"
CHARGES = {"muon": ("piMuon", "KMuon"), "antimuon": ("piAntiMuon", "KAntiMuon")}


def call_root(macro, arguments, log_path):
    expression = f"{macro}({','.join(json.dumps(arg) for arg in arguments)})"
    with log_path.open("w") as log:
        result = subprocess.run(["root", "-l", "-b", "-q", expression],
                                stdout=log, stderr=subprocess.STDOUT)
    if result.returncode:
        raise RuntimeError(f"ROOT failed ({result.returncode}); see {log_path}")
    return expression


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    mode = parser.add_mutually_exclusive_group(required=True)
    mode.add_argument("--proton", type=Path, help="Raw proton MusAirS ROOT file")
    mode.add_argument("--templates", type=Path, help="Existing Michel template interface ROOT file")
    parser.add_argument("--helium", type=Path, help="Raw helium MusAirS ROOT file; required with --proton")
    parser.add_argument("--output-dir", type=Path, default=HERE / "results/end_to_end")
    parser.add_argument("--raw-weights", action="store_true", help="Read raw Weight branches")
    parser.add_argument("--legacy-deduplicate", action=argparse.BooleanOptionalAction, default=True)
    parser.add_argument("--proton-scale", type=float, default=9.704844263299673)
    parser.add_argument("--helium-scale", type=float, default=100.0)
    parser.add_argument("--polarization-bins", type=int, default=50)
    parser.add_argument("--template-events", type=int, default=100_000_000)
    parser.add_argument("--seed", type=int, default=20241025)
    parser.add_argument("--events", type=float, default=10_000_000, help="Experimental events for measurement")
    parser.add_argument("--true-ratio", type=float, default=0.05)
    parser.add_argument("--scan-start", type=float, default=0.01)
    parser.add_argument("--scan-end", type=float, default=0.10)
    parser.add_argument("--points", type=int, default=1000)
    parser.add_argument("--ratio-axis", choices=("physical", "notebook"), default="physical")
    parser.add_argument("--measurement-charge", choices=CHARGES, default="muon")
    parser.add_argument("--exclusion-charge", choices=CHARGES, default="antimuon")
    args = parser.parse_args()
    if bool(args.proton) != bool(args.helium):
        parser.error("Supply both --proton and --helium, or only --templates")
    inputs = [args.templates] if args.templates else [args.proton, args.helium]
    for path in inputs:
        if not path.is_file():
            parser.error(f"Input not found: {path}")
    if not args.templates and not shutil.which("root"):
        parser.error("ROOT executable is required for raw-data mode")
    if (not 1 <= args.template_events <= 10_000_000_000 or not 1 <= args.seed <= 4294967295
            or args.polarization_bins < 1):
        parser.error("Invalid template event count, seed or polarization bin count")
    if (not np.isfinite(args.events) or args.events <= 0 or args.points < 4
            or not 0 <= args.scan_start < args.scan_end
            or not args.scan_start <= args.true_ratio <= args.scan_end):
        parser.error("Invalid experimental event count, scan range or true ratio")
    scales = [args.proton_scale, args.helium_scale]
    if not np.all(np.isfinite(scales)) or min(scales) < 0 or sum(scales) <= 0:
        parser.error("Source scales must be finite, nonnegative and not both zero")
    directory = args.output_dir.resolve()
    if directory.exists() and (not directory.is_dir() or any(directory.iterdir())):
        parser.error(f"Output directory must be empty; refusing to overwrite: {directory}")
    measurement_names = CHARGES[args.measurement_charge]
    exclusion_names = CHARGES[args.exclusion_charge]
    if args.templates:
        # Validate both downstream interfaces before creating any output.
        for names in (measurement_names, exclusion_names):
            load_templates(args.templates, *names)
    directory.mkdir(parents=True, exist_ok=True)
    manifest = {"started_utc": datetime.now(timezone.utc).isoformat(),
                "status": "running", "parameters": {
                    key: str(value.resolve()) if isinstance(value, Path) else value
                    for key, value in vars(args).items()},
                "template_interface": {"measurement": measurement_names, "exclusion": exclusion_names,
                                       "variable": "cos(theta)", "range": [-1, 1]},
                "threshold_convention": "measurement delta WW chi2=1; exclusion WW chi2=9; not coverage-calibrated",
                "root_commands": []}
    manifest_path = directory / "manifest.json"

    def save_manifest():
        manifest_path.write_text(json.dumps(manifest, indent=2) + "\n")

    save_manifest()
    try:
        if args.templates:
            templates = args.templates.resolve()
        else:
            upstream_dir = directory / "upstream"
            expression = call_root(UPSTREAM / "analyze_pmu.C", [
                str(args.proton.resolve()), str(args.helium.resolve()), str(upstream_dir),
                args.polarization_bins, args.legacy_deduplicate, False, args.raw_weights,
                args.proton_scale, args.helium_scale,
            ], directory / "analyze_pmu.log")
            manifest["root_commands"].append(expression)
            expression = call_root(UPSTREAM / "generate_michel.C", [
                str(upstream_dir / "polarizationData.root"), str(upstream_dir),
                args.template_events, args.seed, False, True,
            ], directory / "generate_michel.log")
            manifest["root_commands"].append(expression)
            templates = upstream_dir / "michelData.root"
        manifest["templates"] = str(templates)
        with uproot.open(templates) as source:
            if "provenance" in source:
                record = source["provenance"]
                manifest["template_provenance"] = (
                    record.member("fTitle") if record.classname == "TNamed" else str(record)
                )
        for names in (measurement_names, exclusion_names):
            load_templates(templates, *names)
        measurement = directory / "measurement"
        distributions = measurement / "michelDis.root"
        scan = measurement / "measure_ratio_chi2test.root"
        generate_distributions(templates, distributions, args.true_ratio, args.scan_start,
                               args.scan_end, args.points, *measurement_names)
        scan_distributions(distributions, scan, args.events, args.ratio_axis)
        manifest["measurement_summary"] = plot_scan(scan, measurement)
        exclusion_dir = directory / "exclusion"
        exclusion_result = exclusion_dir / "eigen_sensitiviy_output.root"
        scan_exclusion(templates, exclusion_result, pion_name=exclusion_names[0],
                       kaon_name=exclusion_names[1])
        manifest["exclusion_summary"] = plot_exclusion(exclusion_result, exclusion_dir)
        manifest["status"] = "complete"
        save_manifest()
        print(f"Complete: {directory}\nMeasurement: {measurement}\nExclusion: {exclusion_dir}")
    except Exception as error:
        manifest["status"] = "failed"
        manifest["error"] = str(error)
        save_manifest()
        raise


if __name__ == "__main__":
    main()
