#!/usr/bin/env bash
# Generate polarization and Michel plots with explicit raw event-weight handling.
set -euo pipefail
script_dir="$(cd "$(dirname "$0")" && pwd)"
python="${PYTHON:-python3}"
# shlex.quote makes the assignments safe for spaces, quotes and shell characters.
# Keep assignment separate from eval so argparse failures propagate under set -e.
configuration="$("$python" - "$@" <<'PY'
import argparse
import math
import shlex

parser = argparse.ArgumentParser(
    prog="run_pipeline.sh", description="Generate polarization and independent Michel plots from proton/helium MusAirS ROOT files.")
parser.add_argument("proton", help="proton raw ROOT file")
parser.add_argument("helium", help="helium raw ROOT file")
parser.add_argument("output", nargs="?", default="pmu_output")
parser.add_argument("events", nargs="?", type=int, default=100000000, help="events per Michel source (default: 100M)")
parser.add_argument("seed", nargs="?", type=int, default=20241025)
weights = parser.add_mutually_exclusive_group()
weights.add_argument("--use-weights", dest="use_weights", action="store_true", help="require and apply raw Weight in all input muon trees")
weights.add_argument("--no-weights", dest="use_weights", action="store_false", help="ignore raw Weight, use unit weights (default)")
parser.set_defaults(use_weights=False)
parser.add_argument("--all-muons", action="store_true", help="keep all selected muons instead of historical parent/event deduplication")
parser.add_argument("--proton-scale", type=float, default=9.704844263299673)
parser.add_argument("--helium-scale", type=float, default=100.0)
args = parser.parse_intermixed_args()
if not 1 <= args.events <= 10000000000 or not 1 <= args.seed <= 4294967295:
    parser.error("EVENTS_PER_COMPONENT must be 1..10000000000; SEED must be 1..4294967295")
scales = (args.proton_scale, args.helium_scale)
if any(not math.isfinite(value) or value < 0 for value in scales) or not any(value > 0 for value in scales):
    parser.error("Source scales must be finite, nonnegative and not both zero")
if not args.output:
    parser.error("Output directory must not be empty")
for key, value in dict(proton=args.proton, helium=args.helium, output=args.output,
                       events=args.events, seed=args.seed,
                       use_weights=str(args.use_weights).lower(),
                       legacy_deduplicate=str(not args.all_muons).lower(),
                       proton_scale=args.proton_scale, helium_scale=args.helium_scale).items():
    print(key + "=" + shlex.quote(str(value)))
PY
)"
# --help emits usage rather than assignments; return without evaluating it.
if [[ "$configuration" == usage:* ]]; then
    printf '%s\n' "$configuration"
    exit 0
fi
eval "$configuration"
command -v root >/dev/null || { echo "ROOT is required (root not found)." >&2; exit 1; }
"$python" -c 'import matplotlib, scipy' || {
    echo "Install dependencies from $script_dir/requirements.txt first." >&2
    exit 1
}
# Check all outputs before the expensive raw-data pass. No stage overwrites files.
"$python" - "$output" <<'PY'
from pathlib import Path
import sys
for name in ("polarizationData.root", "polarizationData.csv", "provenance.txt",
             "combinedPmu_revise.pdf", "michelData.root", "michelData.csv",
             "michel_provenance.txt", "michelComponents.pdf"):
    path = Path(sys.argv[1]) / name
    if path.exists():
        sys.exit(f"Refusing to overwrite {path}; choose a new output directory")
PY
root_call() {
    local expression
    # JSON literals escape quotes/backslashes in C++ strings; numeric and bool
    # parameters are converted separately rather than passed as string literals.
    expression="$("$python" - "$@" <<'PY'
import json
import sys
macro, args = sys.argv[1], sys.argv[2:]
if macro.endswith('generate_michel.C'):
    args[2:] = [int(value) for value in args[2:]]
elif macro.endswith('analyze_pmu.C'):
    args[3] = int(args[3])
    args[4:7] = [value == 'true' for value in args[4:7]]
    args[7:] = [float(value) for value in args[7:]]
print(macro + '(' + ','.join(json.dumps(value) for value in args) + ')')
PY
)"
    root -l -b -q "$expression"
}
root_call "$script_dir/analyze_pmu.C" "$proton" "$helium" "$output" 50 "$legacy_deduplicate" false "$use_weights" "$proton_scale" "$helium_scale"
"$python" "$script_dir/plot_pmu.py" "$output/polarizationData.csv"
root_call "$script_dir/generate_michel.C" "$output/polarizationData.root" "$output" "$events" "$seed"
"$python" "$script_dir/plot_michel.py" "$output/michelData.csv"
