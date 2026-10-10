#!/usr/bin/env bash
# Generate both polarization and Michel plots from raw MusAirS data.
set -euo pipefail
if [[ $# -lt 2 || $# -gt 5 ]]; then
    echo "Usage: bash $0 PROTON.root HELIUM.root [OUTPUT_DIRECTORY [EVENTS_PER_COMPONENT [SEED]]]" >&2
    exit 2
fi
script_dir="$(cd "$(dirname "$0")" && pwd)"
python="${PYTHON:-python3}"
command -v root >/dev/null || { echo "ROOT is required (root not found)." >&2; exit 1; }
"$python" -c 'import matplotlib, scipy' || {
    echo "Install dependencies from $script_dir/requirements.txt first." >&2
    exit 1
}
output="${3:-pmu_output}"
events="${4:-100000000}"
seed="${5:-20241025}"
# Check all outputs before the expensive raw-data pass. No stage overwrites files.
"$python" - "$output" "$events" "$seed" <<'PY'
from pathlib import Path
import sys
try:
    events, seed = int(sys.argv[2]), int(sys.argv[3])
    if not 1 <= events <= 10000000000 or not 1 <= seed <= 4294967295:
        raise ValueError("EVENTS_PER_COMPONENT must be 1..10000000000; SEED must be 1..4294967295")
    for name in ("polarizationData.root", "polarizationData.csv", "provenance.txt",
                 "combinedPmu_revise.pdf", "michelData.root", "michelData.csv",
                 "michel_provenance.txt", "michelComponents.pdf"):
        path = Path(sys.argv[1]) / name
        if path.exists():
            raise ValueError(f"Refusing to overwrite {path}; choose a new output directory")
except ValueError as exc:
    sys.exit(str(exc))
PY
root_call() {
    local expression
    # JSON literals escape quotes/backslashes in C++ arguments. Numeric Michel
    # arguments are parsed as integers before being interpolated into the macro.
    expression="$("$python" - "$@" <<'PY'
import json
import sys
macro, args = sys.argv[1], sys.argv[2:]
if macro.endswith('generate_michel.C'):
    args[2:] = [int(value) for value in args[2:]]
print(macro + '(' + ','.join(json.dumps(value) for value in args) + ')')
PY
)"
    root -l -b -q "$expression"
}
root_call "$script_dir/analyze_pmu.C" "$1" "$2" "$output"
"$python" "$script_dir/plot_pmu.py" "$output/polarizationData.csv"
root_call "$script_dir/generate_michel.C" "$output/polarizationData.root" "$output" "$events" "$seed"
"$python" "$script_dir/plot_michel.py" "$output/michelData.csv"
