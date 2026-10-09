#!/usr/bin/env bash
# Run the two stages from any working directory. Dependencies must be installed.
set -euo pipefail
if [[ $# -lt 2 || $# -gt 3 ]]; then
    echo "Usage: bash $0 PROTON.root HELIUM.root [OUTPUT_DIRECTORY]" >&2
    exit 2
fi
script_dir="$(cd "$(dirname "$0")" && pwd)"
python="${PYTHON:-python3}"
command -v root >/dev/null || { echo "ROOT is required (root not found)." >&2; exit 1; }
"$python" -c 'import matplotlib' || {
    echo "Install plotting dependencies from $script_dir/requirements.txt first." >&2
    exit 1
}
output="${3:-pmu_output}"
# JSON string literals safely escape quotes/backslashes in the ROOT macro arguments.
expression="$("$python" - "$script_dir/analyze_pmu.C" "$1" "$2" "$output" <<'PY'
import json
import sys
print(sys.argv[1] + '(' + ','.join(json.dumps(value) for value in sys.argv[2:]) + ')')
PY
)"
root -l -b -q "$expression"
"$python" "$script_dir/plot_pmu.py" "$output/polarizationData.csv"
