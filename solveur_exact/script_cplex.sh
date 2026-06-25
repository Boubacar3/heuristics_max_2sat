
#!/usr/bin/env bash
set -uo pipefail

# Batch runner for CPLEX on all .lp files.
# Usage:
#   ./script_cplex.sh [CPLEX_BIN] [INPUT_DIR] [OUT_DIR] [TIMEOUT_SECONDS] [PRM_FILE]
# Defaults:
#   CPLEX_BIN=${1:-$HOME/cplex}
#   INPUT_DIR=${2:-.}
#   OUT_DIR=${3:-results_cplex}
#   TIMEOUT_SECONDS: optional fourth arg, overrides TIMEOUT env var if present
#   PRM_FILE: optional fifth arg, if provided loads the .prm parameter file

CPLEX_BIN=${1:-"$HOME/cplex"}
INPUT_DIR=${2:-.}
OUT_DIR=${3:-results_cplex}
TIMEOUT_ARG=${4:-}
PRM_FILE=${5:-}

# If a timeout argument was provided, use it; otherwise allow TIMEOUT env var to remain.
if [ -n "$TIMEOUT_ARG" ]; then
  TIMEOUT="$TIMEOUT_ARG"
fi

if [ -n "$PRM_FILE" ] && [ ! -f "$PRM_FILE" ]; then
  echo "Error: parameter file '$PRM_FILE' does not exist." >&2
  exit 1
fi

if [ ! -x "$CPLEX_BIN" ] && [ ! -f "$CPLEX_BIN" ]; then
  echo "Warning: CPLEX binary '$CPLEX_BIN' not found or not executable." >&2
  echo "You can pass the path to CPLEX as the first argument." >&2
fi

if [ ! -d "$INPUT_DIR" ]; then
  echo "Error: input directory '$INPUT_DIR' does not exist." >&2
  exit 1
fi

mkdir -p "$OUT_DIR"

# Use timeout if available and TIMEOUT is set
TIMEOUT_CMD=""
if command -v timeout >/dev/null 2>&1 && [ -n "${TIMEOUT:-}" ]; then
	TIMEOUT_CMD="timeout ${TIMEOUT}"
fi

find "$INPUT_DIR" -type f -name '*.lp' -print0 | sort -z | while IFS= read -r -d '' lp; do
	base=$(basename "$lp")
	out="$OUT_DIR/${base%.lp}.cplex.txt"
	echo "Processing: $lp -> $out"

	{
		echo "### CPLEX run for: $lp"
		echo "Started: $(date --iso-8601=seconds 2>/dev/null || date)"
		if [ -n "$PRM_FILE" ]; then
			echo "Parameter file: $PRM_FILE"
		fi
		printf "read %s\n%s\nset threads 1\noptimize\ndisplay solution variables -\nquit\n" "$lp" "${PRM_FILE:+read $PRM_FILE}"
	} | $TIMEOUT_CMD "$CPLEX_BIN" > "$out" 2>&1 || {
		echo "CPLEX exited with non-zero status for $lp (see $out)" >&2
		continue
	}

	echo "Finished: $(date --iso-8601=seconds 2>/dev/null || date)" >> "$out"
done

echo "All done. Results are in: $OUT_DIR"

