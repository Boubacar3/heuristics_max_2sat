#!/usr/bin/env bash

# script_linear_relaxation.sh
# Solve the linear relaxation of every .lp file in a directory and save a CSV of results.
# Usage: ./script_linear_relaxation.sh <lp-folder> <output-csv>

set -euo pipefail

if [[ $# -ne 2 ]]; then
  cat <<EOF
Usage: $0 <lp-folder> <output-csv>

Example:
  $0 grid/ results_linear_relaxation.csv
EOF
  exit 1
fi

LP_FOLDER=$1
OUTPUT_CSV=$2
SOLVER_CMD=${CPLEX_CMD:-"$HOME/cplex"}

if [[ ! -d "$LP_FOLDER" ]]; then
  echo "Error: folder '$LP_FOLDER' does not exist." >&2
  exit 2
fi

if ! command -v "$SOLVER_CMD" >/dev/null 2>&1; then
  echo "Error: solver command '$SOLVER_CMD' not found. Set CPLEX_CMD or install cplex." >&2
  exit 3
fi

mkdir -p "$(dirname "$OUTPUT_CSV")"

echo "file,objective,status,solve_seconds" > "$OUTPUT_CSV"

for lp in "$LP_FOLDER"/*.lp; do
  if [[ ! -f "$lp" ]]; then
    continue
  fi

  lp_name=$(basename "$lp")
  echo "Solving linear relaxation for: $lp_name"

  tmp_output=$(mktemp)
  start_time=$(date +%s.%N)

  # Run CPLEX in batch mode and capture the relevant output.
  "$SOLVER_CMD" <<EOF >"$tmp_output" 2>&1
read "$lp"
change problem lp
optimize
display solution objective
display solution status
quit
EOF

  end_time=$(date +%s.%N)
  solve_seconds=$(awk "BEGIN { print $end_time - $start_time }" )

  objective=$(grep -Ei 'Objective\s*=|Objective\s+value|Objective value' "$tmp_output" | head -n 1 | sed -E 's/.*([0-9]+(\.[0-9]+)?).*/\1/')
  if [[ -z "$objective" ]]; then
    objective=NA
  fi

  status=$(grep -Ei 'Optimal|Infeasible|Unbounded|Status|Solution status' "$tmp_output" | head -n 1 | sed -E 's/^\s*//; s/\s*$//' | tr ',' ';')
  if [[ -z "$status" ]]; then
    status=UNKNOWN
  fi

  printf '%s,%s,%s,%s\n' "$lp_name" "$objective" "$status" "$solve_seconds" >> "$OUTPUT_CSV"
  rm -f "$tmp_output"
done

echo "Results saved to $OUTPUT_CSV"
