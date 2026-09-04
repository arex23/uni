#!/usr/bin/env bash
# Run the per-sample pipeline across the frozen 16-sample analysis cohort (D6).
#
# One Rscript process per sample per stage, on purpose: the SpaNorm adjustment
# peaks near 12.5 GB on a 15 GB machine, and only process exit reliably returns
# that memory. An in-process loop accumulates and gets OOM-killed partway through.
#
# Usage:
#   ./run_cohort.sh                 # entropy, then stemness, over the whole cohort
#   ./run_cohort.sh entropy         # entropy stage only
#   ./run_cohort.sh stemness        # stemness stage only (needs the .rds files)
#   ./run_cohort.sh entropy sample4 sample21   # named samples only
#
# Logs land in results/cohort_qc/logs/<stage>_<sample>.log. A failing sample is
# reported and skipped so one bad sample does not abort the run; the exit status
# is non-zero if anything failed.
#
# Activate the conda environment first (`conda activate stemness`) -- this runs
# whatever Rscript is on PATH and does not manage the environment itself.

set -uo pipefail
cd "$(dirname "$0")"

STAGE="${1:-all}"
if [[ "$STAGE" == "entropy" || "$STAGE" == "stemness" || "$STAGE" == "all" ]]; then
  shift || true
else
  STAGE="all"
fi

if [[ $# -gt 0 ]]; then
  SAMPLES=("$@")
else
  mapfile -t SAMPLES < <(Rscript -e 'source("R/cohort.R"); cat(cohort_samples(), sep="\n")')
fi

if [[ ${#SAMPLES[@]} -eq 0 ]]; then
  echo "No samples resolved. Is data/ populated?" >&2
  exit 1
fi

LOG_DIR="results/cohort_qc/logs"
mkdir -p "$LOG_DIR"

FAILED=()

run_stage() {
  local script="$1" stage="$2" sample="$3"
  local log="$LOG_DIR/${stage}_${sample}.log"
  printf '[%s] %-8s %-10s ... ' "$(date +%H:%M:%S)" "$stage" "$sample"
  if Rscript "$script" "$sample" >"$log" 2>&1; then
    echo "ok"
  else
    echo "FAILED (see $log)"
    FAILED+=("$stage/$sample")
  fi
}

echo "Cohort: ${#SAMPLES[@]} samples — ${SAMPLES[*]}"
echo "Stage:  $STAGE"
echo

for sample in "${SAMPLES[@]}"; do
  if [[ "$STAGE" == "entropy" || "$STAGE" == "all" ]]; then
    run_stage analyze_entropy.R entropy "$sample"
  fi
  if [[ "$STAGE" == "stemness" || "$STAGE" == "all" ]]; then
    run_stage analyze_stemness.R stemness "$sample"
  fi
done

echo
if [[ ${#FAILED[@]} -gt 0 ]]; then
  echo "Failed: ${FAILED[*]}"
  exit 1
fi

echo "All samples completed. Next: Rscript check_cohort_retention.R"
