#!/bin/bash
PYMOL="/Applications/PyMOL.app/Contents/bin/pymol"
DIR="$(cd "$(dirname "$0")" && pwd)"
mkdir -p "$DIR/output/logs"
TS=$(date +"%Y%m%d_%H%M%S")
LOG="$DIR/output/logs/batch_${TS}.log"
echo "=== GLINT Batch Run ===" | tee "$LOG"
echo "Started: $(date)" | tee -a "$LOG"
for f in run_case1_5FQD.py run_case2_6UAN.py run_case3_5HXB.py run_case4_6BOY.py; do
  echo "" | tee -a "$LOG"
  echo ">>> $f" | tee -a "$LOG"
  S=$(date +%s)
  "$PYMOL" -cq "$DIR/$f" >> "$LOG" 2>&1
  E=$?
  D=$(( $(date +%s) - S ))
  [ $E -eq 0 ] && echo "  ✅ OK (${D}s)" | tee -a "$LOG" || echo "  ❌ FAIL exit=$E (${D}s)" | tee -a "$LOG"
done
echo "" | tee -a "$LOG"
echo "=== Done: $(date) ===" | tee -a "$LOG"
echo "Log: $LOG"
