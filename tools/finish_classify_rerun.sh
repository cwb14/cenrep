#!/usr/bin/env bash
# One-off (2026-09-26): wait for the v0.3.6 launcher (pid 2535367) to exit, then rerun classify onward on 0.3.7 for
# the complexes it was still running. Full 0.3.6 logs are kept as work/cenrep_v0.3.6_full.{log,time}.
# Run detached: setsid nohup bash cenrep/tools/finish_classify_rerun.sh > centromere_finish_classify.log 2>&1 < /dev/null &
set -u
CLEAN=/data2/chris/tetraploid/allo10/clean
PID=${PID:-2535367}
KEEP=/tmp/claude-1007/-data2-chris-tetraploid-allo10-clean/17f72ccc-db98-495e-ac58-1f3555b7e5bc/scratchpad/v036
cd "$CLEAN"
echo "$(date '+%F %T') waiting for pid $PID"
while kill -0 "$PID" 2>/dev/null; do sleep 60; done
ok=()
for c in Carabica Cbursapastoris Msookensis; do
  o=$c/centromere
  if [[ $(cat "$o/cenrep.exit" 2>/dev/null) == 0 ]] && head -1 "$o/cenrep.log" | grep -q 'cenrep 0.3.6'; then
    cp -p "$o/cenrep.log" "$o/work/cenrep_v0.3.6_full.log"; cp -p "$o/cenrep.time" "$o/work/cenrep_v0.3.6_full.time"
    if [[ -d ${KEEP%/*} ]]; then mkdir -p "$KEEP/$c"; cp -p "$o"/*_centromere_summary.tsv "$o/satellite_families.tsv" "$KEEP/$c/"; fi
    ok+=("$c")
  else
    echo "$(date '+%F %T') $c: 0.3.6 run not clean (exit $(cat "$o/cenrep.exit" 2>/dev/null)); not relaunched"
  fi
done
if (( ${#ok[@]} )); then
  echo "$(date '+%F %T') classify rerun: ${ok[*]}"
  PARALLEL=3 CENREP_ARGS="--force classify" bash run_centromere.sh "${ok[@]}"
fi
