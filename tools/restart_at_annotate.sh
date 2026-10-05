#!/usr/bin/env bash
# One-off helper (2026-09-25): runs started on cenrep 0.3.0 have a slow annotate stage (1-2 threads).
# For each named complex, wait until its 0.3.0 run logs "annotate: start", stop that run, and relaunch the
# complex with the current code and CENREP_ARGS="--force annotate" (TRF/families checkpoints are reused).
# Session-independent: launch with  setsid nohup bash restart_at_annotate.sh C1 C2 ... > log 2>&1 < /dev/null &
set -uo pipefail
CLEAN=/data2/chris/tetraploid/allo10/clean
pending=("$@")

stop_run() {                      # kill the cenrep python tree under this complex's /usr/bin/time process
  local c=$1 t
  t=$(ps -eo pid,args | awk -v f="/clean/$c/centromere/cenrep.time" '$2=="/usr/bin/time" && index($0,f) {print $1}')
  [[ -n $t ]] || return 1
  for p in $(pgrep -P "$t"); do
    for q in $(pgrep -P "$p"); do pgrep -P "$q" | xargs -r kill; kill "$q"; done
    kill "$p"
  done
}

while (( ${#pending[@]} )); do
  next=()
  for c in "${pending[@]}"; do
    log=$CLEAN/$c/centromere/cenrep.log
    if ! head -1 "$log" | grep -q "cenrep 0.3.0"; then
      echo "$(date '+%F %T') $c: not a 0.3.0 run; leaving it alone"; continue
    fi
    if [[ -f $CLEAN/$c/centromere/cenrep.exit ]]; then
      echo "$(date '+%F %T') $c: run ended (exit $(cat "$CLEAN/$c/centromere/cenrep.exit")) before annotate; not restarting"; continue
    fi
    if grep -q "INFO annotate: start" "$log"; then
      echo "$(date '+%F %T') $c: reached annotate on 0.3.0; stopping"
      stop_run "$c" || echo "$(date '+%F %T') $c: no running process found"
      for _ in $(seq 60); do [[ -f $CLEAN/$c/centromere/cenrep.exit ]] && break; sleep 2; done
      CENREP_ARGS="--force annotate" PARALLEL=1 STAGGER=0 bash "$CLEAN/run_centromere.sh" "$c"
      echo "$(date '+%F %T') $c: relaunched with --force annotate"
      continue
    fi
    next+=("$c")
  done
  pending=("${next[@]}")
  (( ${#pending[@]} )) && sleep 30
done
echo "$(date '+%F %T') all restarts handled"
