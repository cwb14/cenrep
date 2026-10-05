#!/usr/bin/env bash
# One-off helper (2026-09-25): classify changed in cenrep 0.3.2 (strong-family precedence). For each named
# complex, wait until its current run has ended (cenrep.exit written, no cenrep process left); if that run used
# code older than 0.3.2 and exited 0, relaunch it with CENREP_ARGS="--force classify" (earlier checkpoints reused).
# Session-independent: (setsid nohup bash rerun_classify_when_done.sh C1 C2 ... > log 2>&1 < /dev/null &)
set -uo pipefail
CLEAN=/data2/chris/tetraploid/allo10/clean
NEED=0.3.2
pending=("$@")

older() { [[ $(printf '%s\n%s\n' "$1" "$NEED" | sort -V | head -1) == "$1" && $1 != "$NEED" ]]; }
running() { ps -eo args | awk -v f="/clean/$1/centromere/cenrep.time" '$1=="/usr/bin/time" && index($0,f)' | grep -q .; }

while (( ${#pending[@]} )); do
  next=()
  for c in "${pending[@]}"; do
    out=$CLEAN/$c/centromere
    if [[ ! -f $out/cenrep.exit ]] || running "$c"; then next+=("$c"); continue; fi
    ver=$(head -1 "$out/cenrep.log" | grep -o 'cenrep [0-9.]*' | cut -d' ' -f2)
    rc=$(cat "$out/cenrep.exit")
    if [[ $rc != 0 ]]; then
      echo "$(date '+%F %T') $c: run $ver ended with exit $rc; NOT relaunching (needs a look)"
    elif older "$ver"; then
      CENREP_ARGS="--force classify" PARALLEL=1 STAGGER=0 bash "$CLEAN/run_centromere.sh" "$c"
      echo "$(date '+%F %T') $c: run $ver done; relaunched with --force classify"
    else
      echo "$(date '+%F %T') $c: run $ver is current; nothing to do"
    fi
  done
  pending=("${next[@]}")
  (( ${#pending[@]} )) && sleep 60
done
echo "$(date '+%F %T') all handled"
