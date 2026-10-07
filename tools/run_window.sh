#!/usr/bin/env bash
# Runs run_bio_night.R inside a time window: waits until START, then stops the
# run at END (SIGINT, as Ctrl-C would; targets keeps every finished target and
# the same command resumes the next night). Used for M2.4, which needs several
# nights (developer, 2026-10-07).
#
#   tools/run_window.sh 22:00 09:00 pauvert korhonen_wood ... species
#
# START already past: the run starts at once. END is the first END after the
# start (09:00 the next morning for a 22:00 start). Launch it in a detached
# tmux session so it survives the end of the desktop session:
#
#   tmux new-session -d -s m24 'tools/run_window.sh 22:00 09:00 <datasets>'
#
# Its own log (wait, start, stop) goes to logs/run_window_<date>.log; the run
# writes logs/run_bio_<date>.log as usual.

set -u

cd "$(dirname "$0")/.." || exit 1
mkdir -p logs

start=$1
end=$2
shift 2

log="logs/run_window_$(date +%Y%m%d_%H%M%S).log"
say() {
  echo "$(date '+%Y-%m-%d %H:%M:%S') | $*" | tee -a "$log"
}

now=$(date +%s)
start_s=$(date -d "$start" +%s)
if ((start_s < now)); then
  start_s=$now
fi
end_s=$(date -d "$end" +%s)
while ((end_s <= start_s)); do
  end_s=$((end_s + 86400))
done

say "waiting until $(date -d "@$start_s" '+%Y-%m-%d %H:%M'), stop at $(date -d "@$end_s" '+%Y-%m-%d %H:%M'); datasets: $*"
sleep $((start_s - $(date +%s)))

remaining=$((end_s - $(date +%s)))
say "start, ${remaining} s allowed; heaviest processes now (another job distorts the timings):"
ps -eo pid,rss,etime,args --sort=-rss | head -6 | cut -c1-150 | tee -a "$log"
timeout --signal=INT "${remaining}s" Rscript run_bio_night.R "$@"
status=$?
if ((status == 124)); then
  say "stopped at the end of the window; rerun the same command to continue"
else
  say "run_bio_night.R exited with status $status"
fi

sleep 10
left=$(pgrep -af "crew_worker|vsearch|blastn|sample_tree" | grep -v pgrep)
if [[ -n "$left" ]]; then
  say "processes left after the stop:"
  echo "$left" | tee -a "$log"
fi
