# H36, step 3: take one idle card of this node for this job. The cluster sets no CUDA_VISIBLE_DEVICES and schedules no cards, so a job
# takes one under a node-local lock (one lock file per card, held for as long as the job runs) and only if no process of any user is on
# it (no memory in use, twice, ten seconds apart: a card that looks idle may be a job still loading). Waits up to two hours, then gives up.
# Sourced by the job files; sets CUDA_VISIBLE_DEVICES and tells the dashboard's comment which card.
take_card() {
  local n
  n=$(nvidia-smi --query-gpu=index --format=csv,noheader | wc -l)
  for attempt in $(seq 1 24); do
    for i in $(seq 0 $((n - 1))); do
      exec {H36_LOCK}>"/tmp/h36-saprot-card$i.lock"
      if flock -n "$H36_LOCK"; then
        local a b
        a=$(nvidia-smi --query-gpu=memory.used --format=csv,noheader,nounits -i "$i" | tr -d ' ')
        sleep 10
        b=$(nvidia-smi --query-gpu=memory.used --format=csv,noheader,nounits -i "$i" | tr -d ' ')
        if [ "$a" -lt 500 ] && [ "$b" -lt 500 ]; then
          export CUDA_VISIBLE_DEVICES=$i
          echo "taking card $i ($a and $b MiB in use)"
          scontrol update JobId="$SLURM_JOB_ID" Comment="${H36_COMMENT_PREFIX:-seat=mitoage/lead; what=H36 step 3 SaProt}; eta=${H36_ETA:-?}; cuda=$i" 2>/dev/null || true
          return 0
        fi
        flock -u "$H36_LOCK"
      fi
      exec {H36_LOCK}>&-
    done
    echo "no idle card on $(hostname) (attempt $attempt), waiting"
    sleep 300
  done
  echo "no idle card on $(hostname) after two hours, stopping"
  exit 1
}
