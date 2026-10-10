#!/bin/zsh
# Passive-scalar conservation tests: build, run every test of the matrix, summarise.
#
#   run_all.sh <tree> <label> <outdir> [--quick] [--nobuild]
#              [--tests "advect sod amr tva source restart sink"]
#
#   <tree>    RAMSES source tree to build (any tree: the current scheme shows the errors,
#             a tree with the fix passes)
#   <label>   name of the executable family (e.g. base, cma); with label 'cmaoff' the
#             executables of label 'cma' are run with passive_cma=.false.
#   <outdir>  executables in <outdir>/exe, runs in <outdir>/runs/<label>/<test>/<case>,
#             summary in <outdir>/summary_<label>.txt
#
# Also runs the Python replica of the hydro step (replica/run_replica.py) once per outdir.
# Heavy steps (builds, RAMSES runs) run one at a time. Exit status 0 if every case passes.
setopt null_glob
HERE=${0:A:h}
TREE=${1:A}; LABEL=$2; OUT=${3:A}; shift 3
QUICK=(); BUILD=1; TESTS=(advect sod amr tva source restart sink)
while [ $# -gt 0 ]; do
  case $1 in
    --quick) QUICK=(--quick);;
    --nobuild) BUILD=0;;
    --tests) TESTS=(${=2}); shift;;
  esac; shift
done
PY=${PS_PYTHON:-/Users/currodri/.pyenv/versions/calima/bin/python}
EXE=$OUT/exe; RUNS=$OUT/runs; mkdir -p $EXE $RUNS
SUM=$OUT/summary_$LABEL.txt
log() { echo "$*" | tee -a $SUM }
: > $SUM
log "== passive-scalar tests: $LABEL, tree $TREE ($(git -C $TREE log -1 --format='%h %s' | cut -c1-70)), $(date '+%Y-%m-%d %H:%M')"

PREFIX=$LABEL; [ $LABEL = cmaoff ] && PREFIX=cma
if [ $BUILD = 1 ] && [ $LABEL != cmaoff ]; then
  sets=($($PY -c "import sys; sys.path.insert(0,'$HERE/common'); import ps_cases as C; print(' '.join(C.ALL_SETS))"))
  for s in $sets; do
    $HERE/build.sh $TREE $s 1 $EXE/${PREFIX}_${s}_1d | tee -a $SUM
  done
  if (( ${TESTS[(Ie)sink]} )); then
    $HERE/build.sh $TREE full_nodust 3 $EXE/${PREFIX}_full_nodust_3d INDI_STAR=1 | tee -a $SUM
  fi
fi

st=0
for t in $TESTS; do
  if [ $t = restart ]; then
    $PY $HERE/restart/run_restart.py $EXE $LABEL $RUNS > $RUNS/log_${LABEL}_restart.txt 2>&1
    grep -E "^(PASS|FAIL) |restart " $RUNS/log_${LABEL}_restart.txt | sed "s/^/  /" | tee -a $SUM
    grep -q "^PASS " $RUNS/log_${LABEL}_restart.txt || st=1
    continue
  fi
  if [ $t = sink ]; then
    # MPI with 2 ranks: synchronize_sink_info broadcasts from root 1
    $PY $HERE/sink/run_sink.py $EXE/${PREFIX}_full_nodust_3d $RUNS/$LABEL/sink/sink32 --tend 0.01 --levelmin 5 \
        --np 2 --label $LABEL > $RUNS/log_${LABEL}_sink.txt 2>&1
    grep -E "^(PASS|FAIL) |neutral |rel. error" $RUNS/log_${LABEL}_sink.txt | sed "s/^/  /" | tee -a $SUM
    grep -q "^PASS .*with neutral" $RUNS/log_${LABEL}_sink.txt || st=1
    continue
  fi
  $PY $HERE/common/ps_run.py --test $t --exes $EXE --label $LABEL --out $RUNS $QUICK \
      > $RUNS/log_${LABEL}_$t.txt 2>&1 || st=1
  grep -E "^(PASS|FAIL|SKIP) " $RUNS/log_${LABEL}_$t.txt | sed "s/^/  /" >> $SUM
  log "$(grep '^== ' $RUNS/log_${LABEL}_$t.txt)"
  grep -q "^FAIL " $RUNS/log_${LABEL}_$t.txt && st=1
done

if [ ! -f $OUT/replica/replica_summary.txt ]; then
  $PY $HERE/replica/run_replica.py $QUICK --out $OUT/replica > $OUT/replica.log 2>&1
fi
log "$(tail -1 $OUT/replica.log 2>/dev/null)"
$PY $HERE/common/ps_table.py $RUNS $LABEL >> $SUM 2>&1
log "== done $(date +%H:%M), status $st"
exit $st
