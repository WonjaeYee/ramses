#!/bin/zsh
# Build the RAMSES executable of one species set from a source tree.
#   build.sh <tree> <set> <ndim> <exe_out> [extra MAKEVAR=value ...]
# <set> is a key of SPECIES_SETS (common/ps_cases.py), which gives the Makefile variables.
# The build directory is <tree>/bin_ps_<set>_<ndim>d (untracked); the executable is copied
# to <exe_out> (rm first: overwriting a running executable gets it killed on macOS).
# With SOLVER=mhd among the extra variables the build directory is <tree>/bin_ps_mhd_<set>_<ndim>d,
# and the patch is the tree's mhd/condinit.f90 with case('table') and patch/table_condinit_mhd.f90.
setopt null_glob
TREE=${1:A}; SET=$2; NDIM=$3; OUT=${4:A}; shift 4
HERE=${0:A:h}
PY=${PS_PYTHON:-/Users/currodri/.pyenv/versions/calima/bin/python}
ADAPT=${PS_ADAPT:-/Users/currodri/Documents/RAMSES_dev/tests_final/nightly_suite/bin/adapt_makefile.py}
B=$TREE/bin_ps_${SET}_${NDIM}d
PATCH=../tests/passive_scalars/patch
if (( ${@[(Ie)SOLVER=mhd]} )); then
  B=$TREE/bin_ps_mhd_${SET}_${NDIM}d
  PATCH=$B/patch
  mkdir -p $PATCH
  $PY - $TREE/mhd/condinit.f90 $HERE/patch/table_condinit_mhd.f90 $PATCH/condinit.f90 <<'PYEOF' || exit 1
import re, sys
src, tab, dst = sys.argv[1:4]
s = open(src).read()
m = re.search(r"^\s*select case\s*\(\s*condinit_kind\s*\)\s*\n", s, re.M | re.I)
assert m, "no 'select case (condinit_kind)' in " + src
s = s[:m.end()] + "\n  case('table')\n     call table_condinit(x, q, dx, nn)\n" + s[m.end():]
open(dst, "w").write(s + "\n" + open(tab).read())
PYEOF
fi
mkdir -p $B
vars=($($PY -c "import sys; sys.path.insert(0,'$HERE/common'); import ps_cases as C; print(' '.join(f'{k}={v}' for k,v in C.make_vars('$SET',$NDIM).items()))") "$@")
$PY $ADAPT $TREE $TREE/bin_dust/Makefile $B/Makefile ${vars[@]} PATCH=$PATCH > /dev/null || exit 1
mkdir -p ${OUT:h}
log=${OUT}.build.log
(cd $B && timeout 2400 make > $log 2>&1)
exes=($B/ramses*[123]d(N.))
if [ ${#exes} -eq 0 ]; then
  echo "build $SET ${NDIM}d ($TREE): FAIL $(grep -m1 -E 'Error|No rule|Undefined' $log | cut -c1-160)"; exit 1
fi
rm -f $OUT; cp ${exes[1]} $OUT
echo "build $SET ${NDIM}d ($TREE): ok"
