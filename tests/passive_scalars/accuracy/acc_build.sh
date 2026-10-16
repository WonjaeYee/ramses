#!/bin/bash
# Build the executable of the advection-accuracy tests from a source tree.
#   acc_build.sh <tree> <hydro|mhd> <ndim> <builddir> [make arguments, e.g. COMPILER=INTEL MPIF90=mpiifort]
# Plain build (no RT/RTZ/CALIMA), three passive scalars.  The patch is composed from the tree's own
# <solver>/condinit.f90 (case 'acc_table' added, table_<solver>.f90 appended), so it builds on any
# branch.  MHD uses NENER=1 (the non-thermal energy stays zero): NENER=0 with passive scalars does not
# compile in mhd/godunov_utils.f90 (irad is declared only for NENER>0).  The build directory can be
# anywhere; the executable path is printed.
# bash, not zsh: the nightly suite runs it on machines without zsh
shopt -s nullglob
SOLVER=$2; NDIM=$3
TREE=$(cd "$1" && pwd) || exit 1
mkdir -p "$4/patch" && B=$(cd "$4" && pwd) || exit 1
shift 4
HERE=$(cd "$(dirname "$0")" && pwd)
PY=${PS_PYTHON:-python3}
"$PY" - "$TREE/$SOLVER/condinit.f90" "$HERE/table_$SOLVER.f90" "$B/patch/condinit.f90" <<'PYEOF' || exit 1
import re, sys
src, tab, dst = sys.argv[1:4]
s = open(src).read()
m = re.search(r"^\s*select case\s*\(\s*condinit_kind\s*\)\s*\n", s, re.M | re.I)
assert m, "no 'select case (condinit_kind)' in " + src
s = s[:m.end()] + "\n  case('acc_table')\n     call acc_table_condinit(x, q, dx, nn)\n" + s[m.end():]
open(dst, "w").write(s + "\n" + open(tab).read())
PYEOF
NENER=0; [ "$SOLVER" = mhd ] && NENER=1
vars=(COMPILER=GNU MPIF90=mpif90 "SOLVER=$SOLVER" "NDIM=$NDIM" RT=0 RTZ=0 CALIMA=0 USE_TURB=0 INDI_STAR=0
      "NENER=$NENER" NMETALS=0 NPSCAL=3 CO=0 EXEC=ramses_acc "PATCH=$B/patch")
"$PY" - "$TREE/bin_dust/Makefile" "$B/Makefile" "${vars[@]}" <<'PYEOF' || exit 1
import re, sys
src, dst = sys.argv[1:3]
s = open(src).read()
for a in sys.argv[3:]:
    k, v = a.split("=", 1)
    s, n = re.subn(rf"^{k}\s*\?*=.*$", f"{k} = {v}", s, flags=re.M)
    if n == 0:
        s = f"{k} = {v}\n" + s
tree = src.rsplit("/bin_dust/", 1)[0]
s = re.sub(r"(?<![\w/)])\.\./", tree + "/", s)        # the build directory need not be inside the tree
open(dst, "w").write(s)
PYEOF
(cd "$B" && make "$@" > build.log 2>&1)
exes=("$B"/ramses_acc[12]d)
if [ ${#exes[@]} -eq 0 ] || [ ! -f "${exes[0]}" ]; then
  echo "build $SOLVER ${NDIM}d ($TREE): FAIL $(grep -m1 -E 'Error|No rule|Undefined' "$B/build.log" | cut -c1-160)" >&2; exit 1
fi
echo "${exes[0]}"
