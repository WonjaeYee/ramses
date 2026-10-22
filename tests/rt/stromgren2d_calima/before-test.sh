# Links the CALIMA tables the namelists name locally (calima_data/) and creates SEDtables/, where
# CALIMA writes its mean-opacity tables.
#   RAMSES_RTZ_DATA    rtzdata_for_public_ramses (default: <repo>/rtzdata_for_public_ramses)
rtz=${RAMSES_RTZ_DATA:-../../../rtzdata_for_public_ramses}
[ -d "$rtz/calima_data" ] || { echo "before-test: no CALIMA tables at $rtz/calima_data; set RAMSES_RTZ_DATA"; exit 1; }
ln -sfn "$rtz/calima_data" calima_data
mkdir -p SEDtables
