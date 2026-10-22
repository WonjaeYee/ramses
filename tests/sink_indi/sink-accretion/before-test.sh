# Links the external data the namelist refers to by local names (rtz_data/, mist_sample.unf)
# and creates SEDtables/, where CALIMA writes its mean-opacity tables.
#   RAMSES_RTZ_DATA    rtzdata_for_public_ramses (default: <repo>/rtzdata_for_public_ramses)
#   RAMSES_MIST_SAMPLE MIST sample with SEDs   (default: $RAMSES_RTZ_DATA/260822_mist-sample_v2.5_mcut4.0_with-seds.unf)
rtz=${RAMSES_RTZ_DATA:-../../../rtzdata_for_public_ramses}
mist=${RAMSES_MIST_SAMPLE:-$rtz/260822_mist-sample_v2.5_mcut4.0_with-seds.unf}
[ -d "$rtz/calima_data" ] || { echo "before-test: no RTZ data at $rtz; set RAMSES_RTZ_DATA"; exit 1; }
[ -f "$mist" ] || { echo "before-test: no MIST sample at $mist; set RAMSES_MIST_SAMPLE"; exit 1; }
ln -sfn "$rtz" rtz_data
ln -sf "$mist" mist_sample.unf
mkdir -p SEDtables
