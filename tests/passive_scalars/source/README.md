# source

A uniform 8-cell box at rest with the RTZ chemistry on (`rtz_cooling`, `rtz_include_dust`):
cold dense gas (n=1e3, T~50 K: H2 and CO formation move H into H2 and C, O into CO) and hot
gas (n=1, T~1e6 K, `dust_sputtering`).  No transport, so it isolates the source terms: the
CO exchange and the write-back of cooling_fine (ions times the new element density in the
x*rho_element convention).  PASS: `top`, `child_abs` to round-off.  Limitations: in 1D
`dust_accretion` cannot be used (it forces cmp_sigma_turb, which is 3D-only), and the
sputtering did not change the grains in this setup, so the CALIMA gas-dust exchange itself
is not exercised here.
