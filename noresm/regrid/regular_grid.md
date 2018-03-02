# Converting original grid to regular grid

## Trim duplicate latitude
The model files have 385 latitudes, but the last latitude is a duplicate. Use ferret to remove it:

`save/file=out.nc xxx[j=1:384]`

## Regrid to the 1°x1° grid
`ncremap -m map_tnx1v1_to_woa09_aave_20120501.nc -R '--rgr lat_nm=Y1_384 --rgr lon_nm=X' -i 384_test.nc -o 384_out.nc`

## Apply the coastal fractions
The coastal fraction has to be applied to correct values in partial sea sands. The `frac_b` variable contains this, but it's an inverse value so the reciprocal needs to be applied:

`pco2 * (1 / frac_b)`
