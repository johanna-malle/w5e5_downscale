#!/bin/sh

source activate W5E5_regrid_env

# do bias adjustment for all variables in one program call
python -u w5e5_regrid/downscale_wind.py \
--input_path_w5e5 /home/storage/malle/W5E5_global \
--input_path_wind_atlas /home/storage/malle/CLM5_alps/gwa3_250_wind-speed_10m_alps.tif \
--input_target_grid /home/storage/malle/CLM5_alps/chelsa_input/chelsa-w5e5v1.0_obsclim_mask_30arcsec_global.nc \
--output_loc /home/malle/Desktop/temp_download/test \
--extent_target 42.25,48.75,2.25,17.25 \
--var_interest sfcWind \
--year_start 1979 \
--year_end 1982 \
