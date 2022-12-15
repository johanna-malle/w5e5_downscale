w5e5_regrid

#!/bin/sh

source activate /home/malle/miniconda3/envs/W5E5_regrid_env
cd /home/malle/slfhome/Postdoc3/python-scripts/pythonProject/adhoris_regrid/

# do bias adjustment for all variables in one program call
python -u w5e5_regrid/downscale_all.py \
--input_path_w5e5 /storage/malle/W5E5_global \
--input_path_wind_atlas /storage/malle/CLM5_alps/gwa3_250_wind-speed_10m_alps.tif \
--input_target_grid /storage/malle/CLM5_alps/chelsa_input/chelsa-w5e5v1.0_obsclim_mask_30arcsec_global.nc \
--output_loc /home/malle/test_run \
--extent_target 42.25,48.75,2.25,17.25 \
--var_interest ps,sfcWind,rh, rlds \
--year_start 1981 \
--year_end 1982 \
--input_path_orog /storage/malle/CLM5_alps/chelsa_input/chelsa-w5e5v1.0_obsclim_orog_30arcsec_global.nc \
--input_path_chelsa /storage/brunp/Data/CHELSA_V.2.1/RH/Monthly \
--input_path_tas_fine /home/malle/Desktop/temp_download/test/tas \
--input_path_rh_fine /home/malle/Desktop/temp_download/test/rh \
