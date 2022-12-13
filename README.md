# w5e5_downscale / bias-correction / regridding 


Author: Johanna Malle (<mailto:johanna.malle@slf.ch>)

## Summary / purpose and methods


## Main outputs
* **downscaled/bias corrected Wind speed @ target domain** (stored as .nc-file)
* **downscaled/bias-corrected relative humidity @ target domain** (stored as .nc-file)
* **bias-corrected surface air pressure @ target domain** (stored as .nc-file)


## Requirements
Before running the script you need to download the required input data files:
### Wind data from Global Wind Atlas

`wget https://globalwindatlas.info/api/gis/global/wind-speed/10/`

### Monthly Relative Humidity Timeseries at 1km (CHELSA)

https://envicloud.wsl.ch/#/?prefix=chelsa%2Fchelsa_V2%2FGLOBAL%2Fmonthly%2F

### Wind W5E5 data at 0.5deg

`wget https://files.isimip.org/ISIMIP3a/SecondaryInputData/climate/atmosphere/obsclim/global/daily/historical/W5E5v2.0/sfcWind_W5E5v2.0_19790101-19801231.nc`

`wget https://files.isimip.org/ISIMIP3a/SecondaryInputData/climate/atmosphere/obsclim/global/daily/historical/W5E5v2.0/sfcWind_W5E5v2.0_19810101-19901231.nc`

`wget https://files.isimip.org/ISIMIP3a/SecondaryInputData/climate/atmosphere/obsclim/global/daily/historical/W5E5v2.0/sfcWind_W5E5v2.0_19910101-20001231.nc`

`wget https://files.isimip.org/ISIMIP3a/SecondaryInputData/climate/atmosphere/obsclim/global/daily/historical/W5E5v2.0/sfcWind_W5E5v2.0_20010101-20101231.nc`

`wget https://files.isimip.org/ISIMIP3a/SecondaryInputData/climate/atmosphere/obsclim/global/daily/historical/W5E5v2.0/sfcWind_W5E5v2.0_20110101-20191231.nc`

### RH W5E5 data at 0.5deg

`wget https://files.isimip.org/ISIMIP3a/SecondaryInputData/climate/atmosphere/obsclim/global/daily/historical/W5E5v2.0/hurs_W5E5v2.0_19790101-19801231.nc`

`wget https://files.isimip.org/ISIMIP3a/SecondaryInputData/climate/atmosphere/obsclim/global/daily/historical/W5E5v2.0/hurs_W5E5v2.0_19810101-19901231.nc`

`wget https://files.isimip.org/ISIMIP3a/SecondaryInputData/climate/atmosphere/obsclim/global/daily/historical/W5E5v2.0/hurs_W5E5v2.0_19910101-20001231.nc`

`wget https://files.isimip.org/ISIMIP3a/SecondaryInputData/climate/atmosphere/obsclim/global/daily/historical/W5E5v2.0/hurs_W5E5v2.0_20010101-20101231.nc`

`wget https://files.isimip.org/ISIMIP3a/SecondaryInputData/climate/atmosphere/obsclim/global/daily/historical/W5E5v2.0/hurs_W5E5v2.0_20110101-20191231.nc`

### sea level pressure W5E5 data at 0.5deg

`wget https://files.isimip.org/ISIMIP3a/SecondaryInputData/climate/atmosphere/obsclim/global/daily/historical/W5E5v2.0/psl_W5E5v2.0_19790101-19801231.nc`

`wget https://files.isimip.org/ISIMIP3a/SecondaryInputData/climate/atmosphere/obsclim/global/daily/historical/W5E5v2.0/psl_W5E5v2.0_19810101-19901231.nc`

`wget https://files.isimip.org/ISIMIP3a/SecondaryInputData/climate/atmosphere/obsclim/global/daily/historical/W5E5v2.0/psl_W5E5v2.0_19910101-20001231.nc`

`wget https://files.isimip.org/ISIMIP3a/SecondaryInputData/climate/atmosphere/obsclim/global/daily/historical/W5E5v2.0/psl_W5E5v2.0_20010101-20101231.nc`

`wget https://files.isimip.org/ISIMIP3a/SecondaryInputData/climate/atmosphere/obsclim/global/daily/historical/W5E5v2.0/psl_W5E5v2.0_20110101-20191231.nc`


## Installation and environment set-up
**Python 3.7 is required.**

Tested on: Ubuntu 18.04 ... to do rest

### Environment set-up
#### With Conda package manager (recommended)
##### Create the environment and install all required packages

`conda env create`

##### Activate the environment

Windows: `activate W5E5_regrid_env`

Linux: `source activate W5E5_regrid_env`


## Getting Started

Update the file application_example with your file paths etc. and then just run:
`./application_example.sh`

## Main processing steps

#### Wind speed
* To include orographic effects into daily mean near-surface wind speed (sfcwind) we follow the approach of Brun et al. 2022, and use an aggregation of the Global Wind Atlas 3.0 data (https://globalwindatlas.info) in combination with daily 0.5° sfcwind from W5E5. 

* We first regrid both the Global Wind Atlas data and the W5E5 sfcwind data to the target grid/extent of 30’’ using bilinear interpolation. 

* The Global Wind Atlas data product represents the period 2008-2017, we therefore average daily regridded W5E5 sfcwind data over this time period. We assume surface wind follows a Weibull distribution and log-transform both datasets before computing the difference layer. We add this difference layer to each log-transformed daily W5E5 raster, and back-transform the sum to obtain the final daily mean near-surface wind speed raster.


#### Relative Humidity
* The provided downscaling algorithm combines monthly 30’’ CHELSA-BIOCLIM+ data (https://doi.org/10.16904/envidat.332 , Brun et al. 2022) with daily W5E5 data
* In a first step we regrid daily 0.5° W5E5 hurs to the target grid (30”, 1.5’ or 5’) and extent using bilinear interpolation. We assume relative humidity follows a beta-distribution and logit-transform both regridded monthly-averaged W5E5 and monthly CHELSA-BIOCLIM+ relative humidity data. The difference layer is then added to daily regridded and logit-transformed W5E5 hurs of the respective month, and the final raster is obtained by back-transforming the sum. 

$$ hurs_{dly} = {1 \over (1+exp^{-h})}  $$

$$ h = log({hurs_{dly}^{W5E5} \over {1 - hurs_{dly}^{W5E5}}} +  Δhurs_{mon}$$  

$$ Δhurs_{mon} = log({{hursmonCHELSA} \over {1 - hurs_{mon}^{CHELSA})} - log({{hurs_{mon}^{W5E5}} \over {1- hurs_{mon}^{W5E5}}} $$  

#### Surface air pressure
* uses W5E5 daily mean sea-level pressure and a DEM to calculate surface air pressure via the barometric formula:

$$ps_{dly} = {psl_{dly}^{W5E5}exp^{-(g * orog * M) / (T_{0} * R)}} $$

with $ps_{dly}$ being the regridded 0.5° W5E5 daily mean sea-level pressure (bilinear interpolation), g the gravitational acceleration constant (9.80665 m/s2), orog the height at which air pressure is calculated (CHELSA-W5E5 orog, m), M the molar mass of dry air (0.02896968 kg/mol), R the universal gas constant (8.314462618 J/(mol K)) and T0 the sea level standard temperature (288.16 K).
