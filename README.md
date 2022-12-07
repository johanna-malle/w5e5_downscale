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

### Wind W5E5 data at 0.05deg

`wget https://files.isimip.org/ISIMIP3a/SecondaryInputData/climate/atmosphere/obsclim/global/daily/historical/W5E5v2.0/sfcWind_W5E5v2.0_19790101-19801231.nc`

`wget https://files.isimip.org/ISIMIP3a/SecondaryInputData/climate/atmosphere/obsclim/global/daily/historical/W5E5v2.0/sfcWind_W5E5v2.0_19810101-19901231.nc`

`wget https://files.isimip.org/ISIMIP3a/SecondaryInputData/climate/atmosphere/obsclim/global/daily/historical/W5E5v2.0/sfcWind_W5E5v2.0_19910101-20001231.nc`

`wget https://files.isimip.org/ISIMIP3a/SecondaryInputData/climate/atmosphere/obsclim/global/daily/historical/W5E5v2.0/sfcWind_W5E5v2.0_20010101-20101231.nc`

`wget https://files.isimip.org/ISIMIP3a/SecondaryInputData/climate/atmosphere/obsclim/global/daily/historical/W5E5v2.0/sfcWind_W5E5v2.0_20110101-20191231.nc`

### RH W5E5 data at 0.05deg

`wget https://files.isimip.org/ISIMIP3a/SecondaryInputData/climate/atmosphere/obsclim/global/daily/historical/W5E5v2.0/hurs_W5E5v2.0_19790101-19801231.nc`

`wget https://files.isimip.org/ISIMIP3a/SecondaryInputData/climate/atmosphere/obsclim/global/daily/historical/W5E5v2.0/hurs_W5E5v2.0_19810101-19901231.nc`

`wget https://files.isimip.org/ISIMIP3a/SecondaryInputData/climate/atmosphere/obsclim/global/daily/historical/W5E5v2.0/hurs_W5E5v2.0_19910101-20001231.nc`

`wget https://files.isimip.org/ISIMIP3a/SecondaryInputData/climate/atmosphere/obsclim/global/daily/historical/W5E5v2.0/hurs_W5E5v2.0_20010101-20101231.nc`

`wget https://files.isimip.org/ISIMIP3a/SecondaryInputData/climate/atmosphere/obsclim/global/daily/historical/W5E5v2.0/hurs_W5E5v2.0_20110101-20191231.nc`

### sea level pressure W5E5 data at 0.05deg

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
* **uses  Global
Wind Atlas 3.0, a free, web-based application developed, owned and operated by the Technical University of Denmark
(https://globalwindatlas.info)
* **We assume Wind follows a Weibull distribution, perform log transform on both the high resolution global wind atlas and the temporal mean of the lower resolution w5e5 data, add this diff. layer to each log-transformed time step and back-transform the sum by exponentiating it

#### Relative Humidity
* **uses 1km monthly CHELSA data (https://chelsa-climate.org/)
* **We assume relative humidity follows a Beta distribution, perform logit transform on both the high resolution chelsa data and the monthly mean of the lower resolution w5e5 data, add this diff. layer to each logit-transformed time step and back-transform the sum

#### Surface air pressure
* **uses W5E5 daily sea level pressure and a DEM to calculate surface air pressure via the barometric formula 
