# w5e5_downscale / regridding 


Author: Dr Johanna Malle (<mailto:johanna.malle@slf.ch>)

## Summary / purpose and methods


## Main outputs
* **downscaled Wind speed to target domain** (stored as .nc-file)


## Requirements
Before running the script you need to download the required input data files:
### Wind data from Global Wind Atlas

`wget https://globalwindatlas.info/api/gis/global/wind-speed/10/`


### Wind W5E5 data at 0.05deg

`wget https://files.isimip.org/ISIMIP3a/SecondaryInputData/climate/atmosphere/obsclim/global/daily/historical/W5E5v2.0/sfcWind_W5E5v2.0_19790101-19801231.nc`

`wget https://files.isimip.org/ISIMIP3a/SecondaryInputData/climate/atmosphere/obsclim/global/daily/historical/W5E5v2.0/sfcWind_W5E5v2.0_19810101-19901231.nc`

`wget https://files.isimip.org/ISIMIP3a/SecondaryInputData/climate/atmosphere/obsclim/global/daily/historical/W5E5v2.0/sfcWind_W5E5v2.0_19910101-20001231.nc`

`wget https://files.isimip.org/ISIMIP3a/SecondaryInputData/climate/atmosphere/obsclim/global/daily/historical/W5E5v2.0/sfcWind_W5E5v2.0_20010101-20101231.nc`

`wget https://files.isimip.org/ISIMIP3a/SecondaryInputData/climate/atmosphere/obsclim/global/daily/historical/W5E5v2.0/sfcWind_W5E5v2.0_20110101-20191231.nc`


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


## Common problems



## Getting Started

`cd example`

`python example.py`

## Main processing steps
