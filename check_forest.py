# -*- coding: utf-8 -*-
"""
Desc: script to fetch CHELSA data at all selected EVA locations
Created on 10.10.22 15:16
@author: malle
"""

from pathlib import Path
import pandas as pd
import numpy as np
import urllib.request
import xarray as xr
import datetime

bf_hea_all = Path('/home/malle/slfhome/Postdoc3/SAR_biodiv/EVA_request/EVA_selected')

res = '1km'
arcres = '30arcsec'


bf_lakes = Path('/home/lud11/malle/Adhoris/lakes_sector/simstrat/forcing/lakes_coordinats.csv')
eva_in = pd.read_csv(bf_lakes)
eva_in.rename(columns = {'Lake Name in file name':'name',
                         'latitude (dec deg)':'Latitude',
                         'longitude (dec deg)':'Longitude'}, inplace=True)

lat_all = [49.3, 41.85]
lon_all = [18.32, 13.59]
name_all = ['bily', 'colle']
coords_splot_int = list(zip(lat_all, lon_all, name_all))


bf_out = Path('/home/lud11/malle/Adhoris/forest_sector/check_forcing')
bf_out.mkdir(parents=True, exist_ok=True)

year_in_all = range(1979, 1982)
month_in_all = range(1, 13)

var_in_all = ['tas', 'pr']

# create yearly csvs for each variable -> not very efficient but it works...
for var_in in var_in_all:
    bf_data = Path(
        '/home/malle/storage/CHELSA_global/' + var_in + '_' + res + '/files.isimip.org/ISIMIP3a/SecondaryInputData/climate/'
                                                         'atmosphere/obsclim/global/daily/historical/CHELSA-W5E5v1.0')
    print(var_in)
    for year_in in year_in_all:
        print(year_in)
        TSTART = datetime.datetime.now()
        for month_in in month_in_all:
            savename = bf_data / Path(f'chelsa-w5e5v1.0_obsclim_{var_in}_{arcres}_global_daily_{year_in}{month_in:02}.nc')
            file_temp_in = xr.open_dataset(savename)
            print(month_in)

            for splot in coords_splot_int:
                name = splot[2]
                lat = splot[0]
                lon = splot[1]

                data1_lin = file_temp_in.interp(lon=lon, lat=lat, method="nearest")

                file_out1 = bf_out / str(name) / var_in / arcres / f'{var_in}_{year_in}{month_in:02}.csv'
                Path(bf_out, str(name), var_in, arcres).mkdir(parents=True, exist_ok=True)
                if file_out1.is_file():
                    pass
                else:
                    df = data1_lin.to_dataframe()
                    df1 = df.to_csv(file_out1, columns=[var_in])

            file_temp_in.close()
        TEND1 = datetime.datetime.now()
        print(f' Saved all locations for {var_in}, year:{year_in}; took: {TEND1-TSTART} [HH:MM:SS]')
