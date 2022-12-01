# -*- coding: utf-8 -*-
"""
Desc: Prepares and executes the regridding algorithm
Created on 30.11.22 16:01
@author: malle
"""


from optparse import OptionParser
import numpy as np
import xarray as xr
import xesmf as xe
from pathlib import Path
import rioxarray
import pandas as pd
import os
from math import floor, ceil
import datetime
# import urllib.request
# import matplotlib.pyplot as plt
# import warnings


def main():
    """
    Prepares and executes the regridding algorithm.

    """
    # parse command line options and arguments
    parser = OptionParser()
    parser.add_option('-i', '--input_path_w5e5', action='store',
                      type='string', dest='input_path_w5e5', default='',
                      help='path to stored w5e5 data [str]')
    parser.add_option('-w', '--input_path_wind_atlas', action='store',
                      type='string', dest='input_path_wind_atlas', default='',
                      help='path to stored global wind atlas file, if memory problem then clip global file [str]')
    parser.add_option('-t', '--input_target_grid', action='store',
                      type='string', dest='input_target_grid', default='',
                      help='path to stored target grid, incl. file_name [str]')
    parser.add_option('-o', '--output_loc', action='store',
                      type='string', dest='output_loc', default='',
                      help='path to where files should be saved to [str]')
    parser.add_option('-e', '--extent_target', action='store',
                      type='string', dest='extent_target', default='42.25, 48.75, 2.25, 17.25',
                      help='comma-separated list: min_lat, max_lat, min_lon, max_lon')
    parser.add_option('-v', '--var_interest', action='store',
                      type='string', dest='var_interest', default='sfcWind',
                      help='name of variable to be analysed [str]')
    # parser.add_option('-m', '--months', action='store',
    #                  type='string', dest='months', default='1,2,3,4,5,6,7,8,9,10,11,12',
    #                  help=('comma-separated list of integers from {1,...,12} representing months to downscale'))
    parser.add_option('-y', '--year_start', action='store',
                      type='float', dest='year_start', default=1979,
                      help='year to start')
    parser.add_option('-z', '--year_end', action='store',
                      type='float', dest='year_end', default=2016,
                      help='year to end')

    (options, args) = parser.parse_args()

    extent = list(np.array(options.extent_target.split(',')))
    options_min_lon = float(extent[2])
    options_max_lon = float(extent[3])
    options_min_lat = float(extent[0])
    options_max_lat = float(extent[1])

    bf_out = Path(options.output_loc)  # output data will be saved here
    bf_w5e5 = Path(options.input_path_w5e5)
    var_in = options.var_interest
    (bf_out / var_in).mkdir(parents=True, exist_ok=True)

    target_grid = xr.open_dataset(options.input_target_grid)
    clipped_target_grid = target_grid.loc[
        {'lat': slice(options_min_lat, options_max_lat),
         'lon': slice(options_min_lon, options_max_lon)}]
    ds_target = xr.Dataset({"lat": (["lat"], clipped_target_grid.lat.data),
                            "lon": (["lon"], clipped_target_grid.lon.data)})

    # step 1: cut global wind speed to our extent and regrid (upscale) it to the target grid already
    wind_windatlas_file = Path(options.input_path_wind_atlas)
    # source: https://globalwindatlas.info/en
    wind_windatlas = rioxarray.open_rasterio(wind_windatlas_file, masked_and_scale=True)\
        .rio.clip_box(minx=options_min_lon, miny=options_min_lat, maxx=options_max_lon, maxy=options_max_lat)
    wind_windatlas1 = wind_windatlas.rename({'x': 'lon', 'y': 'lat'}).isel(band=0).to_dataset(name='sfcWind')

    regridder_coarse = xe.Regridder(wind_windatlas1, ds_target, "bilinear")
    dr_out_wind_atlas = regridder_coarse(wind_windatlas1)

    # step 2: cut w5e5 data to our extent and average over time to match global wind atlas timeframe (2008-2017)
    all_files = list(Path(bf_w5e5 / var_in).iterdir())
    w5e5_data = xr.open_mfdataset(all_files)
    clipped_w5e5 = w5e5_data.loc[{'lat': slice(ceil(options_max_lat), floor(options_min_lat)),
                                  'lon': slice(floor(options_min_lon), ceil(options_max_lon))}]
    clipped_w5e5_time = clipped_w5e5.loc[{'time': slice(np.datetime64('2008-01-01'), np.datetime64('2017-12-31'))}]
    w5e5_avg_time = clipped_w5e5_time.mean(dim='time')

    regridder = xe.Regridder(w5e5_avg_time, ds_target, "bilinear")
    dr_out_w5e5_avg = regridder(w5e5_avg_time)

    # step 3: create diff layer:
    log_w5e5 = np.log(dr_out_w5e5_avg.sfcWind.values, out=np.zeros_like(dr_out_w5e5_avg.sfcWind.values),
                      where=(dr_out_w5e5_avg.sfcWind.values != 0))  # avoid  -inf at locations where wind == 0

    log_wind_atlas = np.log(dr_out_wind_atlas.sfcWind.values, out=np.zeros_like(dr_out_wind_atlas.sfcWind.values),
                            where=(dr_out_wind_atlas.sfcWind.values != 0))  # avoid  -inf at locations where wind == 0

    diff_layer = log_wind_atlas - log_w5e5   # potentially I could just calc this layer globally for a few resol.?

    # step 4: add this diff layer to each time step => perform analysis on monthly basis
    datetimes = pd.to_datetime(clipped_w5e5['time'])
    year_all = range(np.int32(options.year_start),np.int32(options.year_end+1))

    for year_int in year_all:
        file_save = bf_out / var_in / Path(var_in + '_downscaled_v1.0_'+str(year_int)+'.nc')
        id_all = datetimes.year == year_int
        clip_w5e5_time = clipped_w5e5.loc[{'time': slice(datetimes[id_all][0], datetimes[id_all][-1])}]
        clip_w5e5_time_regrid = regridder(clip_w5e5_time)
        clip_log_corr = np.log(clip_w5e5_time_regrid.sfcWind.values,
                               out=np.zeros_like(clip_w5e5_time_regrid.sfcWind.values),
                               where=(clip_w5e5_time_regrid.sfcWind.values != 0)) + diff_layer   # avoid  -inf at x=0
        w5e5_cor = np.exp(clip_log_corr, out=np.zeros_like(clip_log_corr), where=(clip_log_corr != 0))
        final_datset = clip_w5e5_time_regrid.copy()
        final_datset.sfcWind.values = w5e5_cor

        final_datset.attrs = {
            'history': 'File was created by ' + os.getlogin() + ', PC ' + os.uname().nodename + ', created on ' +
                       datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S'),
            'institution': 'Swiss Federal Institute for Forest, Snow and Landscape Research WSL',
            'contact': 'J.Malle <johanna.malle@wsl.ch>',
            'title': 'Downscaled W5E5 wind data to custom domain CHELSA data',
            'version': '1.0',
            'project': 'Inter-Sectoral Impact Model Intercomparison Project phase 3 (ISIMIP3), HighRes experiments',
            'based on': 'Global wind atlas (https://globalwindatlas.info/en)'}

        final_datset.time.attrs = clipped_w5e5.time.attrs
        final_datset.lon.attrs = clipped_w5e5.lon.attrs
        final_datset.lat.attrs = clipped_w5e5.lat.attrs
        final_datset.sfcWind.attrs = clipped_w5e5.sfcWind.attrs

        final_datset.to_netcdf(file_save)


if __name__ == '__main__':
    main()
