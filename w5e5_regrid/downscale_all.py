# -*- coding: utf-8 -*-
"""
Desc: Bias correction of Wind, relative humidity, and surface air pressure
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
import time
import glob


def downscale_wind(wind_windatlas_file, ds_target, extent, bf_w5e5, bf_out, years_all, var_in):
    """
    Downscales wind based on global wind atlas data

    Parameters
    ----------
    wind_windatlas_file : Path
        Location of global wind atlas file.
    ds_target : xarray Dataset
        Target grid, as xr Dataset .
    extent : list
        comma seperated list of strings:
        extent of interest [min_lat, max_lat, min_lon, max_lon]
    bf_w5e5 : Path
        Location of w5e5 base folder.
    bf_out : str, optional
        Location of output base folder.
    years_all : ndarray
        int of start and end of regridding [year_start, year_end]
    var_in : str
        which variable are we dealing with? ['sfcWind','ps','rh']

    Returns
    -------

    """

    # step 1: cut global wind speed to our extent and regrid (upscale) it to the target grid already
    wind_windatlas = rioxarray.open_rasterio(wind_windatlas_file, masked_and_scale=True)\
        .rio.clip_box(minx=float(extent[2]), miny=float(extent[0]), maxx=float(extent[1]), maxy=float(extent[1]))
    wind_windatlas1 = wind_windatlas.rename({'x': 'lon', 'y': 'lat'}).isel(band=0).to_dataset(name='sfcWind')

    regridder_coarse = xe.Regridder(wind_windatlas1, ds_target, "bilinear")
    dr_out_wind_atlas = regridder_coarse(wind_windatlas1)

    # step 2: cut w5e5 data to our extent and average over time to match global wind atlas timeframe (2008-2017)
    all_files = list(Path(bf_w5e5 / var_in).iterdir())
    w5e5_data = xr.open_mfdataset(all_files)
    clipped_w5e5 = w5e5_data.loc[{'lat': slice(ceil(float(extent[1])), floor(float(extent[0]))),
                                  'lon': slice(floor(float(extent[2])), ceil(float(extent[3])))}]
    clipped_w5e5_time = clipped_w5e5.loc[{'time': slice(np.datetime64('2008-01-01'), np.datetime64('2017-12-31'))}]
    w5e5_avg_time = clipped_w5e5_time.mean(dim='time')

    regridder = xe.Regridder(w5e5_avg_time, ds_target, "bilinear")
    dr_out_w5e5_avg = regridder(w5e5_avg_time)

    # step 3: create diff layer:
    # assume wind follows weibull distribution => do log transform
    log_w5e5 = np.log(dr_out_w5e5_avg.sfcWind.values, out=np.zeros_like(dr_out_w5e5_avg.sfcWind.values),
                      where=(dr_out_w5e5_avg.sfcWind.values != 0))  # avoid  -inf at locations where wind == 0

    log_wind_atlas = np.log(dr_out_wind_atlas.sfcWind.values, out=np.zeros_like(dr_out_wind_atlas.sfcWind.values),
                            where=(dr_out_wind_atlas.sfcWind.values != 0))  # avoid  -inf at locations where wind == 0

    diff_layer = log_wind_atlas - log_w5e5   # to be added to log-transformed daily
    # potentially I could just calc this layer globally for a few resol.?

    # step 4: add this diff layer to each time step
    datetimes = pd.to_datetime(clipped_w5e5['time'])
    year_all = range(years_all[0], years_all[1]+1)

    for year_int in year_all:  # yearly analysis to avoid files becoming too big
        file_save = bf_out / var_in / Path(var_in + '_corr_v1.0_'+str(year_int)+'.nc')
        id_all = datetimes.year == year_int
        clip_w5e5_time = clipped_w5e5.loc[{'time': slice(datetimes[id_all][0], datetimes[id_all][-1])}]
        clip_w5e5_time_regrid = regridder(clip_w5e5_time)
        clip_log_corr = np.log(clip_w5e5_time_regrid.sfcWind.values,
                               out=np.zeros_like(clip_w5e5_time_regrid.sfcWind.values),
                               where=(clip_w5e5_time_regrid.sfcWind.values != 0)) + diff_layer   # avoid  -inf at x=0
        w5e5_cor = np.exp(clip_log_corr, out=np.zeros_like(clip_log_corr), where=(clip_log_corr != 0))
        final_datset = clip_w5e5_time_regrid.copy()
        final_datset.sfcWind.values = w5e5_cor

        final_datset.time.attrs = clipped_w5e5.time.attrs
        final_datset.lon.attrs = clipped_w5e5.lon.attrs
        final_datset.lat.attrs = clipped_w5e5.lat.attrs
        final_datset.sfcWind.attrs = clipped_w5e5.sfcWind.attrs

        final_datset.attrs = {
            'history': 'File was created by ' + os.getlogin() + ', PC ' + os.uname().nodename + ', created on ' +
                       datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S'),
            'institution': 'Swiss Federal Institute for Forest, Snow and Landscape Research WSL',
            'contact': 'J.Malle <johanna.malle@wsl.ch>',
            'title': 'Downscaled daily W5E5 ' + var_in + ' data to custom domain',
            'version': '1.0',
            'project': 'Inter-Sectoral Impact Model Intercomparison Project phase 3 (ISIMIP3), HighRes experiments',
            'based on': 'Global wind atlas (https://globalwindatlas.info/en)'}

        final_datset.to_netcdf(file_save)
        print('done with year ' + str(year_int))


def downscale_pressure(orog_file, ds_target, extent, bf_w5e5, bf_out, years_all, var_in):
    """
    Calculates surface air pressure at custom grid
    based on daily w5e5 sealevel pressure and altitude,
    using the barometric formula

    Parameters
    ----------
    orog_file : Path
        Location of dem
    ds_target : xarray Dataset
        Target grid, as xr Dataset
    extent : list
        comma seperated list of strings:
        extent of interest [min_lat, max_lat, min_lon, max_lon]
    bf_w5e5 : Path
        Location of w5e5 base folder.
    bf_out : str, optional
        Location of output base folder.
    years_all : ndarray
        int of start and end of regridding [year_start, year_end]
    var_in : str
        which variable are we dealing with? ['sfcWind','ps','rh']

    Returns
    -------

    """

    # constants for barometric formula:
    print('start downscale pressure')
    g = 9.80665  # gravitational acceleration [m/s2]
    M = 0.02896968  # molar mass of dry air [kg/mol]
    r0 = 8.314462618  # universal gas constant [J/(mol·K)]
    T0 = 288.16  # Sea level standard temperature  [K]

    all_files = list(Path(bf_w5e5 / 'psl').iterdir())
    w5e5_data = xr.open_mfdataset(all_files)
    clipped_w5e5 = w5e5_data.loc[{'lat': slice(ceil(float(extent[1])), floor(float(extent[0]))),
                                  'lon': slice(floor(float(extent[2])), ceil(float(extent[3])))}]

    regridder_w5e5 = xe.Regridder(clipped_w5e5, ds_target, "bilinear")

    file_temp_orog = xr.open_dataset(orog_file)
    clipped_orog = file_temp_orog.loc[
        {'lat': slice(floor(float(extent[0])), ceil(float(extent[1]))),
         'lon': slice(floor(float(extent[2])), ceil(float(extent[3])))}]
    del file_temp_orog

    regridder_orog = xe.Regridder(clipped_orog, ds_target, "bilinear")
    dr_out_orog = regridder_orog(clipped_orog)
    del clipped_orog

    datetimes = pd.to_datetime(clipped_w5e5['time'])
    year_all = range(years_all[0], years_all[1]+1)
    print('done with regridding... now applying')

    for year_int in year_all:  # yearly analysis to avoid files becoming too big
        file_save = bf_out / var_in / Path(var_in + '_corr_v1.0_'+str(year_int)+'.nc')
        id_all = datetimes.year == year_int
        clip_w5e5_time = clipped_w5e5.loc[{'time': slice(datetimes[id_all][0], datetimes[id_all][-1])}]
        clip_w5e5_time_regrid = regridder_w5e5(clip_w5e5_time)

        ps_air = clip_w5e5_time_regrid * np.exp(-(g * dr_out_orog.to_array() * M) / (T0 * r0))

        final_datset = clip_w5e5_time_regrid.copy()
        final_datset.psl.values = np.squeeze(ps_air.psl)
        final_datset.time.attrs = clipped_w5e5.time.attrs
        final_datset.lon.attrs = clipped_w5e5.lon.attrs
        final_datset.lat.attrs = clipped_w5e5.lat.attrs
        final_datset.psl.attrs = clipped_w5e5.psl.attrs
        final_datset.psl.attrs = {'long_name': 'Surface Air Pressure', 'standard_name': 'surface_air_pressure'}
        final_datset1 = final_datset.rename({'psl': 'ps'})
        del final_datset

        final_datset1.attrs = {
            'history': 'File was created by ' + os.getlogin() + ', PC ' + os.uname().nodename + ', created on ' +
                       datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S'),
            'institution': 'Swiss Federal Institute for Forest, Snow and Landscape Research WSL',
            'contact': 'J.Malle <johanna.malle@wsl.ch>',
            'title': 'Surface Air Pressure: Downscaled daily W5E5 ' + var_in + ' data to custom domain',
            'version': '1.0',
            'project': 'Inter-Sectoral Impact Model Intercomparison Project phase 3 (ISIMIP3), HighRes experiments',
            'info': 'calculations based on barometric formula'}

        final_datset1.to_netcdf(file_save)
        print('done with year ' + str(year_int))


def downscale_rh(path_monthly_chelsa, ds_target, extent, bf_w5e5, bf_out, years_all, var_in):
    """
    Calculates surface air pressure based on w5e5 sealevel pressure and altitude, using the barometric formula

    Parameters
    ----------
    path_monthly_chelsa : Path
        Location of monthly 1km chelsa data
    ds_target : xarray Dataset
        Target grid, as xr Dataset
    extent : list
        comma seperated list of strings:
        extent of interest [min_lat, max_lat, min_lon, max_lon]
    bf_w5e5 : Path
        Location of w5e5 base folder.
    bf_out : str, optional
        Location of output base folder.
    years_all : ndarray
        int of start and end of regridding [year_start, year_end]
    var_in : str
        which variable are we dealing with? ['sfcWind','ps','rh']

    Returns
    -------

    """

    all_files = list(Path(bf_w5e5 / 'hurs').iterdir())
    w5e5_data = xr.open_mfdataset(all_files)
    clipped_w5e5 = w5e5_data.loc[{'lat': slice(ceil(float(extent[1])), floor(float(extent[0]))),
                                  'lon': slice(floor(float(extent[2])), ceil(float(extent[3])))}]

    regridder_w5e5 = xe.Regridder(clipped_w5e5, ds_target, "bilinear")

    datetimes = pd.to_datetime(clipped_w5e5['time'])
    year_all = range(years_all[0], years_all[1]+1)

    # open first chelsa file to generate regridding weights
    file_in_chelsa = rioxarray.open_rasterio(path_monthly_chelsa / 'CHELSA_hurs_01_1981_V.2.1.tif',
                                             mask_and_scale=True).rio.clip_box(minx=float(extent[2]),
                                                                               miny=float(extent[3]),
                                                                               maxx=float(extent[0]),
                                                                               maxy=float(extent[1]))

    regridder_chelsa = xe.Regridder(file_in_chelsa, ds_target, "bilinear")

    print('done with regridding... now applying to all timesteps')

    for year_int in year_all:  # yearly analysis to avoid files becoming too big
        file_save = bf_out / var_in / Path(var_in + '_corr_v1.0_'+str(year_int)+'.nc')
        id_all = datetimes.year == year_int
        clip_w5e5_time = clipped_w5e5.loc[{'time': slice(datetimes[id_all][0], datetimes[id_all][-1])}]
        clip_w5e5_time_re = regridder_w5e5(clip_w5e5_time)
        monthly_rh_cor = clip_w5e5_time_re.copy()  # values to be overwritten
        datetimes_year = pd.to_datetime(clip_w5e5_time['time'])

        for month_int in range(1, 13):
            id_all = datetimes_year.month == month_int
            month_w5e5 = clip_w5e5_time.loc[{'time': slice(datetimes_year[id_all][0], datetimes_year[id_all][-1])}]

            file_in_chelsa = glob.glob(str(path_monthly_chelsa) + '/*' + f"{month_int:02}" + '_' + str(year_int) + '*')
            month_chelsa = rioxarray.open_rasterio(file_in_chelsa[0], mask_and_scale=True).rio.clip_box(
                                                                               minx=float(extent[2]),
                                                                               miny=float(extent[3]),
                                                                               maxx=float(extent[0]),
                                                                               maxy=float(extent[1]))

            dr_out_chelsa = regridder_chelsa(month_chelsa) * 0.01  # before transformation: change rh to vary b/w 0 1
            dr_out_chelsa_tr = np.log(dr_out_chelsa / (1 - dr_out_chelsa))  # assume beta distrib. => logit transform

            dr_out_w5e5 = regridder_w5e5(month_w5e5) * 0.01  # before transformation: change rh to vary b/w 0 and 1
            dr_out_w5e5_mean = dr_out_w5e5.mean(dim='time')  # monthly average for diff layer b/c chelsa == monthly
            dr_out_w5e5_tr = np.log(dr_out_w5e5 / (1 - dr_out_w5e5))  # assume beta distribuation => logit transform
            dr_out_w5e5_mean_tr = np.log(dr_out_w5e5_mean / (1 - dr_out_w5e5_mean))  # logit transform

            monthly_diff_layer = dr_out_chelsa_tr - dr_out_w5e5_mean_tr
            new_transf_rh1 = dr_out_w5e5_tr + monthly_diff_layer
            monthly_rh_cor1 = (1 / (1 + np.exp(-new_transf_rh1))) * 100  # back transform

            id_all_month = np.logical_and(datetimes.year == year_int, datetimes.month == month_int)
            monthly_rh_cor['hurs'].loc[{'time': slice(datetimes[id_all_month][0], datetimes[id_all_month][-1])}] = \
                monthly_rh_cor1.hurs.isel(band=0).squeeze().values

        monthly_rh_cor.attrs = {
            'history': 'File was created by ' + os.getlogin() + ', PC ' + os.uname().nodename + ', created on ' +
                       datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S'),
            'institution': 'Swiss Federal Institute for Forest, Snow and Landscape Research WSL',
            'contact': 'J.Malle <johanna.malle@wsl.ch>',
            'title': 'Downscaled daily W5E5 ' + var_in + ' data to custom domain',
            'version': '1.0',
            'project': 'Inter-Sectoral Impact Model Intercomparison Project phase 3 (ISIMIP3), HighRes experiments',
            'based on': 'monthly chelsa data'}

        monthly_rh_cor.to_netcdf(file_save)
        print('done with year ' + str(year_int))


def downscale_longwave(path_rh_fine, path_tas_fine, ds_target, extent, bf_w5e5, bf_out, years_all, var_in):
    """
    Calculates downscaled longwave radiation based on w5e5 longwave etc

    Parameters
    ----------
    path_rh_fine : Path
        Location of rh high resolution data
    path_tas_fine : Path
        Location of tas high resolution data
    ds_target : xarray Dataset
        Target grid, as xr Dataset
    extent : list
        comma seperated list of strings:
        extent of interest [min_lat, max_lat, min_lon, max_lon]
    bf_w5e5 : Path
        Location of w5e5 base folder.
    bf_out : str, optional
        Location of output base folder.
    years_all : ndarray
        int of start and end of regridding [year_start, year_end]
    var_in : str
        which variable are we dealing with? ['sfcWind','ps','rh']

    Returns
    -------

    """

    x1 = 0.43
    x2 = 5.7
    sbc = 5.67 * 10E-8   # stefan boltzman constant [Js−1 m−2 K−4]

    es0 = 6.11  # reference saturation vapour pressure
    T0 = 273.15
    lv = 2.5E6  # latent heat of vaporization of water
    Rv = 461.5  # gas constant for water vapour [J K kg-1]

    all_files = list(Path(bf_w5e5 / 'hurs').iterdir())
    w5e5_data = xr.open_mfdataset(all_files)
    clipped_w5e5_rh = w5e5_data.loc[{'lat': slice(ceil(float(extent[1])), floor(float(extent[0]))),
                                     'lon': slice(floor(float(extent[2])), ceil(float(extent[3])))}]

    all_files = list(Path(bf_w5e5 / 'tas').iterdir())
    w5e5_data = xr.open_mfdataset(all_files)
    clipped_w5e5_tas = w5e5_data.loc[{'lat': slice(ceil(float(extent[1])), floor(float(extent[0]))),
                                      'lon': slice(floor(float(extent[2])), ceil(float(extent[3])))}]

    all_files = list(Path(bf_w5e5 / 'rlds').iterdir())
    w5e5_data = xr.open_mfdataset(all_files)
    clipped_w5e5_rlds = w5e5_data.loc[{'lat': slice(ceil(float(extent[1])), floor(float(extent[0]))),
                                       'lon': slice(floor(float(extent[2])), ceil(float(extent[3])))}]

    regridder_w5e5 = xe.Regridder(clipped_w5e5_rh, ds_target, "bilinear")

    datetimes = pd.to_datetime(clipped_w5e5_rh['time'])
    year_all = range(years_all[0], years_all[1]+1)

    # we also need target grid rh and temp:
    file_in_tas_fine = glob.glob(str(path_tas_fine) + '/*' + str(year_all[0]) + '*')
    tas_in_fine = xr.open_dataset(file_in_tas_fine[0])
    clipped_tas_fine = tas_in_fine.loc[{'lat': slice(floor(float(extent[0])), ceil(float(extent[1]))),
                                        'lon': slice(floor(float(extent[2])), ceil(float(extent[3])))}]
    del tas_in_fine
    regridder_chelsa_tas = xe.Regridder(clipped_tas_fine, ds_target, "bilinear")

    file_in_rh_fine = glob.glob(str(path_rh_fine) + '/*' + str(year_all[0]) + '*')
    rh_in_fine = xr.open_dataset(file_in_rh_fine[0])
    clipped_rh_fine = rh_in_fine.loc[{'lat': slice(floor(float(extent[0])), ceil(float(extent[1]))),
                                      'lon': slice(floor(float(extent[2])), ceil(float(extent[3])))}]
    regridder_chelsa_rh = xe.Regridder(clipped_rh_fine, ds_target, "bilinear")

    print('done with regridding... now applying to all timesteps')

    for year_int in year_all:  # yearly analysis to avoid files becoming too big
        file_save = bf_out / var_in / Path(var_in + '_corr_v1.0_'+str(year_int)+'.nc')
        (bf_out / var_in).mkdir(parents=True, exist_ok=True)
        id_all = datetimes.year == year_int

        clip_w5e5_time = clipped_w5e5_rh.loc[{'time': slice(datetimes[id_all][0], datetimes[id_all][-1])}]
        rh_coarse = regridder_w5e5(clip_w5e5_time)

        clip_w5e5_time = clipped_w5e5_tas.loc[{'time': slice(datetimes[id_all][0], datetimes[id_all][-1])}]
        temp_coarse = regridder_w5e5(clip_w5e5_time)

        clip_w5e5_rlds = clipped_w5e5_rlds.loc[{'time': slice(datetimes[id_all][0], datetimes[id_all][-1])}]
        lw_coarse = regridder_w5e5(clip_w5e5_rlds)

        file_in_rh_fine = glob.glob(str(path_rh_fine) + '/*' + str(year_int) + '*')
        file_rh_fine = xr.open_mfdataset(file_in_rh_fine)  # should be just 1 file
        clipped_rh_fine = file_rh_fine.loc[{'lat': slice(floor(float(extent[0])), ceil(float(extent[1]))),
                                            'lon': slice(floor(float(extent[2])), ceil(float(extent[3])))}]
        rh_fine = regridder_chelsa_rh(clipped_rh_fine)

        file_in_tas_fine = glob.glob(str(path_tas_fine) + '/*' + str(year_int) + '*')
        file_tas_fine = xr.open_mfdataset(file_in_tas_fine)  # multiple files
        clipped_tas_fine = file_tas_fine.loc[{'lat': slice(floor(float(extent[0])), ceil(float(extent[1]))),
                                              'lon': slice(floor(float(extent[2])), ceil(float(extent[3])))}]
        temp_fine = regridder_chelsa_tas(clipped_tas_fine)

        # now ready for calculation:
        es_coarse = es0 * np.exp((lv / Rv) * (1 / T0 - 1 / temp_coarse.tas.data))
        pV_coarse = (rh_coarse.hurs.data * es_coarse) / 100  # water vapur pressure

        es_fine = es0 * np.exp((lv / Rv) * (1 / T0 - 1 / temp_fine.tas.data))
        pV_fine = (rh_fine.hurs.data * es_fine) / 100  # water vapur pressure

        e_cl_coarse = 0.23 + x1 * (pV_coarse / temp_coarse.tas.data) ** (1 / x2)  # clear-sky emissivity w5e5
        e_cl_fine = 0.23 + x1 * (pV_fine / temp_fine.tas.data) ** (1 / x2)  # clear-sky emissivity target grid

        e_as_coarse = lw_coarse.rlds.data / (sbc * temp_coarse.tas.data ** 4)  # all-sky emissivity w5e5
        delta_e = e_as_coarse - e_cl_coarse  # cloud-based component of emissivity w5e5

        lw_fine = (e_cl_fine + delta_e) * sbc * temp_fine.tas.data ** 4  # downscaled lwr! assume cloud e is the same

        final_datset = lw_coarse.copy()
        final_datset.rlds.values = lw_fine

        final_datset.time.attrs = clip_w5e5_rlds.time.attrs
        final_datset.lon.attrs = clip_w5e5_rlds.lon.attrs
        final_datset.lat.attrs = clip_w5e5_rlds.lat.attrs
        final_datset.rlds.attrs = clip_w5e5_rlds.rlds.attrs

        final_datset.attrs = {
            'history': 'File was created by ' + os.getlogin() + ', PC ' + os.uname().nodename + ', created on ' +
                       datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S'),
            'institution': 'Swiss Federal Institute for Forest, Snow and Landscape Research WSL',
            'contact': 'J.Malle <johanna.malle@wsl.ch>',
            'title': 'Downscaled daily W5E5 ' + var_in + ' data to custom domain',
            'version': '1.0',
            'project': 'Inter-Sectoral Impact Model Intercomparison Project phase 3 (ISIMIP3a), HighRes experiments'}

        final_datset.to_netcdf(file_save)
        print('done with year ' + str(year_int))


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
    parser.add_option('-d', '--input_path_orog', action='store',
                      type='string', dest='input_path_orog', default='',
                      help='path to stored elevational height [str]')
    parser.add_option('-c', '--input_path_chelsa', action='store',
                      type='string', dest='input_path_chelsa', default='',
                      help='path to stored elevational height [str]')
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
                      help='name of variable to be analysed, comma separated if more than 1 [str]')
    parser.add_option('-v', '--input_path_tas_fine', action='store',
                      type='string', dest='input_path_tas_fine', default='',
                      help='path to highres tas for lwr computation [str]')
    parser.add_option('-v', '--input_path_rh_fine', action='store',
                      type='string', dest='input_path_rh_fine', default='',
                      help='path to highres rh for lwr computation [str]')
    parser.add_option('-y', '--year_start', action='store',
                      type='int', dest='year_start', default=1979,
                      help='year to start')
    parser.add_option('-z', '--year_end', action='store',
                      type='int', dest='year_end', default=2016,
                      help='year to end')
    # parser.add_option('-m', '--months', action='store',
    #                  type='string', dest='months', default='1,2,3,4,5,6,7,8,9,10,11,12',
    #                  help=('comma-separated list of integers from {1,...,12} representing months to downscale'))

    (options, args) = parser.parse_args()

    tt = time.time()
    timeit2 = 'Done with all timesteps for this variable [{:.3f}s]'
    timeit3 = 'Saved all requested variables and timesteps [total time = {:.3f}s]'

    extent = list(np.array(options.extent_target.split(',')))
    options_min_lon = float(extent[2])
    options_max_lon = float(extent[3])
    options_min_lat = float(extent[0])
    options_max_lat = float(extent[1])

    target_grid = xr.open_dataset(options.input_target_grid)
    clipped_target_grid = target_grid.loc[
        {'lat': slice(options_min_lat, options_max_lat),
         'lon': slice(options_min_lon, options_max_lon)}]
    ds_target = xr.Dataset({"lat": (["lat"], clipped_target_grid.lat.data),
                            "lon": (["lon"], clipped_target_grid.lon.data)})

    bf_out = Path(options.output_loc)  # output data will be saved here
    bf_w5e5 = Path(options.input_path_w5e5)
    var_in_all = options.var_interest.split(',')

    years_all = np.array([options.year_start, options.year_end])

    for var_in in var_in_all:
        tt1 = time.time()
        (bf_out / var_in).mkdir(parents=True, exist_ok=True)
        print('Starting with variable ' + var_in)
        if var_in == 'sfcWind':
            wind_windatlas_file = Path(options.input_path_wind_atlas)  # source: https://globalwindatlas.info/en
            downscale_wind(wind_windatlas_file, ds_target, extent, bf_w5e5, bf_out, years_all, var_in)
        elif var_in == 'ps':
            orog_file = Path(options.input_path_orog)
            downscale_pressure(orog_file, ds_target, extent, bf_w5e5, bf_out, years_all, var_in)
        elif var_in == 'rh':
            path_monthly_chelsa = Path(options.input_path_chelsa)
            downscale_rh(path_monthly_chelsa, ds_target, extent, bf_w5e5, bf_out, years_all, var_in)
        elif var_in == 'rlds':
            path_rh_fine = Path(options.input_path_rh_fine)
            path_tas_fine = Path(options.input_path_tas_fine)
            downscale_longwave(path_rh_fine, path_tas_fine, ds_target, extent, bf_w5e5, bf_out, years_all, var_in)

        print(timeit2.format(time.time() - tt1))

    print(timeit3.format(time.time() - tt))


if __name__ == '__main__':
    main()
