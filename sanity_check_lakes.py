# -*- coding: utf-8 -*-
"""
Desc:
Created on 07.05.23 16:02
@author: malle
"""

import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
import xarray as xr
import matplotlib
import platform
import os
from pathlib import Path
import pandas as pd
import seaborn as sns


bf = Path('/home/malle/Desktop/test_package/output_lakes/')

var_all = ['sfcWind', 'ps', 'rlds', 'rh']
res_all = ['1800arcsec', '300arcsec', '90arcsec', '30arcsec']
all_lakes = ['allequash', 'alqueva', 'annie', 'argyle', 'biel', 'big-muskellunge', 'black-oak', 'burley-griffin', 'crystal-lake', 'crystal-bog', 'delavan', 'dickie', 'eagle', 'ekoln', 'erken', 'esthwaite-water', 'falling-creek', 'feeagh', 'fish', 'great', 'green', 'harp', 'kilpisjarvi', 'kinneret', 'kivu', 'klicava', 'kuivajarvi', 'langtjern', 'laramie', 'lower-zurich', 'mendota', 'monona', 'mozhaysk', 'mt-bold', 'mueggelsee', 'neuchatel', 'ngoring', 'nohipalo-mustjaerv', 'nohipalo-valgejaerv', 'okauchee', 'paaijarvi', 'rappbode', 'rimov', 'rotorua', 'sammamish', 'sau', 'sparkling', 'stechlin', 'sunapee', 'tahoe', 'tarawera', 'toolik', 'trout-lake', 'trout-bog', 'two-sisters', 'vendyurskoe', 'vortsjaerv', 'washington', 'windermere', 'wingra', 'zlutice', 'taihu', 'chao', 'hulun', 'rappbode', 'hassel', 'arendsee', 'scharmutzelsee', 'zurich', 'thun', 'murten', 'bryrup', 'bosumtwi']
bf_out_figures = bf / 'figures_sanity_check'
bf_out_figures.mkdir(parents=True, exist_ok=True)
col_res = ['orange', (161 / 255, 218 / 255, 180 / 255),
           (65 / 255, 182 / 255, 196 / 255), (34 / 255, 94 / 255, 168 / 255)]

bf_ref = Path('/home/malle/Desktop/test_package/output_lakes/ref')

for lake_in in all_lakes:
    print(lake_in)
    bf_30 = bf / '30arcsec' / lake_in
    bf_90 = bf / '90arcsec' / lake_in
    bf_300 = bf / '300arcsec' / lake_in
    bf_1800 = bf / '1800arcsec' / lake_in

    var_in = 'sfcWind'
    file_30 = bf_30 / Path(var_in + '_corr_v1.0_1979_2016.nc')
    file_90 = bf_90 / Path(var_in + '_corr_v1.0_1979_2016.nc')
    file_300 = bf_300 / Path(var_in + '_corr_v1.0_1979_2016.nc')
    file_1800 = bf_1800 / Path(var_in + '_corr_v1.0_1979_2016.nc')

    if  file_30.is_file()==False:
        print(var_in+'file 1km does not exist')
        wind_30 = xr.open_dataset(bf_ref / Path(var_in + '_corr_v1.0_1979_2016.nc'))
        wind_30['sfcWind'] = np.nan
        wind_30['time'] = np.nan
    else:
        wind_30 = xr.open_dataset(file_30)
        if np.sum(wind_30.sfcWind.isnull()).values > 0 :
            print('problem with 1km '+ var_in)
            os.remove(file_30)
        elif ((wind_30.sfcWind.values < 0) | (wind_30.sfcWind.values > 40)).any():
            print('out of bounds problem 1km ' + var_in)

    if  file_90.is_file()==False:
        print(var_in + 'file 3km does not exist')
        wind_90 = xr.open_dataset(bf_ref / Path(var_in + '_corr_v1.0_1979_2016.nc'))
        wind_90['sfcWind'] = np.nan
        wind_90['time'] = np.nan
    else:
        wind_90 = xr.open_dataset(bf_90 / Path(var_in + '_corr_v1.0_1979_2016.nc'))
        if np.sum(wind_90.sfcWind.isnull()).values > 0 :
            print('problem with 3km '+ var_in)
            os.remove(file_90)
        elif ((wind_90.sfcWind.values < 0) | (wind_90.sfcWind.values > 40)).any():
            print('out of bounds problem 3km ' + var_in)

    if  file_300.is_file()==False:
        print(var_in + 'file 10km does not exist')
        wind_300 = xr.open_dataset(bf_ref / Path(var_in + '_corr_v1.0_1979_2016.nc'))
        wind_300['sfcWind'] = np.nan
        wind_300['time'] = np.nan
    else:
        wind_300 = xr.open_dataset(bf_300 / Path(var_in + '_corr_v1.0_1979_2016.nc'))
        if np.sum(wind_300.sfcWind.isnull()).values > 0 :
            print('problem with 10km '+ var_in)
            os.remove(file_300)
        elif ((wind_300.sfcWind.values < 0) | (wind_300.sfcWind.values > 40)).any():
            print('out of bounds problem 10km ' + var_in)

    if  file_1800.is_file()==False:
        print(var_in + 'file 60km does not exist')
        wind_1800 = xr.open_dataset(bf_ref / Path(var_in + '_corr_v1.0_1979_2016.nc'))
        wind_1800['sfcWind'] = np.nan
        wind_1800['time'] = np.nan
    else:
        wind_1800 = xr.open_dataset(bf_1800 / Path(var_in + '_corr_v1.0_1979_2016.nc'))
        if np.sum(wind_1800.sfcWind.isnull()).values > 0 :
            print('problem with 60km '+ var_in)
            os.remove(file_1800)
        elif ((wind_1800.sfcWind.values < 0) | (wind_1800.sfcWind.values > 40)).any():
            print('out of bounds problem 60km ' + var_in)

    var_in = 'ps'
    file_30 = bf_30 / Path(var_in + '_corr_v1.0_1979_2016.nc')
    file_90 = bf_90 / Path(var_in + '_corr_v1.0_1979_2016.nc')
    file_300 = bf_300 / Path(var_in + '_corr_v1.0_1979_2016.nc')
    file_1800 = bf_1800 / Path(var_in + '_corr_v1.0_1979_2016.nc')

    if  file_30.is_file()==False:
        print(var_in + 'file 1km does not exist')
        ps_30 = xr.open_dataset(bf_ref / Path(var_in + '_corr_v1.0_1979_2016.nc'))
        ps_30['ps'] = np.nan
        ps_30['time'] = np.nan
    else:
        ps_30 = xr.open_dataset(bf_30 / Path(var_in + '_corr_v1.0_1979_2016.nc'))
        if np.sum(ps_30.ps.isnull()).values > 0:
            print('problem with 1km ' + var_in)
            os.remove(file_30)
        elif ((ps_30.ps.values < 50000) | (ps_30.ps.values > 107000)).any():
            print('out of bounds problem 1km ' + var_in)

    if  file_90.is_file()==False:
        print(var_in + 'file 3km does not exist')
        ps_90 = xr.open_dataset(bf_ref / Path(var_in + '_corr_v1.0_1979_2016.nc'))
        ps_90['ps'] = np.nan
        ps_90['time'] = np.nan
    else:
        ps_90 = xr.open_dataset(bf_90 / Path(var_in + '_corr_v1.0_1979_2016.nc'))
        if np.sum(ps_90.ps.isnull()).values > 0:
            print('problem with 3km ' + var_in)
            os.remove(file_90)
        elif ((ps_90.ps.values < 50000) | (ps_90.ps.values > 107000)).any():
            print('out of bounds problem 3km ' + var_in)

    if  file_300.is_file()==False:
        print(var_in + 'file 10km does not exist')
        ps_300 = xr.open_dataset(bf_ref / Path(var_in + '_corr_v1.0_1979_2016.nc'))
        ps_300['ps'] = np.nan
        ps_300['time'] = np.nan
    else:
        ps_300 = xr.open_dataset(bf_300 / Path(var_in + '_corr_v1.0_1979_2016.nc'))
        if np.sum(ps_300.ps.isnull()).values > 0:
            print('problem with 10km ' + var_in)
            os.remove(file_300)
        elif ((ps_300.ps.values < 50000) | (ps_300.ps.values > 107000)).any():
            print('out of bounds problem 10km ' + var_in)

    if  file_1800.is_file()==False:
        print(var_in + 'file 60km does not exist')
        ps_1800 = xr.open_dataset(bf_ref / Path(var_in + '_corr_v1.0_1979_2016.nc'))
        ps_1800['ps'] = np.nan
        ps_1800['time'] = np.nan
    else:
        ps_1800 = xr.open_dataset(bf_1800 / Path(var_in + '_corr_v1.0_1979_2016.nc'))
        if np.sum(ps_1800.ps.isnull()).values > 0:
            print('problem with 60km ' + var_in)
            os.remove(file_1800)
        elif ((ps_1800.ps.values < 50000) | (ps_1800.ps.values > 107000)).any():
            print('out of bounds problem 60km ' + var_in)
            print(np.min(ps_1800.ps.values))
            print(np.max(ps_1800.ps.values))

    var_in = 'rlds'
    file_30 = bf_30 / Path(var_in + '_corr_v1.0_1979_2016.nc')
    file_90 = bf_90 / Path(var_in + '_corr_v1.0_1979_2016.nc')
    file_300 = bf_300 / Path(var_in + '_corr_v1.0_1979_2016.nc')
    file_1800 = bf_1800 / Path(var_in + '_corr_v1.0_1979_2016.nc')

    if file_30.is_file()==False:
        print(var_in + 'file 1km does not exist')
        rlds_30 = xr.open_dataset(bf_ref / Path(var_in + '_corr_v1.0_1979_2016.nc'))
        rlds_30['rlds'] = np.nan
        rlds_30['time'] = np.nan
    else:
        rlds_30 = xr.open_dataset(bf_30 / Path(var_in + '_corr_v1.0_1979_2016.nc'))
        if np.sum(rlds_30.rlds.isnull()).values > 0:
            print('problem with 1km ' + var_in)
            os.remove(file_30)
        elif ((rlds_30.rlds.values < 100) | (rlds_30.rlds.values > 600)).any():
            print('out of bounds problem 1km ' + var_in)
            print(np.min(rlds_30.rlds.values))
            print(np.max(rlds_30.rlds.values))
            print(np.mean(rlds_30.rlds.values))

    if  file_90.is_file()==False:
        print(var_in + 'file 3km does not exist')
        rlds_90 = xr.open_dataset(bf_ref / Path(var_in + '_corr_v1.0_1979_2016.nc'))
        rlds_90['rlds'] = np.nan
        rlds_90['time'] = np.nan
    else:
        rlds_90 = xr.open_dataset(bf_90 / Path(var_in + '_corr_v1.0_1979_2016.nc'))
        if np.sum(rlds_90.rlds.isnull()).values > 0:
            print('problem with 3km ' + var_in)
            os.remove(file_90)
        elif ((rlds_90.rlds.values < 100) | (rlds_90.rlds.values > 600)).any():
            print('out of bounds problem 3km ' + var_in)
            print(np.min(rlds_90.rlds.values))
            print(np.max(rlds_90.rlds.values))

    if  file_300.is_file()==False:
        print(var_in + 'file 10km does not exist')
        rlds_300 = xr.open_dataset(bf_ref / Path(var_in + '_corr_v1.0_1979_2016.nc'))
        rlds_300['rlds'] = np.nan
        rlds_300['time'] = np.nan
    else:
        rlds_300 = xr.open_dataset(bf_300 / Path(var_in + '_corr_v1.0_1979_2016.nc'))
        if np.sum(rlds_300.rlds.isnull()).values > 0:
            print('problem with 10km ' + var_in)
            os.remove(file_300)
        elif ((rlds_300.rlds.values < 100) | (rlds_300.rlds.values > 600)).any():
            print('out of bounds problem 10km ' + var_in)
            print(np.min(rlds_300.rlds.values))
            print(np.max(rlds_300.rlds.values))

    if  file_1800.is_file()==False:
        print(var_in + 'file 60km does not exist')
        rlds_1800 = xr.open_dataset(bf_ref / Path(var_in + '_corr_v1.0_1979_2016.nc'))
        rlds_1800['rlds'] = np.nan
        rlds_1800['time'] = np.nan
    else:
        rlds_1800 = xr.open_dataset(bf_1800 / Path(var_in + '_corr_v1.0_1979_2016.nc'))
        if np.sum(rlds_1800.rlds.isnull()).values > 0:
            print('problem with 60km ' + var_in)
            os.remove(file_1800)
        elif ((rlds_1800.rlds.values < 100) | (rlds_1800.rlds.values > 600)).any():
            print('out of bounds problem 60km ' + var_in)
            print(np.min(rlds_1800.rlds.values))
            print(np.max(rlds_1800.rlds.values))
            print(np.mean(rlds_1800.rlds.values))

    var_in = 'rh'
    file_30 = bf_30 / Path(var_in + '_corr_v1.0_1979_2016.nc')
    file_90 = bf_90 / Path(var_in + '_corr_v1.0_1979_2016.nc')
    file_300 = bf_300 / Path(var_in + '_corr_v1.0_1979_2016.nc')
    file_1800 = bf_1800 / Path(var_in + '_corr_v1.0_1979_2016.nc')

    if  file_30.is_file()==False:
        print(var_in + 'file 1km does not exist')
        rh_30 = xr.open_dataset(bf_ref / Path(var_in + '_corr_v1.0_1979_2016.nc'))
        rh_30['hurs'] = np.nan
        rh_30['time'] = np.nan
    else:
        rh_30 = xr.open_dataset(bf_30 / Path(var_in + '_corr_v1.0_1979_2016.nc'))
        if np.sum(rh_30.hurs.isnull()).values > 0:
            print('problem with 1km ' + var_in)
            os.remove(file_30)
        elif ((rh_30.hurs.values < 0) | (rh_30.hurs.values > 100)).any():
            print('out of bounds problem 1km ' + var_in)
        elif np.max(rh_30.hurs.values) < 50:
            print('out of bounds problem 1km ' + var_in)

    if  file_90.is_file()==False:
        print(var_in + 'file 3km does not exist')
        rh_90 = xr.open_dataset(bf_ref / Path(var_in + '_corr_v1.0_1979_2016.nc'))
        rh_90['hurs'] = np.nan
        rh_90['time'] = np.nan
    else:
        rh_90 = xr.open_dataset(bf_90 / Path(var_in + '_corr_v1.0_1979_2016.nc'))
        if np.sum(rh_90.hurs.isnull()).values > 0:
            print('problem with 3km ' + var_in)
            os.remove(file_90)
        elif ((rh_90.hurs.values < 0) | (rh_90.hurs.values > 100)).any():
            print('out of bounds problem 3km ' + var_in)
        elif np.max(rh_90.hurs.values) < 50:
            print('out of bounds problem 3km ' + var_in)

    if  file_300.is_file()==False:
        print(var_in + 'file 10km does not exist')
        rh_300 = xr.open_dataset(bf_ref / Path(var_in + '_corr_v1.0_1979_2016.nc'))
        rh_300['hurs'] = np.nan
        rh_300['time'] = np.nan
    else:
        rh_300 = xr.open_dataset(bf_300 / Path(var_in + '_corr_v1.0_1979_2016.nc'))
        if np.sum(rh_300.hurs.isnull()).values > 0:
            print('problem with 10km ' + var_in)
            os.remove(file_300)
        elif ((rh_300.hurs.values < 0) | (rh_300.hurs.values > 100)).any():
            print('out of bounds problem 10km ' + var_in)
        elif np.max(rh_300.hurs.values) < 50:
            print('out of bounds problem 10km ' + var_in)

    if  file_1800.is_file()==False:
        print(var_in + 'file 60km does not exist')
        rh_1800 = xr.open_dataset(bf_ref / Path(var_in + '_corr_v1.0_1979_2016.nc'))
        rh_1800['hurs'] = np.nan
        rh_1800['time'] = np.nan
    else:
        rh_1800 = xr.open_dataset(bf_1800 / Path(var_in + '_corr_v1.0_1979_2016.nc'))
        if np.sum(rh_1800.hurs.isnull()).values > 0:
            print('problem with 60km ' + var_in)
            os.remove(file_1800)
        elif ((rh_1800.hurs.values < 0) | (rh_1800.hurs.values > 100)).any():
            print('out of bounds problem 60km ' + var_in)
        elif np.max(rh_1800.hurs.values) < 50:
            print('out of bounds problem 60km ' + var_in)

    rh_all = rh_30.rename({'hurs':'1km'})
    rh_all['3km'] = rh_90.hurs
    rh_all['10km'] = rh_300.hurs
    rh_all['60km'] = rh_1800.hurs
    rh_all_df = rh_all.to_dataframe()
    rh_all_df_1 = pd.melt(rh_all_df[['60km', '10km', '3km', '1km']], var_name="resolution", value_name='rh')
    rh_all_df_1['lake'] = lake_in

    lw_all = rlds_30.rename({'rlds':'1km'})
    lw_all['3km'] = rlds_90.rlds
    lw_all['10km'] = rlds_300.rlds
    lw_all['60km'] = rlds_1800.rlds
    lw_all_df = lw_all.to_dataframe()
    lw_all_df_1 = pd.melt(lw_all_df[['60km', '10km', '3km', '1km']], var_name="resolution", value_name='lwr')
    lw_all_df_1['lake'] = lake_in

    ps_all = ps_30.rename({'ps':'1km'})
    ps_all['3km'] = ps_90.ps
    ps_all['10km'] = ps_300.ps
    ps_all['60km'] = ps_1800.ps
    ps_all_df = ps_all.to_dataframe()
    ps_all_df_1 = pd.melt(ps_all_df[['60km', '10km', '3km', '1km']], var_name="resolution", value_name='ps')
    ps_all_df_1['lake'] = lake_in

    wind_all = wind_30.rename({'sfcWind':'1km'})
    wind_all['3km'] = wind_90.sfcWind
    wind_all['10km'] = wind_300.sfcWind
    wind_all['60km'] = wind_1800.sfcWind
    wind_all_df = wind_all.to_dataframe()
    wind_all_df_1 = pd.melt(wind_all_df[['60km', '10km', '3km', '1km']], var_name="resolution", value_name='wind')
    wind_all_df_1['lake'] = lake_in
    my_pal = {"1km": 'orange', "3km": (161 / 255, 218 / 255, 180 / 255),
              "10km": (65 / 255, 182 / 255, 196 / 255), "60km": (34 / 255, 94 / 255, 168 / 255)}
    flierprops = dict(marker='o', markersize=4.5, markeredgecolor='gray', markerfacecolor='silver', alpha=0.45)

    fig, axs = plt.subplots(4, sharex=True, sharey=False, figsize=(9, 12))
    sns.set_style("whitegrid")
    sns.boxplot(x="resolution", y='lwr', palette=my_pal, data=lw_all_df_1, flierprops=flierprops, showmeans=True,
                meanprops={"marker":"o", "markerfacecolor":"white", "markeredgecolor":"black", "markersize":"4"},
                medianprops=dict(color="grey", alpha=0.85, linewidth=1.9, linestyle='-'), linewidth=1.6, saturation=0.9, showfliers = True, ax=axs[0])
    axs[0].set_ylabel('Longwave [W m-2]')
    sns.boxplot(x="resolution", y='rh', palette=my_pal, data=rh_all_df_1, flierprops=flierprops, showmeans=True,
                meanprops={"marker": "o", "markerfacecolor": "white", "markeredgecolor": "black", "markersize": "4"},
                medianprops=dict(color="grey", alpha=0.85, linewidth=1.9, linestyle='-'), linewidth=1.6, saturation=0.9,
                showfliers=True, ax=axs[1])
    axs[1].set_ylabel('Relative Humidity [%]')
    sns.boxplot(x="resolution", y='ps', palette=my_pal, data=ps_all_df_1, flierprops=flierprops, showmeans=True,
                meanprops={"marker": "o", "markerfacecolor": "white", "markeredgecolor": "black", "markersize": "4"},
                medianprops=dict(color="grey", alpha=0.85, linewidth=1.9, linestyle='-'), linewidth=1.6, saturation=0.9,
                showfliers=True, ax=axs[2])
    axs[2].set_ylabel('Surface Air Pressure [Pa]')
    sns.boxplot(x="resolution", y='wind', palette=my_pal, data=wind_all_df_1, flierprops=flierprops, showmeans=True,
                meanprops={"marker": "o", "markerfacecolor": "white", "markeredgecolor": "black", "markersize": "4"},
                medianprops=dict(color="grey", alpha=0.85, linewidth=1.9, linestyle='-'), linewidth=1.6, saturation=0.9,
                showfliers=True, ax=axs[3])
    axs[3].set_ylabel('Wind Speed [m s-1]')
    axs[0].set_xlabel("")
    axs[1].set_xlabel("")
    axs[2].set_xlabel("")
    axs[3].set_xlabel("")
    fig.supxlabel('Resolution', fontsize=13)
    fig.suptitle(lake_in, fontsize=13)
    plt.tight_layout()
    (bf_out_figures / 'boxplots').mkdir(parents=True, exist_ok=True)
    fig.savefig(bf_out_figures / 'boxplots' / Path(lake_in + '_boxplot.png'), facecolor='white', transparent=False)

    fig, axs = plt.subplots(4, sharex=True, sharey=False, figsize=(8, 8))
    fig.supxlabel('Time', fontsize=13)
    fig.suptitle(lake_in, fontsize=13)
    line_w = 1.3
    axs[0].grid(True)
    axs[0].plot(rlds_30.time, rlds_30['rlds'].squeeze(), color=col_res[0], linewidth=line_w, label='1km')
    axs[0].plot(rlds_90.time, rlds_90['rlds'].squeeze(), color=col_res[1], linewidth=line_w, label='3km')
    axs[0].plot(rlds_300.time, rlds_300['rlds'].squeeze(), color=col_res[2], linewidth=line_w, label='10km')
    axs[0].plot(rlds_1800.time, rlds_1800['rlds'].squeeze(), color=col_res[3], linewidth=line_w, label='60km')
    axs[0].set_ylabel('Longwave [W m-2]')
    handles, labels = axs[0].get_legend_handles_labels()
    fig.legend(handles, labels, loc="upper left", bbox_to_anchor=(0.19, 0.96), fancybox=True, shadow=False,
               ncol=4)
    #axs[0].get_legend().remove()
    axs[0].set_xlabel("")

    axs[1].grid(True)
    axs[1].plot(rh_30.time, rh_30['hurs'].squeeze(), color=col_res[0], linewidth=line_w, label='1km')
    axs[1].plot(rh_90.time, rh_90['hurs'].squeeze(), color=col_res[1], linewidth=line_w, label='3km')
    axs[1].plot(rh_300.time, rh_300['hurs'].squeeze(), color=col_res[2], linewidth=line_w, label='10km')
    axs[1].plot(rh_1800.time, rh_1800['hurs'].squeeze(), color=col_res[3], linewidth=line_w, label='60km')
    axs[1].set_ylabel('Relative Humidity [%]')
    #axs[1].get_legend().remove()
    axs[1].set_xlabel("")

    axs[2].grid(True)
    axs[2].plot(ps_30.time, ps_30['ps'].squeeze(), color=col_res[0], linewidth=line_w, label='1km')
    axs[2].plot(ps_90.time, ps_90['ps'].squeeze(), color=col_res[1], linewidth=line_w, label='3km')
    axs[2].plot(ps_300.time, ps_300['ps'].squeeze(), color=col_res[2], linewidth=line_w, label='10km')
    axs[2].plot(ps_1800.time, ps_1800['ps'].squeeze(), color=col_res[3], linewidth=line_w, label='60km')
    axs[2].set_ylabel('Surface Air Pressure [Pa]')
    #axs[2].get_legend().remove()
    axs[2].set_xlabel("")

    axs[3].grid(True)
    axs[3].plot(wind_30.time, wind_30['sfcWind'].squeeze(), color=col_res[0], linewidth=line_w, label='1km')
    axs[3].plot(wind_90.time, wind_90['sfcWind'].squeeze(), color=col_res[1], linewidth=line_w, label='3km')
    axs[3].plot(wind_300.time, wind_300['sfcWind'].squeeze(), color=col_res[2], linewidth=line_w, label='10km')
    axs[3].plot(wind_1800.time, wind_1800['sfcWind'].squeeze(), color=col_res[3], linewidth=line_w, label='60km')
    axs[3].set_ylabel('Wind Speed [m s-1]')
    #axs[3].get_legend().remove()
    axs[3].set_xlabel("")

    axs[0].set_xlim(pd.Timestamp('2010-03-01'), pd.Timestamp('2010-09-01'))
    axs[1].set_xlim(pd.Timestamp('2010-03-01'), pd.Timestamp('2010-09-01'))
    axs[2].set_xlim(pd.Timestamp('2010-03-01'), pd.Timestamp('2010-09-01'))
    axs[3].set_xlim(pd.Timestamp('2010-03-01'), pd.Timestamp('2010-09-01'))

    plt.tight_layout()
    #plt.subplots_adjust(hspace=0.3)
    (bf_out_figures / 'time_series').mkdir(parents=True, exist_ok=True)
    fig.savefig(bf_out_figures / 'time_series' / Path(lake_in + '_comp_2010.png'), facecolor='white', transparent=False)

