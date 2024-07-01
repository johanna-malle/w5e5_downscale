# -*- coding: utf-8 -*-
"""
Desc:
Created on 10.05.23 13:48
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


bf = Path('/home/lud11/malle/Adhoris/forest_sector/caraib/forcing')
var_all = ['rh', 'sfcWind']

bf_30 = bf / 'res_30arcsec'
bf_90 = bf / 'res_90arcsec'
bf_300 = bf / 'res_300arcsec'
bf_1800 = bf / 'res_1800arcsec'

my_pal = {"1km": 'orange', "3km": (161 / 255, 218 / 255, 180 / 255),
          "10km": (65 / 255, 182 / 255, 196 / 255), "60km": (34 / 255, 94 / 255, 168 / 255)}

my_pal = {"1km": 'orange', "3km": (161 / 255, 218 / 255, 180 / 255),
          "10km": (65 / 255, 182 / 255, 196 / 255), "60km": (34 / 255, 94 / 255, 168 / 255)}

col_res = ['orange', (161 / 255, 218 / 255, 180 / 255),
           (65 / 255, 182 / 255, 196 / 255), (34 / 255, 94 / 255, 168 / 255)]

flierprops = dict(marker='o', markersize=4.5, markeredgecolor='gray', markerfacecolor='silver', alpha=0.45)

for var_in in var_all:
    file_30=xr.open_dataset(bf_30 / var_in / Path(var_in+'_corr_v1.0_2005.nc'))
    file_90=xr.open_dataset(bf_90 / var_in / Path(var_in+'_corr_v1.0_2005.nc'))
    file_300=xr.open_dataset(bf_300 / var_in / Path(var_in+'_corr_v1.0_2005.nc'))
    file_1800=xr.open_dataset(bf_1800 / var_in / Path(var_in+'_corr_v1.0_2005.nc'))

    if var_in =='rh':
        var_in='hurs'
    min_all=np.min([file_30[var_in].isel(time=0).min().values,file_90[var_in].isel(time=0).min().values,
                    file_300[var_in].isel(time=0).min().values,file_1800[var_in].isel(time=0).min().values])
    max_all= np.max([file_30[var_in].isel(time=0).max().values,file_90[var_in].isel(time=0).max().values,
                    file_300[var_in].isel(time=0).max().values,file_1800[var_in].isel(time=0).max().values])


    fig, axs = plt.subplots(2,2, sharex=False, sharey=False, figsize=(11, 9))
    file_30[var_in].isel(time=0).plot(ax=axs[0,0], cbar_kwargs={'label': ''},vmin=min_all,vmax=max_all)
    file_90[var_in].isel(time=0).plot(ax=axs[0,1],vmin=min_all,vmax=max_all)
    file_300[var_in].isel(time=0).plot(ax=axs[1,0], cbar_kwargs={'label': ''},vmin=min_all,vmax=max_all)
    file_1800[var_in].isel(time=0).plot(ax=axs[1,1],vmin=min_all,vmax=max_all)
    axs[0,0].set_title('30arcsec')
    axs[0,1].set_title('90arcsec')
    axs[0,0].set_xlabel('')
    axs[0,1].set_xlabel('')
    axs[0,1].set_ylabel('')
    axs[1,0].set_title('300arcsec')
    axs[1,1].set_title('1800arcsec')
    axs[1,1].set_ylabel('')
    fig.suptitle(file_30[var_in].time[0].values, fontsize=13)
    fig.savefig(bf /  Path(var_in + '_comp_20050101.png'), facecolor='white', transparent=False)


    #
    # data_all = file_30.rename({var_in: '1km'})
    # data_all['3km'] = file_90[var_in]
    # data_all['10km'] = file_300[var_in]
    # data_all['60km'] = file_1800[var_in]
    # data_all_df = data_all.to_dataframe()
    # data_all_df_1 = pd.melt(data_all_df[['60km', '10km', '3km', '1km']], var_name="resolution", value_name=var_in)
    #
    # fig, axs = plt.subplots(2, sharex=False, sharey=False, figsize=(9, 12))
    # #sns.set_style("whitegrid")
    # sns.boxplot(x="resolution", y=var_in, palette=my_pal, data=data_all_df_1, flierprops=flierprops, showmeans=True,
    #             meanprops={"marker": "o", "markerfacecolor": "white", "markeredgecolor": "black", "markersize": "4"},
    #             medianprops=dict(color="grey", alpha=0.85, linewidth=1.9, linestyle='-'), linewidth=1.6, saturation=0.9,
    #             showfliers=True, ax=axs[0])
    #
    # line_w = 1.3
    # # axs[1].grid(True)
    # # axs[1].plot(file_30.time, file_30[var_in].squeeze(), color=col_res[0], linewidth=line_w, label='1km')
    # # axs[1].plot(file_90.time, file_90[var_in].squeeze(), color=col_res[1], linewidth=line_w, label='3km')
    # # axs[1].plot(file_300.time, file_300[var_in].squeeze(), color=col_res[2], linewidth=line_w, label='10km')
    # # axs[1].plot(data_all.time, data_all[var_in].squeeze(), color=col_res[3], linewidth=line_w, label='60km')
    # plt.show()
