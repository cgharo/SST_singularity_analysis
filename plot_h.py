#!/usr/bin/env python
# coding: utf-8



import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import numpy as np
import xarray as xr
import os as os
import pandas as pd
import cartopy.feature as cfeature
import glob
from datetime import datetime as dt
import numpy.ma as ma
import matplotlib.path as mpath
import matplotlib.ticker as mticker
import sys; sys.path.append('../')


dir_graph='../figures/agulhas/
dir_data='../data/'

#####################
#Script for plotting Figures 1 and 2 of González-Haro et al. Structural and dynamical quality assessment of gap-filled sea surface temperature products
#
# data of singularity exponents of all SST datasets are available at https://digital.csic.es/handle/10261/309548
# proposed directory 
# /data/
#	/SST_CMC
#	/SST_AMSR2_REMSS
#	/SST_CCI
#	/OSTIA
#	SST_MUR
###################


def plot_h(dataset,dir_data,file, dpi=None,region='Agulhas'):
    infile=dir_data+dataset+'/'+file

    ds=xr.open_dataset(infile)
    if(region=='Agulhas'):
        latmin=-45
        latmax=-30
        lonmin=10
        lonmax=35
        ds=ds.select(lat=slice(latmin,latmax),lon=slice(lonmin,lonmax))
        img_extent=[lonmin,lonmax,latmin,latmax]

    h_dum = np.asarray(ds.h)
    print(np.max(h_dum), np.min(h_dum))
    h_dum[(h_dum > 5)] = np.nan

    proj_carr=ccrs.PlateCarree()#choose of projection
    plt.rcParams.update({'font.size': 20})
    print('Plotting h')
    fig = plt.figure(figsize=(10, 5),dpi=dpi)
    vmin = -0.5
    vmax = 0.5

    levels = np.linspace(vmin, vmax, num=256)
    ticks = np.linspace(vmin, vmax, num=11)
    ax = plt.axes(projection=proj_carr)

    cs= ax.pcolormesh(ds.lon,
                        ds.lat,
                        h_dum[0,:,:],
                        transform=proj_carr, cmap='Greys', vmin = vmin, vmax = vmax)
    ## with cmap attribute we select the color palette, more info: https://matplotlib.org/stable/gallery/color/colormap_reference.html

    cbar = fig.colorbar(cs, orientation='vertical', ax=ax, label=r'h',ticks=ticks, pad=0.05)


    ax.coastlines(resolution='50m', color='black', linewidth=1)
    #for selecting the land mask
    ax.add_feature(cfeature.LAND, facecolor = cfeature.COLORS['land'])
    #ax.stock_img()
    ax.set_extent(img_extent, crs=proj_carr)
    gl = ax.gridlines(draw_labels = True)
    gl.top_labels =False
    gl.ylabels_right =False

    plt.title(dataset)
    #plt.tight_layout()
    plt.subplots_adjust(left=0.1, right=0.95, top=0.9, bottom=0.1)
    str_time='2016-09-29'
    plt.savefig(dir_graph+'sigularity_exponents_'+dataset+'_'+str_time+'.png',facecolor='white',dpi=dpi)
    
    print('Figure saved in ',dir_graph+'sigularity_exponents_'+dataset+'_'+str_time+'.png')
    return None


dataset='SST_AMSR2_REMSS'
file='Singularity_exponentsRSS_AMSR2_ocean_L3_3day_2016-09-29_v08.nc'
plot_h(dataset,dir_data,file, dpi=600)

#dataset='SST_CMC'
#file='Singularity_exponents20160929120000-CMC-L4_GHRSST-SSTfnd-CMC0.nc'
#plot_h(dataset,dir_data,file,dpi=600)

dataset='OSTIA'
file='Singularity_exponents20160929120000-UKMO-L4_GHRSST-SSTfnd-OSTIA-GLOB-v02.nc'
plot_h(dataset,dir_data,file,dpi=600)

dataset='SST_CCI'
file='Singularity_exponents20160929120000-ESACCI-L4_GHRSST-SSTdepth-OSTIA-GLOB_CDR2.nc' 
plot_h(dataset,dir_data,file,dpi=600)

dataset='SST_MUR'
file='Singularity_exponents20160702090000-JPL-L4_GHRSST-SSTfnd-MUR-GLOB-v02.nc'
plot_h(dataset,dir_data,file,dpi=600)
