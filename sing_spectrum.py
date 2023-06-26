"""
 Script for plotting Figures 4a,4b and 5b of González-Haro et al. Structural and dynamical quality assessment of gap-filled sea surface temperature products.
 Input data of singularity exponents of all SST datasets are available at https://digital.csic.es/handle/10261/309548
 proposed directory (dir_data)
------------------------------------------------------------------------
* Input parameters
dataset: Label for identifying dataset according to input folder (i.e: SST_AMSR2_REMSS,SST_CMC,OSTIA, SST_MUR, SST_CCI)
------------------------------------------------------------------------
* Help, how to run it
python 
------------------------------------
* /data/
	/SST_CMC
	/SST_AMSR2_REMSS
	/SST_CCI
	/OSTIA
	/SST_MUR

written by Cristina González-Haro
"""



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
import multifractal as mfr
# Each time a module is modified it is releaded



plt.rcParams.update({'font.size': 20})

#dataset='SST_CCI'

bec13= True

# total arguments
n = len(sys.argv)
print("Total arguments passed:", n)
 
# Arguments passed
print("\nName of Python script:", sys.argv[0])
 

dataset = sys.argv[1]
dir_data=sys.argv[2]
dir_graph=sys.argv[3]
print(dataset)

#if(bec13):
    #dir_data='/mnt/lustre/users/cgonzalez/singularity/outputs/resample_python/'
    #dir_graph='/mnt/lustre/users/cgonzalez/singularity/figures/spectrum/resample_python/'
#else:
#    dir_data='/home/cgonzalez/CGH/Projects/SST_GHRSST_Analysis/singularity/outputs/'
#    dir_data='/home/cgonzalez/CGH/Projects/SST_GHRSST_Analysis/singularity/figures/spectrum/'

print(dir_data+dataset+'/*/*.nc')
in_files = sorted(glob.glob(dir_data+dataset+'/*/*.nc'))

ds=xr.open_mfdataset(in_files)
#nt,nlat,nlon=ds.h.shape
nt,nlat,nlon=ds.exponents.shape
hmin=-0.5
hmax=1.5
#hbin=0.02 #initial value
hbin = 0.08

if (dataset=='SST_AMSR2_REMSS'):
    t=0
    #h=np.asarray(ds.h[t,:,:].data)
    h=np.asarray(ds.exponents[t,:,:].data)
    h[h>5]=np.nan
    #modified
    h = ma.masked_values(h, np.nan)
    h[h.mask==True]=-999
    h.fill_value=-999  

    # Compute singualrity spectra

    hx, hist, dh, edh = mfr.singspec(h[h.mask==False], hmin=hmin, hmax=hmax, hbin=hbin, dmax=2)
    shift = mfr.transinv(h, dmax=2, hbin=hbin,hmin=hmin, hmax=hmax,  return_all=False)
    nbin_h=hx.shape[0]
    spect_hx_amsre=np.zeros((nt,nbin_h))
    shift_h=np.zeros(nt)
    spect_hx_amsre.shape
else:
    
    ds_med=xr.open_dataset('/mnt/lustre/users/cgonzalez/singularity/outputs/SST_AMSR2_REMSS/median/median_spect_sing_SST_AMSR2.nc')
    hx_median=np.asarray(ds_med.median_spect.data)

fig = plt.figure(figsize=(10, 5))

for t in range(nt):
   #modified to check resampled netcdf
   #h=np.asarray(ds.h[t,:,:].data)
   h=np.asarray(ds.exponents[t,:,:].data)
   h[h>5]=np.nan
   #modified
   h = ma.masked_values(h, np.nan)
   h[h.mask==True]=-999
   h.fill_value=-999  

   # Compute singualrity spectra

   hx, hist, dh, edh = mfr.singspec(h[h.mask==False], hmin=hmin, hmax=hmax, hbin=hbin, dmax=2)
   shift = mfr.transinv(h, dmax=2, hbin=hbin,hmin=hmin, hmax=hmax,  return_all=False)

   plt.plot(hx+shift,dh,c='gray')
   if (dataset=='SST_AMSR2_REMSS'):
       spect_hx_amsre[t,:]=dh
       shift_h[t]=shift

if (dataset=='SST_AMSR2_REMSS'):
    hx_median=np.median(spect_hx_amsre,axis=0)
plt.plot(hx,hx_median,color='red')
plt.grid()
plt.xlabel('h')
plt.ylabel('D(h)')
plt.title(dataset)
plt.ylim([0,2.05])
plt.xlim([hmin,hmax])
plt.tight_layout()
plt.savefig(dir_graph+'Dh_singularity_spectrum_'+dataset+'.png')
print('Figure saved')
