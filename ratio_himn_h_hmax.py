
"""
 Script for plotting Figures 4a,4b and 5b of González-Haro et al. Structural and dynamical quality assessment of gap-filled sea surface temperature products.
 Input data of singularity exponents of all SST datasets are available at https://digital.csic.es/handle/10261/309548
 proposed directory (dir_data)
------------------------------------------------------------------------
* Input parameters
dataset: Label for identifying dataset according to input folder (i.e: SST_AMSR2_REMSS,SST_CMC,OSTIA, SST_MUR, SST_CCI)
------------------------------------------------------------------------
* Help, how to run it
python plot_h.py "SST_CMC"
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
import cartopy.feature as cfeature
import glob
from datetime import datetime as dt
import numpy.ma as ma
import sys; sys.path.append('../')
#sys.path.append('/home/cgonzalez/CGH/tool_scripts/pyhton_script/pyseidon')




plt.rcParams.update({'font.size': 20})




dir_graph='../figures/agulhas/
dir_data='../data/'



#total arguments
n = len(sys.argv)
print("Total arguments passed:", n)

# Arguments passed
print("\nName of Python script:", sys.argv[0])


dataset = sys.argv[1]


in_files = sorted(glob.glob(dir_data+dataset+'/*.nc'))



ds=xr.open_mfdataset(in_files)

nt,nlat,nlon=ds.h.shape
print(nt,nlat,nlon)




mask_n_prob_h=np.zeros([nt,nlat,nlon],dtype='int')

#Definition of singularity exponents range
hmin=0.3
hmax=0.6


# In[48]:


#land mask index
inan=np.where(np.isnan(ds.h[0,:,:]) )


# In[24]:


for t in range(nt):
    h_dum=np.asarray(ds.h[t,:,:].values)
    mask_dum=np.zeros([nlat,nlon])
    mask_dum[inan]=np.nan#land pixels set to nan
    indx=((h_dum > hmin) & (h_dum < hmax))
    mask_dum[indx]=1
    mask_n_prob_h[t,:,:]=mask_dum
    



#Probability of having a singularity exponents within the range [hmin,hmax]
prob=mask_n_prob_h.sum(axis=0, dtype='float')/nt


#Land mask
prob[inan]=np.nan
ind=(prob == 0.)
prob[ind]=np.nan

proj_carr=ccrs.PlateCarree()#choose of projection
plt.rcParams.update({'font.size': 20})

fig = plt.figure(figsize=(10, 5),dpi=600)
vmin = 0
vmax = 0.2

levels = np.linspace(vmin, vmax, num=256)
ticks = np.linspace(vmin, vmax, num=11)
ax = plt.axes(projection=proj_carr)

cs= ax.pcolormesh(ds.lon,
                  ds.lat,
                  prob,
                  transform=proj_carr, cmap='Blues', vmin = vmin, vmax = vmax)
## with cmap attribute we select the color palette, more info: https://matplotlib.org/stable/gallery/color/colormap_reference.html

cbar = fig.colorbar(cs, orientation='vertical', ax=ax, label=r'fraction',ticks=ticks, pad=0.05)


ax.coastlines(resolution='50m', color='black', linewidth=1)
#for selecting the land mask
ax.add_feature(cfeature.LAND,facecolor=cfeature.COLORS['land'])
gl = ax.gridlines(draw_labels = True)
gl.top_labels =False
gl.ylabels_right =False
plt.title(r'prob '+str(hmin)+'$< h <$'+str(hmax)+' ('+dataset+')')
plt.savefig(dir_graph+'Ratio_hmin_'+str(hmin)+'h_hmax_'+str(hmax)+'_'+dataset+'.png')
print('Figure saved in ',dir_graph+'Ratio_hmin_h_hmax__'+dataset+'.png')



###Save to netCDF

lon=ds.lon
lat=ds.lat

now = dt.now()
dt_string = now.strftime("%d/%m/%Y %H:%M:%S")
filename='prob_h_min_'+str(hmin)+'max_'+str(hmax)+'_'+dataset
    


attrs_lat ={'standard_name':'latitude','long_name':'Latitude','units':'degrees_north','axis':'Y'}
attrs_lon ={'standard_name':'longitude','long_name':'Longitude','units':'degrees_east','axis':'X'}
attrs_h={'units':'1', 'long_name': 'Prob hmin<h<hmax'}
attrs_glo={'source_file':filename,'title':'Map of probability of hmin<h<hmax','file_created':dt_string,'hmin':hmin,'hmax':hmax}


ds_prob= xr.Dataset(
data_vars=dict(
    prob_h= ([ "lat", "lon"],  prob,attrs_h),
),
coords=dict(
    lat = (["lat"], lat,attrs_lat),
    lon = (["lon"], lon,attrs_lon),
),
attrs=attrs_glo
)



file_name = 'Prob_'+ filename+'.nc'
ds_prob.to_netcdf(dir_out+file_name,encoding = {"prob_h": {"dtype": "int16", "scale_factor": 0.01, "zlib": True}})
print('Saving file:',dir_out+file_name)








