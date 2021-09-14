#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Dec  2 17:35:26 2019

@author: ldiaz

"""

#%%

#Plot precipitation trends


#%%

#Load requested libraries
import numpy as np   # Manejar arrays
import xarray as xr  # Manejar arrays (especializado nc)
from scipy import stats # LLamo a la fun stats de scipy 

from matplotlib import pyplot as plt #Graficar (ploteos basicos)
from matplotlib import gridspec 
import matplotlib

import cartopy.crs as ccrs	# Graficar mapas 
import cartopy.feature 	# Notacion de punto (escribir todo)
from cartopy.mpl.ticker import LongitudeFormatter, LatitudeFormatter

#%%

#Definefunction to compute trend for xarray (with lon and lat as dimensions)

def trend_xarray(ds,VAR) :

    #Compute trends for each gridpoint
    ds_stack=ds.stack(points=['lat', 'lon']) 
    trends = np.empty_like(ds_stack[VAR][0,:]) # Arrays vacios
    pval = np.empty_like(ds_stack[VAR][0,:])
    for k in range(ds_stack[VAR].shape[1]):
        y = ds_stack[VAR][:, k]
        [trends[k], interc, r_va, pval[k], z] = stats.linregress(np.arange(len(ds_stack['time'])), y)       

    # Save trends of all ensemble members
    #Change to mm/day/summer/decade
    tend = np.reshape(trends, (len(ds['lat']), len(ds['lon'])))
    pv_tend = np.reshape(pval, (len(ds['lat']), len(ds['lon'])))
    return tend,pv_tend

#%%
#Define path to save figure

FIG_PATH = '/home/dalia/Proyecto/Salidas'

#%%

#Region to plot. 
LATMIN_PLOT = -36
LATMAX_PLOT = -25
LONMIN_PLOT = 293
LONMAX_PLOT = 300

#Region to calculate.
LATMIN = LATMIN_PLOT-3
LATMAX = LATMAX_PLOT+3
LONMIN = LONMIN_PLOT-3
LONMAX = LONMAX_PLOT+3

#%%

# 1950-2014 Available period
YEARMIN=1901
YEARMAX=2005

TIMEMIN = str(YEARMIN)+'-01-01'
TIMEMAX = str(YEARMAX)+'-12-31'

#%%

#Path to observational datasets

FILE_OBS = '/home/dalia/Proyecto/BasesDatos/precip.mon.total.v2018.nc'

#%%

#Compute seasonal means

#Open observed data with xarray
ds_obs = xr.open_dataset(FILE_OBS, decode_coords=False) #data set
#Select latitudinal, longitudinal and temporal range 
ds_reg_obs = ds_obs.sel(lat=slice(LATMAX, LATMIN), lon=slice(LONMIN,LONMAX), time=slice(TIMEMIN, TIMEMAX))
#Compute seasonal accumulated (3-month) (DJF season) (from mm/season to mm/day)
ds_reg_3m_obs = ds_reg_obs.resample(time='QS-DEC').sum(skipna=False)         


#Select summer accumulated (DJF season)
ds_djf = ds_reg_3m_obs.sel(time=ds_reg_3m_obs['time.month']==12)
#Exclude JF of first year and december of the last
ds_djf=ds_djf.sel(time=slice(ds_djf['time'][1],ds_djf['time'][-2]))
#Select autumn accumulated (MAM season)
ds_mam = ds_reg_3m_obs.sel(time=ds_reg_3m_obs['time.month']==3)   
#Select winter accumulated (JJA season)
ds_jja = ds_reg_3m_obs.sel(time=ds_reg_3m_obs['time.month']==6)   
#Select spring accumulated (SON season)
ds_son = ds_reg_3m_obs.sel(time=ds_reg_3m_obs['time.month']==9)   

#%%

#Compute trends for each season

tend_djf,pv_tend_djf=trend_xarray(ds_djf,'precip')
tend_mam,pv_tend_mam=trend_xarray(ds_mam,'precip')
tend_jja,pv_tend_jja=trend_xarray(ds_jja,'precip')
tend_son,pv_tend_son=trend_xarray(ds_son,'precip')

#Trend for point nearest to Las Bre?as
#No funca
LBaproxtrend=tend_djf[np.argmin(np.abs(ds_djf['lat']-(-27.1))),np.argmin(np.abs(ds_djf['lon']-(360-61.1)))]*10

print(LBaproxtrend)

#%%

#Define figure
fig, ax = plt.subplots(figsize=(2*3.3,2*5.3))

#Define grid for subplots
gs = gridspec.GridSpec(2,2)     

#Define levels and colors
levels = [-12,-9,-6,-3,0,3,6,9,12]
colors=['#67001f','#b2182b','#d6604d','#f4a582','#fddbc7','#d1e5f0','#92c5de','#4393c3','#2166ac','#053061']
cmap, norm = matplotlib.colors.from_levels_and_colors(levels, colors,extend='both')

#latitudes and longitudes to plot
lons, lats = np.meshgrid(ds_reg_obs['lon'], ds_reg_obs['lat'])

#Plot DJF trends

ax = plt.subplot(gs[0],projection=ccrs.PlateCarree(central_longitude=180))

crs_latlon = ccrs.PlateCarree()
ax.set_extent([LONMIN_PLOT, LONMAX_PLOT, LATMIN_PLOT, LATMAX_PLOT], crs=crs_latlon)
ax.add_feature(cartopy.feature.COASTLINE)
ax.add_feature(cartopy.feature.BORDERS, linestyle='-', alpha=.5)

# Create a feature for States/Admin 1 regions at 1:50m from Natural Earth
states_provinces = cartopy.feature.NaturalEarthFeature(
    category='cultural',
    name='admin_1_states_provinces_lines',
    scale='10m',
    facecolor='none')

ax.add_feature(states_provinces, edgecolor='gray')

im=ax.contourf(lons, lats, tend_djf*10,levels=levels,cmap=cmap,norm=norm,extend='both',transform=crs_latlon)
ax.contour(lons, lats, pv_tend_djf,levels=[0.1],colors='r',linewidths=0.5 , transform=crs_latlon) 
#plt.colorbar(im,orientation='horizontal')

ax.set_xticks([295, 300], crs=crs_latlon)
ax.set_yticks([-35,-30, -25], crs=crs_latlon)
ax.grid(which='both', linewidth=0.3, linestyle='-')
ax.tick_params(axis='both', which='major', labelsize=10)
lon_formatter = LongitudeFormatter(zero_direction_label=True)
lat_formatter = LatitudeFormatter()
ax.xaxis.set_major_formatter(lon_formatter)
ax.yaxis.set_major_formatter(lat_formatter)    

# Add meteorological stations location and study territories

ax.plot(360-62,-26,color='green', marker='s',markersize=5,transform=crs_latlon,label='Copo')
ax.plot(360-65.4,-31.33,color='orange', marker='s',markersize=5,transform=crs_latlon,label='Chancani')

ax.plot(360-61.1,-27.1,color='blue', marker='o',markersize=3,transform=crs_latlon,label='INTA - Las Breñas (EMC)')
ax.plot(360-(65+10/100),-(31+57/60),color='red', marker='o',markersize=3,transform=crs_latlon,label='SMN - Villa Dolores AERO')

ax.set_title('Verano (DEF)')

#Plot MAM trends

ax = plt.subplot(gs[1],projection=ccrs.PlateCarree(central_longitude=180))

crs_latlon = ccrs.PlateCarree()
ax.set_extent([LONMIN_PLOT, LONMAX_PLOT, LATMIN_PLOT, LATMAX_PLOT], crs=crs_latlon)
ax.add_feature(cartopy.feature.COASTLINE)
ax.add_feature(cartopy.feature.BORDERS, linestyle='-', alpha=.5)

# Create a feature for States/Admin 1 regions at 1:50m from Natural Earth
states_provinces = cartopy.feature.NaturalEarthFeature(
    category='cultural',
    name='admin_1_states_provinces_lines',
    scale='10m',
    facecolor='none')

ax.add_feature(states_provinces, edgecolor='gray')

im=ax.contourf(lons, lats, tend_mam*10,levels=levels,cmap=cmap,norm=norm,extend='both',transform=crs_latlon)
ax.contour(lons, lats, pv_tend_mam,levels=[0.1],colors='r',linewidths=0.5 , transform=crs_latlon) 
#plt.colorbar(im,orientation='horizontal')

ax.set_xticks([295, 300], crs=crs_latlon)
ax.set_yticks([-35,-30, -25], crs=crs_latlon)
ax.grid(which='both', linewidth=0.3, linestyle='-')
ax.tick_params(axis='both', which='major', labelsize=10)
lon_formatter = LongitudeFormatter(zero_direction_label=True)
lat_formatter = LatitudeFormatter()
ax.xaxis.set_major_formatter(lon_formatter)
ax.yaxis.set_major_formatter(lat_formatter)    

# Add meteorological stations location and study territories

ax.plot(360-62,-26,color='green', marker='s',markersize=5,transform=crs_latlon,label='Copo')
ax.plot(360-65.4,-31.33,color='orange', marker='s',markersize=5,transform=crs_latlon,label='Chancani')

ax.plot(360-61.1,-27.1,color='blue', marker='o',markersize=3,transform=crs_latlon,label='INTA - Las Breñas (EMC)')
ax.plot(360-(65+10/100),-(31+57/60),color='red', marker='o',markersize=3,transform=crs_latlon,label='SMN - Villa Dolores AERO')

ax.set_title('Otoño (MAM)')

#Plot MAM trends

ax = plt.subplot(gs[2],projection=ccrs.PlateCarree(central_longitude=180))

crs_latlon = ccrs.PlateCarree()
ax.set_extent([LONMIN_PLOT, LONMAX_PLOT, LATMIN_PLOT, LATMAX_PLOT], crs=crs_latlon)
ax.add_feature(cartopy.feature.COASTLINE)
ax.add_feature(cartopy.feature.BORDERS, linestyle='-', alpha=.5)

# Create a feature for States/Admin 1 regions at 1:50m from Natural Earth
states_provinces = cartopy.feature.NaturalEarthFeature(
    category='cultural',
    name='admin_1_states_provinces_lines',
    scale='10m',
    facecolor='none')

ax.add_feature(states_provinces, edgecolor='gray')

im=ax.contourf(lons, lats, tend_jja*10,levels=levels,cmap=cmap,norm=norm,extend='both',transform=crs_latlon)
ax.contour(lons, lats, pv_tend_jja,levels=[0.1],colors='r',linewidths=0.5 , transform=crs_latlon) 
#plt.colorbar(im,orientation='horizontal')

ax.set_xticks([295, 300], crs=crs_latlon)
ax.set_yticks([-35,-30, -25], crs=crs_latlon)
ax.grid(which='both', linewidth=0.3, linestyle='-')
ax.tick_params(axis='both', which='major', labelsize=10)
lon_formatter = LongitudeFormatter(zero_direction_label=True)
lat_formatter = LatitudeFormatter()
ax.xaxis.set_major_formatter(lon_formatter)
ax.yaxis.set_major_formatter(lat_formatter)    

# Add meteorological stations location and study territories

ax.plot(360-62,-26,color='green', marker='s',markersize=5,transform=crs_latlon,label='Copo')
ax.plot(360-65.4,-31.33,color='orange', marker='s',markersize=5,transform=crs_latlon,label='Chancani')

ax.plot(360-61.1,-27.1,color='blue', marker='o',markersize=3,transform=crs_latlon,label='INTA - Las Breñas (EMC)')
ax.plot(360-(65+10/100),-(31+57/60),color='red', marker='o',markersize=3,transform=crs_latlon,label='SMN - Villa Dolores AERO')

ax.set_title('Invierno (JJA)')

#Plot MAM trends

ax = plt.subplot(gs[3],projection=ccrs.PlateCarree(central_longitude=180))

crs_latlon = ccrs.PlateCarree()
ax.set_extent([LONMIN_PLOT, LONMAX_PLOT, LATMIN_PLOT, LATMAX_PLOT], crs=crs_latlon)
ax.add_feature(cartopy.feature.COASTLINE)
ax.add_feature(cartopy.feature.BORDERS, linestyle='-', alpha=.5)

# Create a feature for States/Admin 1 regions at 1:50m from Natural Earth
states_provinces = cartopy.feature.NaturalEarthFeature(
    category='cultural',
    name='admin_1_states_provinces_lines',
    scale='10m',
    facecolor='none')

ax.add_feature(states_provinces, edgecolor='gray')

im=ax.contourf(lons, lats, tend_son*10,levels=levels,cmap=cmap,norm=norm,extend='both',transform=crs_latlon)
ax.contour(lons, lats, pv_tend_son,levels=[0.1],colors='r',linewidths=0.5 , transform=crs_latlon) 
#plt.colorbar(im,orientation='horizontal')

ax.set_xticks([295, 300], crs=crs_latlon)
ax.set_yticks([-35,-30, -25], crs=crs_latlon)
ax.grid(which='both', linewidth=0.3, linestyle='-')
ax.tick_params(axis='both', which='major', labelsize=10)
lon_formatter = LongitudeFormatter(zero_direction_label=True)
lat_formatter = LatitudeFormatter()
ax.xaxis.set_major_formatter(lon_formatter)
ax.yaxis.set_major_formatter(lat_formatter)    

# Add meteorological stations location and study territories

ax.plot(360-62,-26,color='green', marker='s',markersize=5,transform=crs_latlon,label='Copo')
ax.plot(360-65.4,-31.33,color='orange', marker='s',markersize=5,transform=crs_latlon,label='Chancani')

ax.plot(360-61.1,-27.1,color='blue', marker='o',markersize=3,transform=crs_latlon,label='INTA - Las Breñas (EMC)')
ax.plot(360-(65+10/100),-(31+57/60),color='red', marker='o',markersize=3,transform=crs_latlon,label='SMN - Villa Dolores AERO')

ax.set_title('Primavera (SON)')

# Ajustar espacio que dejas en los costados, arriba y abajo, y el espacio entre los subplots
fig.subplots_adjust(left=0.02, bottom=0.03, right=0.93,top=0.93, wspace=0.15 ,hspace=0.2)

fig.suptitle("Tendencia precipitación GPCC "+str(YEARMIN)+"-"+str(YEARMAX)+" [mm/estación/decada]",fontsize=14,y=0.98)

#Plot colorbar
cbar_ax = fig.add_axes([0.96, 0.2, 0.03, 0.6])
fig.colorbar(im, cax=cbar_ax,orientation='vertical')


fig.savefig(FIG_PATH + 'TendPR_'+str(YEARMIN)+'_'+str(YEARMAX)+'.png', dpi=300, bbox_inches='tight')

