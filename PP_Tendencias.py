#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Aug  6 17:50:06 2021

@author: Dalia
"""

#%%
# HACER TENDENCIAS LINEALES PARA DOS PERÍODOS DISTINTOS: 1940-1970 Y 1980-2010
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

    #Los parametros de la función son:ds instancia de Dataset y VAR la variable
    #que quiero del dataset
    
    #Apilo latitudes y longitudes
    ds_stack=ds.stack(points=['lat', 'lon']) 
    
    #Arrays vacios para todos los tiempos y el primer punto de grilla
    trends = np.empty_like(ds_stack[VAR][0,:]) # Arrays vacios
    pval = np.empty_like(ds_stack[VAR][0,:])
    
    #Ciclo para completar los array vacios. Indice k recorre el rango de los
    #tiempos de ds_stack (shape[1] es para seleccionar la dim del tiempo)
    for k in range(ds_stack[VAR].shape[1]):
        y = ds_stack[VAR][:, k] #y es todos los puntos de grilla del tiempo k
        [trends[k], interc, r_va, pval[k], z] = stats.linregress(np.arange(len(ds_stack['time'])), y)    

    # Le regreso la forma de latxlon
    tend = np.reshape(trends, (len(ds['lat']), len(ds['lon']))) 
    pv_tend = np.reshape(pval, (len(ds['lat']), len(ds['lon'])))
    return tend,pv_tend

#%%

# Declaro variables

# Constantes (en mayuscula)
LAT_SUR=-60
LAT_NOR=-20
LON_ESTE=320
LON_OESTE=280

# Periodo 1
PER1_ANIO_MIN=1940
PER1_ANIO_MAX=1970

TIEMPO1_MIN = str(PER1_ANIO_MIN)+'-01-01'
TIEMPO1_MAX= str(PER1_ANIO_MAX)+'-12-31'

# Periodo 2
PER2_ANIO_MIN=1980
PER2_ANIO_MAX=2010

TIEMPO2_MIN=str(PER1_ANIO_MIN)+'-01-01'
TIEMPO2_MAX= str(PER2_ANIO_MAX)+'-12-31'

# Directorio al archivo
DATOS='/home/dalia/Proyecto/BasesDatos/precip.mon.total.v2018.nc'

# Directorio a salidas
SALIDAS='/home/dalia/Proyecto/Salidas/'

#%%

# Abro el archivo usando xarray
data=xr.open_dataset(DATOS)

# Hago recorte para el periodo 1
data_per1 = data.sel(lat=slice(LAT_NOR, LAT_SUR), lon=slice(LON_OESTE,LON_ESTE),
                        time=slice(TIEMPO1_MIN,TIEMPO1_MAX))
# Hago recorte para el período 2
data_per2 = data.sel(lat=slice(LAT_NOR, LAT_SUR), lon=slice(LON_OESTE,LON_ESTE),
                        time=slice(TIEMPO2_MIN,TIEMPO2_MAX))


#%%
#Calculo la tendencia usando la función antes definida
 
tend_per1,pv_tend_per1=trend_xarray(data_per1,'precip')
tend_per2,pv_tend_per2=trend_xarray(data_per2,'precip')

#%%

#Grafico 

#Define figure
fig, ax = plt.subplots(figsize=(5.3,3.3))

#Define grid for subplots
gs = gridspec.GridSpec(1,2)     

#latitudes and longitudes to plot
lons, lats = np.meshgrid(data_per1['lon'], data_per1['lat'])

#Tendencias periodo 1

ax = plt.subplot(gs[0],projection=ccrs.PlateCarree(central_longitude=180))

crs_latlon = ccrs.PlateCarree()
ax.set_extent([LON_OESTE, LON_ESTE, LAT_SUR, LAT_NOR], crs=crs_latlon)
ax.add_feature(cartopy.feature.COASTLINE)
ax.add_feature(cartopy.feature.BORDERS, linestyle='-', alpha=.5)

# Create a feature for States/Admin 1 regions at 1:50m from Natural Earth
states_provinces = cartopy.feature.NaturalEarthFeature(
    category='cultural',
    name='admin_1_states_provinces_lines',
    scale='10m',
    facecolor='none')

ax.add_feature(states_provinces, edgecolor='darkslategrey', linewidths=0.3)

im=ax.contourf(lons, lats, tend_per1*10,cmap="viridis_r",extend='both',transform=crs_latlon)
ax.contour(lons, lats, pv_tend_per1,levels=[0.1],colors='r',linewidths=0.5 , transform=crs_latlon) 

ax.set_xticks(np.arange(LON_OESTE,LON_ESTE,10), crs=crs_latlon)
ax.set_yticks(np.arange(LAT_SUR,LAT_NOR,10), crs=crs_latlon)
ax.grid(which='both', linewidth=0.3, linestyle='-')
ax.tick_params(axis='both', which='major', labelsize=5)
lon_formatter = LongitudeFormatter(zero_direction_label=True)
lat_formatter = LatitudeFormatter()
ax.xaxis.set_major_formatter(lon_formatter)
ax.yaxis.set_major_formatter(lat_formatter)    

ax.set_title(str(PER1_ANIO_MIN)+"-"+str(PER1_ANIO_MAX))


#Tendencias periodo 2

ax = plt.subplot(gs[1],projection=ccrs.PlateCarree(central_longitude=180))

crs_latlon = ccrs.PlateCarree()
ax.set_extent([LON_OESTE, LON_ESTE, LAT_SUR, LAT_NOR], crs=crs_latlon)
ax.add_feature(cartopy.feature.COASTLINE)
ax.add_feature(cartopy.feature.BORDERS, linestyle='-', alpha=.5)

# Create a feature for States/Admin 1 regions at 1:50m from Natural Earth
states_provinces = cartopy.feature.NaturalEarthFeature(
    category='cultural',
    name='admin_1_states_provinces_lines',
    scale='10m',
    facecolor='none')

ax.add_feature(states_provinces, edgecolor='darkslategrey', linewidths=0.3)

im=ax.contourf(lons, lats, tend_per2*10,cmap="viridis_r",extend='both',transform=crs_latlon)
ax.contour(lons, lats, pv_tend_per2,levels=[0.1],colors='r',linewidths=0.5 , transform=crs_latlon) 

ax.set_xticks(np.arange(LON_OESTE,LON_ESTE,10), crs=crs_latlon)
ax.set_yticks(np.arange(LAT_SUR,LAT_NOR,10), crs=crs_latlon)
ax.grid(which='both', linewidth=0.3, linestyle='-')
ax.tick_params(axis='both', which='major', labelsize=5)
lon_formatter = LongitudeFormatter(zero_direction_label=True)
lat_formatter = LatitudeFormatter()
ax.xaxis.set_major_formatter(lon_formatter)
ax.yaxis.set_major_formatter(lat_formatter)   

ax.set_title(str(PER2_ANIO_MIN)+"-"+str(PER2_ANIO_MAX))

#Plot colorbar
cbar_ax = fig.add_axes([0.98, 0.1, 0.03, 0.8])

fig.colorbar(im, cax=cbar_ax,orientation='vertical') 



fig.suptitle("Tendencia lineal de precipitación \n GPCC",fontsize=14)

fig.savefig(SALIDAS+"TendenciasLineales.png", dpi=300, bbox_inches='tight')

#%%

