#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Aug 14 19:08:28 2021

@author: dalia
"""

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

def corr_xarray(ds,prom,VAR) :

    #Compute trends for each gridpoint
    ds_stack=ds.stack(points=['lat', 'lon']) 
    corr = np.empty_like(ds_stack[VAR][0,:]) # Arrays vacios
    pval = np.empty_like(ds_stack[VAR][0,:])
    for k in range(ds_stack[VAR].shape[1]):
        y = ds_stack[VAR][:, k]
        [corr[k], pval[k]] = stats.pearsonr(prom[VAR], y)       

    # Save trends of all ensemble members
    #Change to mm/day/summer/decade
    corr = np.reshape(corr, (len(ds['lat']), len(ds['lon'])))
    pval = np.reshape(pval, (len(ds['lat']), len(ds['lon'])))
    return corr,pval

#%%

# Para serie temporal CENTRO ARGENTINA
LAT_SUR_M=-36
LAT_NOR_M=-25
LON_OESTE_M=293
LON_ESTE_M=300

#Para calculos de correlación
LAT_SUR=-42
LAT_NOR=-16
LON_OESTE=287
LON_ESTE=306

#Para graficar 
LAT_SUR_G=-39
LAT_NOR_G=-19
LON_OESTE_G=290
LON_ESTE_G=303


ANIO_MIN=1920
ANIO_MAX=2020

TIMEMIN = str(ANIO_MIN)+'-01-01'
TIMEMAX = str(ANIO_MAX)+'-12-31'

# Directorio al archivo
# Directorio al archivo
DATOS_GISS='/home/dalia/Proyecto/BasesDatos/air.2x2.250.mon.anom.land.nc'
DATOS_CRU='/home/dalia/Proyecto/BasesDatos/air.mon.anom.nc'

# Directorio a salidas
SALIDAS = '/home/dalia/Proyecto/Salidas/'

#%%
    
# Abro el archivo usando xarray
ds=xr.open_dataset(DATOS_GISS)

# Anomalía de temperatura en kelvin respecto de la media 1951-1980


# Hago recorte
ds_recorte = ds.sel(lat=slice(LAT_NOR, LAT_SUR), lon=slice(LON_OESTE,LON_ESTE),
                        time=slice(TIMEMIN,TIMEMAX))

#%%

#Remuestreo los datos dividiendo por estaciones (DEF-MAM-JJA-SON)
ds_est = ds_recorte.resample(time="QS-DEC").mean()


#Divido en cuartos (no me gusta esta forma de separar, muchas variables
#capaz se puede mejorar)
ds_djf=ds_est.sel(time=ds_est["time.month"]==12)
ds_mam=ds_est.sel(time=ds_est["time.month"]==3)
ds_jja=ds_est.sel(time=ds_est["time.month"]==6)
ds_son=ds_est.sel(time=ds_est["time.month"]==9)

#De los datos de verano me quedo desde el segundo hasta el anteultimo ya que
#se agregaron meses vacios para poder completar el primer y ultimo cuarto.
ds_djf=ds_djf.sel(time=slice(ds_djf["time"][1],ds_djf["time"][-2]))


#%%
#Guardo los datos en un diccionario

ds_est_anom = {"Verano (DJF)":ds_djf,
               "Otoño (MAM)":ds_mam,
               "Invierno (JJA)":ds_jja,
               "Primavera (SON)":ds_son}


#%%

#Promedio regional con pesos segun latitud
ds_PromReg=[]
weights = np.cos(np.deg2rad(ds_recorte.lat))

for key in ds_est_anom.keys():
    ds_weighted = ds_est_anom[key].sel(lat=slice(LAT_NOR_M,LAT_SUR_M),
                                      lon=slice(LON_OESTE_M,LON_ESTE_M)).weighted(weights)
    ds_PromReg.append(ds_weighted.mean(("lon", "lat")))

                      
#Paso de list a dict
ds_PromReg= dict(zip(ds_est_anom.keys(), ds_PromReg))  

#Al calcular la media en cada punto de grilla debo considerar de que el area 
#en el cual se aplica el promedio es considerablemente mas grande en el Ecuador
#que en los Polos. Por eso, puedo aplicar una corrección, multiplicando la 
#media por el coseno de la latitud. Así se tendrá en cuenta que a menores latitudes
#el espacio es mayor y por ende el valor de pp pesa mas.


#%%

fig, ax = plt.subplots(figsize=(2*5.5,2*4.8))

gs = gridspec.GridSpec(2,2)

i=0
for key in ds_PromReg.keys():
    
    ax = plt.subplot(gs[i])
    bars = ax.bar(ds_PromReg[key]["air"]["time.year"],
                  ds_PromReg[key]["air"])
    
    for j in range(len(bars)):
    
        if bars[j].get_height()>0:
            bars[j].set_fc("orangered")
        else:
            bars[j].set_fc("lightseagreen")
            
    plt.ylim(-3,3)
    ax.set_ylabel("°C")
    ax.set_xlabel("Año")
    
    ax.set_title(str(key))
    
    i=i+1

fig.suptitle("Anomalía de temperatura", fontsize=18)

fig.savefig(SALIDAS+"AnomaliaTemp_CentroArgentina_GISS.png", dpi=300,  bbox_inches='tight')
plt.show()


#%%

#CORRELACION para cada punto de reticula con el promedio regional


#%%
corr=[]
pval=[]

for key in ds_PromReg.keys():
    cor,pv=corr_xarray(ds_est_anom[key], ds_PromReg[key], 'air')
    corr.append(cor)
    pval.append(pv)

corr= dict(zip(ds_est.keys(), corr))  
pval= dict(zip(ds_est.keys(), pval))  

#%%

#Tamaño para mapas del centro de Argentina
fig, ax = plt.subplots(figsize=(2*2.7,2*3.6))

#Define grid for subplots
gs = gridspec.GridSpec(2,2)     

#latitudes and longitudes to plot
lons, lats = np.meshgrid(ds_recorte['lon'], ds_recorte['lat'])

#Definimos los niveles para los contornos (inferior, superior, longitud de los intervalos)
#Si utilizamos las funciones ds.min ds.max podemos darnos una idea de cómo definirlos
lev_inf = -1
lev_sup = 1
lev_int =0.1
clevs = np.arange(lev_inf, lev_sup+lev_int, lev_int)

i=0
for key in corr.keys():
    
    ax = plt.subplot(gs[i],projection=ccrs.PlateCarree(central_longitude=180))
    crs_latlon = ccrs.PlateCarree()
    ax.set_extent([LON_OESTE_G, LON_ESTE_G, LAT_SUR_G, LAT_NOR_G], crs=crs_latlon)
    ax.add_feature(cartopy.feature.OCEAN, zorder=100, edgecolor='k')
    ax.add_feature(cartopy.feature.COASTLINE)
    ax.add_feature(cartopy.feature.BORDERS, linestyle='-', alpha=.5)
    # Create a feature for States/Admin 1 regions at 1:50m from Natural Earth
    states_provinces = cartopy.feature.NaturalEarthFeature(
        category='cultural',
        name='admin_1_states_provinces_lines',
        scale='10m',
        facecolor='none')
    
    ax.add_feature(states_provinces, edgecolor='darkslategrey', linewidths=0.3)
    
    im=ax.contourf(lons, lats, corr[key],cmap="coolwarm", levels=clevs,extend='both',transform=crs_latlon)
    ax.contour(lons, lats, pval[key],levels=[0.1],colors='k',linewidths=0.5 , transform=crs_latlon) 
    
    ax.set_xticks(np.arange(LON_OESTE_G,LON_ESTE_G,10), crs=crs_latlon)
    ax.set_yticks(np.arange(LAT_SUR_G,LAT_NOR_G,10), crs=crs_latlon)
    ax.grid(which='both', linewidth=0.3, linestyle='-')
    ax.tick_params(axis='both', which='major', labelsize=8)
    lon_formatter = LongitudeFormatter(zero_direction_label=True)
    lat_formatter = LatitudeFormatter()
    ax.xaxis.set_major_formatter(lon_formatter)
    ax.yaxis.set_major_formatter(lat_formatter)    
    
    ax.set_title(str(key))
    
    i=i+1

cbar_ax = fig.add_axes([0.98, 0.1, 0.03, 0.8])

fig.colorbar(im, cax=cbar_ax,orientation='vertical') 
fig.suptitle("Correlación con promedio regional",fontsize=14)

fig.savefig(SALIDAS+"Correlacion_Temp_GISS.png", dpi=300, bbox_inches='tight')


#%%
# CRU #

#%%

# Para serie temporal CENTRO ARGENTINA
LON_OESTE_M=-67
LON_ESTE_M=-60

#Para calculos de correlación
LON_OESTE=-73
LON_ESTE=-54

#Para graficar 
LON_OESTE_G=-70
LON_ESTE_G=-57


#%%
    
# Abro el archivo usando xarray
ds=xr.open_dataset(DATOS_CRU)

# Anomalía de temperatura en kelvin respecto de la media 1951-1980


# Hago recorte
ds_recorte = ds.sel(lat=slice(LAT_NOR, LAT_SUR), lon=slice(LON_OESTE,LON_ESTE),
                        time=slice(TIMEMIN,TIMEMAX))



#%%

#Remuestreo los datos dividiendo por estaciones (DEF-MAM-JJA-SON)
ds_est = ds_recorte.resample(time="QS-DEC").mean()


#Divido en cuartos (no me gusta esta forma de separar, muchas variables
#capaz se puede mejorar)
ds_djf=ds_est.sel(time=ds_est["time.month"]==12)
ds_mam=ds_est.sel(time=ds_est["time.month"]==3)
ds_jja=ds_est.sel(time=ds_est["time.month"]==6)
ds_son=ds_est.sel(time=ds_est["time.month"]==9)

#De los datos de verano me quedo desde el segundo hasta el anteultimo ya que
#se agregaron meses vacios para poder completar el primer y ultimo cuarto.
ds_djf=ds_djf.sel(time=slice(ds_djf["time"][1],ds_djf["time"][-2]))


#%%
#Guardo los datos en un diccionario

ds_est_anom = {"Verano (DJF)":ds_djf,
               "Otoño (MAM)":ds_mam,
               "Invierno (JJA)":ds_jja,
               "Primavera (SON)":ds_son}


#%%

#Promedio regional con pesos segun latitud
ds_PromReg=[]
weights = np.cos(np.deg2rad(ds_recorte.lat))

for key in ds_est_anom.keys():
    ds_weighted = ds_est_anom[key].sel(lat=slice(LAT_NOR_M,LAT_SUR_M),
                                      lon=slice(LON_OESTE_M,LON_ESTE_M)).weighted(weights)
    ds_PromReg.append(ds_weighted.mean(("lon", "lat")))

                      
#Paso de list a dict
ds_PromReg= dict(zip(ds_est_anom.keys(), ds_PromReg))  

#Al calcular la media en cada punto de grilla debo considerar de que el area 
#en el cual se aplica el promedio es considerablemente mas grande en el Ecuador
#que en los Polos. Por eso, puedo aplicar una corrección, multiplicando la 
#media por el coseno de la latitud. Así se tendrá en cuenta que a menores latitudes
#el espacio es mayor y por ende el valor de pp pesa mas.


#%%

fig, ax = plt.subplots(figsize=(2*5.5,2*4.8))

gs = gridspec.GridSpec(2,2)

i=0
for key in ds_PromReg.keys():
    
    ax = plt.subplot(gs[i])
    bars = ax.bar(ds_PromReg[key]["air"]["time.year"],
                  ds_PromReg[key]["air"])
    
    for j in range(len(bars)):
    
        if bars[j].get_height()>0:
            bars[j].set_fc("orangered")
        else:
            bars[j].set_fc("lightseagreen")
            
    plt.ylim(-3,3)
    ax.set_ylabel("°C")
    ax.set_xlabel("Año")
    
    ax.set_title(str(key))
    
    i=i+1

fig.suptitle("Anomalía de temperatura", fontsize=18)

fig.savefig(SALIDAS+"AnomaliaTemp_CentroArgentina_CRU.png", dpi=300,  bbox_inches='tight')
plt.show()


#%%

#CORRELACION para cada punto de reticula con el promedio regional


#%%
corr=[]
pval=[]

for key in ds_PromReg.keys():
    cor,pv=corr_xarray(ds_est_anom[key], ds_PromReg[key], 'air')
    corr.append(cor)
    pval.append(pv)

corr= dict(zip(ds_est.keys(), corr))  
pval= dict(zip(ds_est.keys(), pval))  

#%%

#Tamaño para mapas del centro de Argentina
fig, ax = plt.subplots(figsize=(2*2.7,2*3.6))

#Define grid for subplots
gs = gridspec.GridSpec(2,2)     

#latitudes and longitudes to plot
lons, lats = np.meshgrid(ds_recorte['lon'], ds_recorte['lat'])

#Definimos los niveles para los contornos (inferior, superior, longitud de los intervalos)
#Si utilizamos las funciones ds.min ds.max podemos darnos una idea de cómo definirlos
lev_inf = -1
lev_sup = 1
lev_int =0.1
clevs = np.arange(lev_inf, lev_sup+lev_int, lev_int)

i=0
for key in corr.keys():
    
    ax = plt.subplot(gs[i],projection=ccrs.PlateCarree(central_longitude=180))
    crs_latlon = ccrs.PlateCarree()
    ax.set_extent([LON_OESTE_G, LON_ESTE_G, LAT_SUR_G, LAT_NOR_G], crs=crs_latlon)
    ax.add_feature(cartopy.feature.OCEAN, zorder=100, edgecolor='k')
    ax.add_feature(cartopy.feature.COASTLINE)
    ax.add_feature(cartopy.feature.BORDERS, linestyle='-', alpha=.5)
    # Create a feature for States/Admin 1 regions at 1:50m from Natural Earth
    states_provinces = cartopy.feature.NaturalEarthFeature(
        category='cultural',
        name='admin_1_states_provinces_lines',
        scale='10m',
        facecolor='none')
    
    ax.add_feature(states_provinces, edgecolor='darkslategrey', linewidths=0.3)
    
    im=ax.contourf(lons, lats, corr[key],cmap="coolwarm", levels=clevs,extend='both',transform=crs_latlon)
    ax.contour(lons, lats, pval[key],levels=[0.1],colors='k',linewidths=0.5 , transform=crs_latlon) 
    
    ax.set_xticks(np.arange(LON_OESTE_G,LON_ESTE_G,10), crs=crs_latlon)
    ax.set_yticks(np.arange(LAT_SUR_G,LAT_NOR_G,10), crs=crs_latlon)
    ax.grid(which='both', linewidth=0.3, linestyle='-')
    ax.tick_params(axis='both', which='major', labelsize=8)
    lon_formatter = LongitudeFormatter(zero_direction_label=True)
    lat_formatter = LatitudeFormatter()
    ax.xaxis.set_major_formatter(lon_formatter)
    ax.yaxis.set_major_formatter(lat_formatter)    
    
    ax.set_title(str(key))
    
    i=i+1

cbar_ax = fig.add_axes([0.98, 0.1, 0.03, 0.8])

fig.colorbar(im, cax=cbar_ax,orientation='vertical') 
fig.suptitle("Correlación con promedio regional",fontsize=14)

fig.savefig(SALIDAS+"Correlacion_Temp_CRU.png", dpi=300, bbox_inches='tight')