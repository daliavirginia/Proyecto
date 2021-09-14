#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Aug 16 00:01:29 2021

@author: dalia
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

# Para graficar
LAT_SUR_G=-36
LAT_NOR_G=-25
LON_OESTE_G=293
LON_ESTE_G=300

# Para calculos
LAT_SUR=-39
LAT_NOR=-22
LON_ESTE=303
LON_OESTE=290


# Periodo 1
PER1_ANIO_MIN=1930
PER1_ANIO_MAX=2020

TIEMPO1_MIN = str(PER1_ANIO_MIN)+'-01-01'
TIEMPO1_MAX= str(PER1_ANIO_MAX)+'-12-31'

# Periodo 2
PER2_ANIO_MIN=1980
PER2_ANIO_MAX=2020

TIEMPO2_MIN=str(PER2_ANIO_MIN)+'-01-01'
TIEMPO2_MAX= str(PER2_ANIO_MAX)+'-12-31'

# Directorio al archivo
DATOS_GISS='/home/dalia/Proyecto/BasesDatos/air.2x2.250.mon.anom.land.nc'
DATOS_CRU='/home/dalia/Proyecto/BasesDatos/air.mon.anom.nc'
# Directorio a salidas
SALIDAS='/home/dalia/Proyecto/Salidas/'

#%%

# Abro el archivo usando xarray
data=xr.open_dataset(DATOS_GISS)

# Hago recorte para el periodo 1
ds_per1 = data.sel(lat=slice(LAT_NOR, LAT_SUR), lon=slice(LON_OESTE,LON_ESTE),
                        time=slice(TIEMPO1_MIN,TIEMPO1_MAX))
# Hago recorte para el período 2
ds_per2 = data.sel(lat=slice(LAT_NOR, LAT_SUR), lon=slice(LON_OESTE,LON_ESTE),
                        time=slice(TIEMPO2_MIN,TIEMPO2_MAX))


#%%

#Remuestreo los datos dividiendo por estaciones (DEF-MAM-JJA-SON)
ds_per1_est = ds_per1.resample(time="QS-DEC").mean()
ds_per2_est = ds_per2.resample(time="QS-DEC").mean()

#Divido en cuartos (no me gusta esta forma de separar, muchas variables
#capaz se puede mejorar)
ds_djf_per1=ds_per1_est.sel(time=ds_per1_est["time.month"]==12)
ds_mam_per1=ds_per1_est.sel(time=ds_per1_est["time.month"]==3)
ds_jja_per1=ds_per1_est.sel(time=ds_per1_est["time.month"]==6)
ds_son_per1=ds_per1_est.sel(time=ds_per1_est["time.month"]==9)

ds_djf_per2=ds_per2_est.sel(time=ds_per2_est["time.month"]==12)
ds_mam_per2=ds_per2_est.sel(time=ds_per2_est["time.month"]==3)
ds_jja_per2=ds_per2_est.sel(time=ds_per2_est["time.month"]==6)
ds_son_per2=ds_per2_est.sel(time=ds_per2_est["time.month"]==9)

#De los datos de verano me quedo desde el segundo hasta el anteultimo ya que
#se agregaron meses vacios para poder completar el primer y ultimo cuarto.
ds_djf_per1=ds_djf_per1.sel(time=slice(ds_djf_per1["time"][1],ds_djf_per1["time"][-2]))
ds_djf_per2=ds_djf_per2.sel(time=slice(ds_djf_per2["time"][1],ds_djf_per2["time"][-2]))

#%%
#Guardo los datos en un diccionario

ds_est1= {"Verano (DJF)":ds_djf_per1,
          "Otoño (MAM)":ds_mam_per1,
          "Invierno (JJA)":ds_jja_per1,
          "Primavera (SON)":ds_son_per1}

ds_est2= {"Verano (DJF)":ds_djf_per2,
          "Otoño (MAM)":ds_mam_per2,
          "Invierno (JJA)":ds_jja_per2,
          "Primavera (SON)":ds_son_per2}

#%%
#Calculo la tendencia por estacion usando la función antes definida

#Defino listas vacias dinámicas
tends_per1=[]
tends_per2=[]
pv_per1=[]
pv_per2=[]

#Con un ciclo for itero dentro de los diccionarios y calculo 
#la tendencia para cada estacion. Los resultados se guardan en las 
#listas definidas arriba
for key in ds_est1.keys():
    tend,pv=trend_xarray(ds_est1[key],'air')
    tends_per1.append(tend)
    pv_per1.append(pv)
    
for key in ds_est2.keys():
    tend,pv=trend_xarray(ds_est2[key],'air')
    tends_per2.append(tend)
    pv_per2.append(pv)

#%%

#Transformo las listas a diccionarios usanndo las keys de ds_est1 

tends_per1 = dict(zip(ds_est1.keys(), tends_per1))
tends_per2 = dict(zip(ds_est2.keys(), tends_per2))
pv_per1 = dict(zip(ds_est2.keys(), pv_per1))
pv_per2 = dict(zip(ds_est2.keys(), pv_per2))
#%%

#Grafico  período 1

#Define figure

#Tamaño para mapas del centro de Argentina
fig, ax = plt.subplots(figsize=(2*2.7,2*4))

#Tamaño para sur de sudamerica
#fig, ax = plt.subplots(figsize=(2*3.3,2*3.5))

#Define grid for subplots
gs = gridspec.GridSpec(2,2)     

#latitudes and longitudes to plot
lons, lats = np.meshgrid(ds_per1['lon'], ds_per1['lat'])

#Definimos los niveles para los contornos (inferior, superior, longitud de los intervalos)
#Si utilizamos las funciones ds.min ds.max podemos darnos una idea de cómo definirlos
lev_inf = -0.4
lev_sup = 0.4
lev_int =0.05
clevs = np.arange(lev_inf, lev_sup+lev_int, lev_int)

i=0
for key in tends_per1.keys():
    
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
    
    im=ax.contourf(lons, lats, tends_per1[key]*10,cmap="coolwarm", levels=clevs,extend='both',transform=crs_latlon)
    ax.contour(lons, lats, pv_per1[key],levels=[0.1],colors='r',linewidths=0.5 , transform=crs_latlon) 
    
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
cbar_ax.set_title('°C/decada')
fig.suptitle("Tendencia lineal de temperatura\n" + str(PER1_ANIO_MIN)+"-"+str(PER1_ANIO_MAX),fontsize=14)

fig.savefig(SALIDAS+"TendenciaTemp_Per1_GISS.png", dpi=300, bbox_inches='tight')
#%%

#Grafico período 2

#Define figure

#Tamaño para mapas del centro de Argentina
fig, ax = plt.subplots(figsize=(2*2.7,2*4))

#Tamaño para sur de sudamerica
#fig, ax = plt.subplots(figsize=(2*3.3,2*3.5))

#Define grid for subplots
gs = gridspec.GridSpec(2,2)     

#latitudes and longitudes to plot
lons, lats = np.meshgrid(ds_per2['lon'], ds_per2['lat'])

#Definimos los niveles para los contornos (inferior, superior, longitud de los intervalos)
#Si utilizamos las funciones ds.min ds.max podemos darnos una idea de cómo definirlos
lev_inf = -0.4
lev_sup = 0.4
lev_int =0.05
clevs = np.arange(lev_inf, lev_sup+lev_int, lev_int)

i=0
for key in tends_per1.keys():
    
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
    
    im=ax.contourf(lons, lats, tends_per2[key]*10,cmap="coolwarm", levels=clevs,extend='both',transform=crs_latlon)
    ax.contour(lons, lats, pv_per2[key],levels=[0.1],colors='r',linewidths=0.5 , transform=crs_latlon) 
    
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
cbar_ax.set_title('°C/decada')
fig.suptitle("Tendencia lineal de temperatura\n" + str(PER2_ANIO_MIN)+"-"+str(PER2_ANIO_MAX),fontsize=14)

fig.savefig(SALIDAS+"TendenciaTemp_Per2_GISS.png", dpi=300, bbox_inches='tight')

#%%

# CRU #

#%%

LON_ESTE_G=-60
LON_OESTE_G=-67

LON_ESTE=-57
LON_OESTE=-70

#%%

# Abro el archivo usando xarray
data=xr.open_dataset(DATOS_CRU)

# Hago recorte para el periodo 1
ds_per1 = data.sel(lat=slice(LAT_NOR, LAT_SUR), lon=slice(LON_OESTE,LON_ESTE),
                        time=slice(TIEMPO1_MIN,TIEMPO1_MAX))
# Hago recorte para el período 2
ds_per2 = data.sel(lat=slice(LAT_NOR, LAT_SUR), lon=slice(LON_OESTE,LON_ESTE),
                        time=slice(TIEMPO2_MIN,TIEMPO2_MAX))


#%%

#Remuestreo los datos dividiendo por estaciones (DEF-MAM-JJA-SON)
ds_per1_est = ds_per1.resample(time="QS-DEC").mean()
ds_per2_est = ds_per2.resample(time="QS-DEC").mean()

#Divido en cuartos (no me gusta esta forma de separar, muchas variables
#capaz se puede mejorar)
ds_djf_per1=ds_per1_est.sel(time=ds_per1_est["time.month"]==12)
ds_mam_per1=ds_per1_est.sel(time=ds_per1_est["time.month"]==3)
ds_jja_per1=ds_per1_est.sel(time=ds_per1_est["time.month"]==6)
ds_son_per1=ds_per1_est.sel(time=ds_per1_est["time.month"]==9)

ds_djf_per2=ds_per2_est.sel(time=ds_per2_est["time.month"]==12)
ds_mam_per2=ds_per2_est.sel(time=ds_per2_est["time.month"]==3)
ds_jja_per2=ds_per2_est.sel(time=ds_per2_est["time.month"]==6)
ds_son_per2=ds_per2_est.sel(time=ds_per2_est["time.month"]==9)

#De los datos de verano me quedo desde el segundo hasta el anteultimo ya que
#se agregaron meses vacios para poder completar el primer y ultimo cuarto.
ds_djf_per1=ds_djf_per1.sel(time=slice(ds_djf_per1["time"][1],ds_djf_per1["time"][-2]))
ds_djf_per2=ds_djf_per2.sel(time=slice(ds_djf_per2["time"][1],ds_djf_per2["time"][-2]))

#%%
#Guardo los datos en un diccionario

ds_est1= {"Verano (DJF)":ds_djf_per1,
          "Otoño (MAM)":ds_mam_per1,
          "Invierno (JJA)":ds_jja_per1,
          "Primavera (SON)":ds_son_per1}

ds_est2= {"Verano (DJF)":ds_djf_per2,
          "Otoño (MAM)":ds_mam_per2,
          "Invierno (JJA)":ds_jja_per2,
          "Primavera (SON)":ds_son_per2}

#%%
#Calculo la tendencia por estacion usando la función antes definida

#Defino listas vacias dinámicas
tends_per1=[]
tends_per2=[]
pv_per1=[]
pv_per2=[]

#Con un ciclo for itero dentro de los diccionarios y calculo 
#la tendencia para cada estacion. Los resultados se guardan en las 
#listas definidas arriba
for key in ds_est1.keys():
    tend,pv=trend_xarray(ds_est1[key],'air')
    tends_per1.append(tend)
    pv_per1.append(pv)
    
for key in ds_est2.keys():
    tend,pv=trend_xarray(ds_est2[key],'air')
    tends_per2.append(tend)
    pv_per2.append(pv)

#%%

#Transformo las listas a diccionarios usanndo las keys de ds_est1 

tends_per1 = dict(zip(ds_est1.keys(), tends_per1))
tends_per2 = dict(zip(ds_est2.keys(), tends_per2))
pv_per1 = dict(zip(ds_est2.keys(), pv_per1))
pv_per2 = dict(zip(ds_est2.keys(), pv_per2))
#%%

#Grafico  período 1

#Define figure

#Tamaño para mapas del centro de Argentina
fig, ax = plt.subplots(figsize=(2*2.7,2*4))


#Tamaño para sur de sudamerica
#fig, ax = plt.subplots(figsize=(2*3.3,2*3.5))

#Define grid for subplots
gs = gridspec.GridSpec(2,2)     

#latitudes and longitudes to plot
lons, lats = np.meshgrid(ds_per1['lon'], ds_per1['lat'])

#Definimos los niveles para los contornos (inferior, superior, longitud de los intervalos)
#Si utilizamos las funciones ds.min ds.max podemos darnos una idea de cómo definirlos
lev_inf = -0.4
lev_sup = 0.4
lev_int =0.05
clevs = np.arange(lev_inf, lev_sup+lev_int, lev_int)

i=0
for key in tends_per1.keys():
    
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
    
    im=ax.contourf(lons, lats, tends_per1[key]*10,cmap="coolwarm", levels=clevs,extend='both',transform=crs_latlon)
    ax.contour(lons, lats, pv_per1[key],levels=[0.1],colors='r',linewidths=0.5 , transform=crs_latlon) 
    
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
cbar_ax.set_title('°C/decada')
fig.suptitle("Tendencia lineal de temperatura\n" + str(PER1_ANIO_MIN)+"-"+str(PER1_ANIO_MAX),fontsize=14)

fig.savefig(SALIDAS+"TendenciaTemp_Per1_CRU.png", dpi=300, bbox_inches='tight')
#%%

#Grafico período 2

#Define figure

#Tamaño para mapas del centro de Argentina
fig, ax = plt.subplots(figsize=(2*2.7,2*4))

#Tamaño para sur de sudamerica
#fig, ax = plt.subplots(figsize=(2*3.3,2*3.5))

#Define grid for subplots
gs = gridspec.GridSpec(2,2)     

#latitudes and longitudes to plot
lons, lats = np.meshgrid(ds_per2['lon'], ds_per2['lat'])

#Definimos los niveles para los contornos (inferior, superior, longitud de los intervalos)
#Si utilizamos las funciones ds.min ds.max podemos darnos una idea de cómo definirlos
lev_inf = -0.4
lev_sup = 0.4
lev_int =0.05
clevs = np.arange(lev_inf, lev_sup+lev_int, lev_int)

i=0
for key in tends_per1.keys():
    
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
    
    im=ax.contourf(lons, lats, tends_per2[key]*10,cmap="coolwarm", levels=clevs,extend='both',transform=crs_latlon)
    ax.contour(lons, lats, pv_per2[key],levels=[0.1],colors='r',linewidths=0.5 , transform=crs_latlon) 
    
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
cbar_ax.set_title('°C/decada')
fig.suptitle("Tendencia lineal de temperatura\n" + str(PER2_ANIO_MIN)+"-"+str(PER2_ANIO_MAX),fontsize=14)

fig.savefig(SALIDAS+"TendenciaTemp_Per2_CRU.png", dpi=300, bbox_inches='tight')
