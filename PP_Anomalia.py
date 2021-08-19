#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Aug  8 18:43:27 2021

@author: dalia
"""
# PROMEDIO CENTRO DEL PAIS Y EVOLUCION TEMPORAL

#%%

#Load requested libraries
import numpy as np   # Manejar arrays
import xarray as xr  # Manejar arrays (especializado nc)
from scipy import stats # LLamo a la fun stats de scipy 

from matplotlib import pyplot as plt #Graficar (ploteos basicos)
from matplotlib import gridspec 
import matplotlib

#import cartopy.crs as ccrs	# Graficar mapas 
#import cartopy.feature 	# Notacion de punto (escribir todo)
#from cartopy.mpl.ticker import LongitudeFormatter, LatitudeFormatter

#%%

LAT_SUR=-36
LAT_NOR=-25
LON_OESTE=293
LON_ESTE=300

ANIO_MIN=1930
ANIO_MAX=2016

TIMEMIN = str(ANIO_MIN)+'-01-01'
TIMEMAX = str(ANIO_MAX)+'-12-31'

# Directorio al archivo
DATOS='/home/dalia/Proyecto/BasesDatos/precip.mon.total.v2018.nc'

# Directorio a salidas
SALIDAS='/home/dalia/Proyecto/Salidas/'

#%%

# Abro el archivo usando xarray
ds=xr.open_dataset(DATOS)

# Hago recorte
ds_recorte = ds.sel(lat=slice(LAT_NOR, LAT_SUR), lon=slice(LON_OESTE,LON_ESTE),
                        time=slice(TIMEMIN,TIMEMAX))

ds_anual = ds_recorte.resample(time="Y").sum()


ds_media = ds_anual.mean()

ds_anual_anomalia = ds_anual.mean("lat").mean("lon")-ds_media

#Al calcular la media en cada punto de grilla debo considerar de que el area 
#en el cual se aplica el promedio es considerablemente mas grande en el Ecuador
#que en los Polos. Por eso, puedo aplicar una corrección, multiplicando la 
#media por el coseno de la latitud. Así se tendrá en cuenta que a menores latitudes
#el espacio es mayor y por ende el valor de pp pesa mas.

weights = np.cos(np.deg2rad(ds_anual.lat))

#Aplicamos los pesos para el dataset que estamos usando con la función "weighted"
ds_anual_weighted = ds_anual.weighted(weights)
#Ahora sí, hacemos el promedio para las latitudes y longitudes pero considerando los pesos
ds_anual_PromReg=ds_anual_weighted.mean(("lon", "lat"))

ds_anom_corr = ds_anual_PromReg-ds_media

#%%

fig, ax = plt.subplots()

#ax.plot(ds_anual_anomalia["precip"]["time.year"], [0]*len(ds_anual_anomalia["precip"]["time.year"]), "k-")
bars = ax.bar(ds_anual_anomalia["precip"]["time.year"],ds_anual_anomalia["precip"])

for i in range(len(bars)):

    if bars[i].get_height()>0:
        bars[i].set_fc("darkturquoise")
    else:
        bars[i].set_fc("tomato")
        
plt.ylim(-400,400)
ax.set_ylabel("mm")
ax.set_xlabel("Año")

plt.title("Anomalía de precipitación para el centro de Argentina \n (25°S-36°S, 67°O-60°O)")

fig.savefig(SALIDAS+"AnomaliaPP_CentroArgentina.png", dpi=300)
plt.show()

#%%

fig, ax = plt.subplots()

#ax.plot(ds_anual_anomalia["precip"]["time.year"], [0]*len(ds_anual_anomalia["precip"]["time.year"]), "k-")
bars = ax.bar(ds_anom_corr["precip"]["time.year"],ds_anom_corr["precip"])

for i in range(len(bars)):

    if bars[i].get_height()>0:
        bars[i].set_fc("darkturquoise")
    else:
        bars[i].set_fc("tomato")
        
plt.ylim(-400,400)
ax.set_ylabel("mm")
ax.set_xlabel("Año")

plt.title("Anomalía de precipitación para el centro de Argentina \n (25°S-36°S, 67°O-60°O)")

fig.savefig(SALIDAS+"AnomaliaPP_Corregida_CentroArgentina.png", dpi=300)
plt.show()




