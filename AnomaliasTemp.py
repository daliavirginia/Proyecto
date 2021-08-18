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

LAT_SUR=-36
LAT_NOR=-25
LON_OESTE=-67
LON_ESTE=-60

ANIO_MIN=1880
ANIO_MAX=2019

TIMEMIN = str(ANIO_MIN)+'-01-01'
TIMEMAX = str(ANIO_MAX)+'-12-31'

# Directorio al archivo
DATOS = '/home/dalia/Proyecto/BasesDatos/gistemp250_GHCNv4.nc'

# Directorio a salidas
SALIDAS = '/home/dalia/Proyecto/Salidas/'

#%%
    
# Abro el archivo usando xarray
ds=xr.open_dataset(DATOS)

print(ds["time"])
print(ds.var)
print(ds["tempanomaly"])
print(ds.attrs)
print(ds.lat)
print(ds.lon)
# Anomalía de temperatura en kelvin respecto de la media 1951-1980


# Hago recorte
ds_recorte = ds.sel(lat=slice(LAT_SUR, LAT_NOR), lon=slice(LON_OESTE,LON_ESTE),
                        time=slice(TIMEMIN,TIMEMAX))

print(ds_recorte)

ds_anual = ds_recorte.resample(time = "Y").mean()

#%%

### CORRECCIÓN DEL CALCULO ###

#Al calcular la media en cada punto de grilla debo considerar de que el area 
#en el cual se aplica el promedio es considerablemente mas grande en el Ecuador
#que en los Polos. Por eso, puedo aplicar una corrección, multiplicando la 
#media por el coseno de la latitud. Así se tendrá en cuenta que a menores latitudes
#el espacio es mayor y por ende el valor de pp pesa mas.

weights = np.cos(np.deg2rad(ds_anual.lat))

#Aplicamos los pesos para el dataset que estamos usando con la función "weighted"
ds_weighted = ds_anual.weighted(weights)
#Ahora sí, hacemos el promedio para las latitudes y longitudes pero considerando los pesos
#ds_PromEspacial=ds_weighted.mean(("lon", "lat"))

ds_PromReg=ds_anual.mean(("lon","lat"))
#%%

fig, ax = plt.subplots()

#ax.plot(ds_anual_anomalia["precip"]["time.year"], [0]*len(ds_anual_anomalia["precip"]["time.year"]), "k-")
bars = ax.bar(ds_PromReg["tempanomaly"]["time.year"],ds_PromReg["tempanomaly"])

for i in range(len(bars)):

    if bars[i].get_height()>0:
        bars[i].set_fc("darkturquoise")
    else:
        bars[i].set_fc("tomato")
        
plt.ylim(-2,2)
ax.set_ylabel("°C")


plt.title("Anomalía de temperatura para el centro de Argentina \n (25°S-36°S, 67°O-60°O)")

fig.savefig(SALIDAS+"AnomaliaTemp_CentroArgentina.png", dpi=300)
plt.show()
