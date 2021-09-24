#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Sep 21 17:49:22 2021

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
def corr_xarray(ds,prom) :

    #Compute trends for each gridpoint
    ds_stack=ds.stack(points=['lat', 'lon']) 
    corr = np.empty_like(ds_stack[0,:]) # Arrays vacios
    pval = np.empty_like(ds_stack[0,:])
    for k in range(ds_stack.shape[1]):
        y = ds_stack[:, k]
        [corr[k], pval[k]] = stats.pearsonr(prom, y)       

    # Save trends of all ensemble members
    #Change to mm/day/summer/decade
    corr = np.reshape(corr, (len(ds['lat']), len(ds['lon'])))
    pval = np.reshape(pval, (len(ds['lat']), len(ds['lon'])))
    return corr,pval

#%%
# Valores constantes

# Para serie temporal de anomalias CENTRO ARGENTINA
LAT_SUR_M=-35
LAT_NOR_M=-26
LON_OESTE_M=294
LON_ESTE_M=299

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

ANIO_MIN=1930
ANIO_MAX=2020

TIMEMIN = str(ANIO_MIN)+'-01-01'
TIMEMAX = str(ANIO_MAX)+'-12-31'

TIMEMIN_MEDIA = '1951-01-01'
TIMEMAX_MEDIA = '1980-12-31'

DATOS_SST='/home/dalia/Proyecto/BasesDatos/sst.mnmean.nc'
DATOS_PP_GPCC='/home/dalia/Proyecto/BasesDatos/precip.mon.total.v2018.nc'
DATOS_PP_CRU='/home/dalia/Proyecto/BasesDatos/cru_ts4.05.1901.2020.pre.dat.nc'

#%%

#Leo los datos de sea surface temperature

# Abro el archivo usando xarray
ds=xr.open_dataset(DATOS_SST)

#Recorte temporal
ds_recorte = ds.sel(time=slice(TIMEMIN,TIMEMAX))

#%%

#Calculo las anomalias de tsm mensuales para cada punto de reticula
#http://xarray.pydata.org/en/stable/generated/xarray.DataArray.groupby.html
tsm_anom=mm.groupby('time.month')-mm.groupby('time.month').mean('time')

#%%

#Remuestreo de los datos por estación
tsm_anom_est=tsm_anom.resample(time='QS-DEC').mean()['sst']

#%%

#Guardo en diccionarios
tsm_anom={ 'Verano (DJF)':tsm_anom_est.sel(time=tsm_anom_est["time.month"]==12),
           'Otoño (MAM)':tsm_anom_est.sel(time=tsm_anom_est["time.month"]==3),
           'Invierno (JJA)':tsm_anom_est.sel(time=tsm_anom_est["time.month"]==6),
           'Primavera (SON)':tsm_anom_est.sel(time=tsm_anom_est["time.month"]==9)}

#Exclude JF of first year and december of the last
tsm_anom['Verano (DJF)']=tsm_anom['Verano (DJF)'].sel(time=slice(tsm_anom['Verano (DJF)']['time'][1],
                                                                 tsm_anom['Verano (DJF)']['time'][-2]))
#De los datos de verano me quedo desde el segundo hasta el anteultimo ya que
#se agregaron meses vacios para poder completar el primer y ultimo cuarto.
# ds_djf=ds_djf.sel(time=slice(ds_djf["time"][1],ds_djf["time"][-2]))

#%%

#Ahora calculo la serie de anomalias de pp 

#%%

# Abro el archivo usando xarray
ds=xr.open_dataset(DATOS_PP_GPCC)

# Hago recorte
ds_recorte = ds.sel(lat=slice(LAT_NOR, LAT_SUR), lon=slice(LON_OESTE,LON_ESTE),
                        time=slice(TIMEMIN,TIMEMAX))

#%%

#Remuestreo los datos dividiendo por estaciones (DEF-MAM-JJA-SON)
ds_est = ds_recorte.resample(time="QS-DEC").sum()


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

ds_est= {"Verano (DJF)":ds_djf,
         "Otoño (MAM)":ds_mam,
         "Invierno (JJA)":ds_jja,
         "Primavera (SON)":ds_son}


#%%

#Anomalía para cada punto de reticula
ds_est_anom=[]

for key in ds_est.keys():
    ds_media = ds_est[key].sel(lat=slice(LAT_NOR_M,LAT_SUR_M),
                               lon=slice(LON_OESTE_M,LON_ESTE_M),
                               time=slice(TIMEMIN_MEDIA,TIMEMAX_MEDIA)).mean()
    ds_est_anom.append(ds_est[key]-ds_media)

#Paso de list a dict 
ds_est_anom= dict(zip(ds_est.keys(), ds_est_anom))

#%%

#Promedio regional con pesos segun latitud
ds_PromReg=[]
weights = np.cos(np.deg2rad(ds_recorte.lat))

for key in ds_est_anom.keys():
    ds_weighted = ds_est_anom[key].sel(lat=slice(LAT_NOR_M,LAT_SUR_M),
                                      lon=slice(LON_OESTE_M,LON_ESTE_M)).weighted(weights)
    ds_PromReg.append(ds_weighted.mean(("lon", "lat")))

                      
#Paso de list a dict
pp_anom_PromReg= dict(zip(ds_est.keys(), ds_PromReg))  

#%%

   