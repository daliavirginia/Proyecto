#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Oct  2 19:29:07 2021

@author: dalia
"""
#%%

#Load required libraries

import xarray as xr
import matplotlib.pyplot as plt
import numpy as np
from matplotlib import gridspec
import cartopy.crs as ccrs	
import cartopy.feature 	
from cartopy.mpl.ticker import LongitudeFormatter, LatitudeFormatter
import matplotlib
import matplotlib.path as mpath
from scipy import stats # LLamo a la fun stats de scipy 
#%%
def corr_xarray(ds3D,ds1D) :

    #Compute Pearson correlation for each gridpoint
    ds3D_stack=ds3D.stack(points=['lat', 'lon']) 
    corr = np.empty_like(ds3D_stack[0,:]) # Arrays vacios
    pval = np.empty_like(ds3D_stack[0,:])
    for k in range(ds3D_stack.shape[1]):
        x = ds1D
        y = ds3D_stack[:, k]
        [corr[k], pval[k]] = stats.pearsonr(x,y)       

    # Save trends of all ensemble members
    #Change to mm/day/summer/decade
    corr = np.reshape(corr, (len(ds3D['lat']), len(ds3D['lon'])))
    pval = np.reshape(pval, (len(ds3D['lat']), len(ds3D['lon'])))
    return corr,pval

#%%
# Valores constantes

# Para serie temporal de anomalias CENTRO ARGENTINA
LAT_SUR_M=-35
LAT_NOR_M=-26
LON_OESTE_M=294
LON_ESTE_M=299

#Para calculos 
LAT_SUR=-42
LAT_NOR=-16
LON_OESTE=287
LON_ESTE=306

#Para graficar 
LAT_SUR_G=-39
LAT_NOR_G=-19
LON_OESTE_G=290
LON_ESTE_G=303

#Periodo 
ANIO_MIN=1930
ANIO_MAX=2016

#Periodo en formato
TIMEMIN = str(ANIO_MIN)+'-01-01'
TIMEMAX = str(ANIO_MAX)+'-12-31'

#Periodo en el cual se calculará la media para las anomalías
TIMEMIN_MEDIA = '1951-01-01'
TIMEMAX_MEDIA = '1980-12-31'

#Path a los datos
DATOS_SST='/home/dalia/Proyecto/BasesDatos/sst.mnmean.nc'
DATOS_PP_GPCC='/home/dalia/Proyecto/BasesDatos/precip.mon.total.v2018.nc'
DATOS_PP_CRU='/home/dalia/Proyecto/BasesDatos/cru_ts4.05.1901.2020.pre.dat.nc'

KEYS=['Verano (DJF)', 'Otoño (MAM)', 'Invierno (JJA)', 'Primavera (SON)']
#%%

#Leo los datos de sea surface temperature

# Abro el archivo usando xarray
ds_tsm=xr.open_dataset(DATOS_SST)
print(len(ds_tsm.time))

#Recorte temporal
tsm_recorte = ds_tsm.sel(time=slice(TIMEMIN,TIMEMAX))['sst']
print(len(tsm_recorte.time))
#Agrupando por meses, le resto la media climatologica mensual a los datos. Asi
#obtengo anomalías de tsm. 
tsm_anom=tsm_recorte.groupby('time.month')-tsm_recorte.groupby('time.month').mean('time')
print(len(tsm_anom))
#%%

#Remuestreo de los datos por estación (DJF-MAM-JJA-SON)
tsm_anom_est=tsm_anom.resample(time='QS-DEC').mean()
print(len(tsm_anom_est))
print(tsm_anom_est)
#Guardo los datos separados por estación en un diccionario
tsm_anom_est={KEYS[0]:tsm_anom.sel(time=tsm_anom["time.month"]==12),
              KEYS[1]:tsm_anom.sel(time=tsm_anom["time.month"]==3),
              KEYS[2]:tsm_anom.sel(time=tsm_anom["time.month"]==6),
              KEYS[3]:tsm_anom.sel(time=tsm_anom["time.month"]==9)}
 
print(len(tsm_anom_est[KEYS[0]].time))
print(tsm_anom_est[KEYS[0]].time)
#Excluyo el primer y ultimo datos de verano (Tienen meses creados que no estan en los datos)
tsm_anom_est[KEYS[0]]=tsm_anom_est[KEYS[0]].sel(time=slice(tsm_anom_est[KEYS[0]]['time'][0],
                                                           tsm_anom_est[KEYS[0]]['time'][-2]))

print(len(tsm_anom_est[KEYS[0]].time))
#%%

##### ANOMALÍAS DE PP EN EL CENTRO DE ARGENTINA #########

#%%

# Abro el archivo usando xarray
ds_pp=xr.open_dataset(DATOS_PP_GPCC)
print(len(ds_pp.time))
#Recorte temporal y espacial
pp_recorte = ds_pp.sel(lat=slice(LAT_NOR, LAT_SUR), lon=slice(LON_OESTE,LON_ESTE),
                        time=slice(TIMEMIN,TIMEMAX))['precip']
print(len(pp_recorte))
#%%

#Remuestreo los datos dividiendo por estaciones (DEF-MAM-JJA-SON)
pp_est = pp_recorte.resample(time="QS-DEC").sum()
print(len(pp_est.time))
print(pp_est.time)

pp_est={KEYS[0]:pp_est.sel(time=pp_est["time.month"]==12),
        KEYS[1]:pp_est.sel(time=pp_est["time.month"]==3),
        KEYS[2]:pp_est.sel(time=pp_est["time.month"]==6),
        KEYS[3]:pp_est.sel(time=pp_est["time.month"]==9)}
print(len(pp_est[KEYS[0]].time))
#Excluyo el primer y ultimo datos de verano (Tienen meses creados que no estan en los datos)
pp_est[KEYS[0]]=pp_est[KEYS[0]].sel(time=slice(pp_est[KEYS[0]]['time'][1],
                                                pp_est[KEYS[0]]['time'][-2]))
print(len(pp_est[KEYS[0]].time))


#%%

#Anomalía para cada punto de reticula

#Inicializo lista vacia 
pp_est_anom=[]

#Ciclo para iterar en el diccionario de datos de pp
for key in KEYS:
    #Hago un recorte espacial para el centro de Arg y el periodo 51-80
    pp_media = pp_est[key].sel(lat=slice(LAT_NOR_M,LAT_SUR_M),
                               lon=slice(LON_OESTE_M,LON_ESTE_M),
                               time=slice(TIMEMIN_MEDIA,TIMEMAX_MEDIA)).mean()
    #A la lista vacía le agrego la anomalía para la estación correspondiente
    pp_est_anom.append(pp_est[key]-pp_media)
    
    del pp_media

#Paso de list a dict 
pp_est_anom= dict(zip(KEYS, pp_est_anom))

#%%

#Promedio regional con pesos segun latitud

#Inicializo lista vacía
pp_anom_est_PromReg=[]

#Data array con valores entre 0 y 1 para hacer promedio con peso x latitud
weights = np.cos(np.deg2rad(pp_recorte.lat))

for key in KEYS:
    ds_weighted = pp_est_anom[key].sel(lat=slice(LAT_NOR_M,LAT_SUR_M),
                                      lon=slice(LON_OESTE_M,LON_ESTE_M)).weighted(weights)
    pp_anom_est_PromReg.append(ds_weighted.mean(("lon", "lat")))

                      
#Paso de list a dict
pp_anom_est_PromReg= dict(zip(KEYS, pp_anom_est_PromReg))

#%%

### CORRELACIÓN ENTRE ANOMALÍAS DE TSM Y ANOMALÍAS DE PP EN EL CENTRO DE ARG ###

#%%

#pp_anom_est_PromReg --> 1 dimensión, tiempo
#tsm_anom_est --> 3 dimensiones, tiempo, lon y lat

#Con la función definida al principio

#Inicializo lista vacia
pp_tsm_corr=[]

for key in KEYS:
    
    tsm=tsm_anom_est[key]
    pp=pp_anom_est_PromReg[key]
    
    corr=corr_xarray(tsm,pp)
    
    pp_tsm_corr.append(corr)
    
    del corr,tsm,pp

#Hizo algo :o

#%%

#Con xarray.corr()

for key in KEYS:
    
    tsm=tsm_anom_est[key]
    pp=pp_anom_est_PromReg[key]
    
    corr=xr.corr(tsm,pp,'time') #No se si entiende que lo tiene que hacer
                                #para cada pto de grilla
    
    pp_tsm_corr.append(corr)
    
    del corr,tsm,pp

#

#%%

#Rellenando los nans con ceros 
