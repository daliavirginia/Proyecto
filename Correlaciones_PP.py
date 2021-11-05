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

def trend_xarray(darray) :

    #Los parametros de la función son:darray instancia de Dataset y VAR la variable
    #que quiero del dataset
    
    #Apilo latitudes y longitudes
    darray_stack=darray.stack(points=['lat', 'lon']) 
    
    #Arrays vacios para todos los tiempos y el primer punto de grilla
    trendarray = np.empty_like(darray_stack[0,:]) # Arrays vacios
    pval = np.empty_like(darray_stack[0,:])
    b = np.empty_like(darray_stack[0,:])
    
    #Ciclo para completar los array vacios. Indice k recorre el rango de los
    #tiempos de darray_stack (shape[1] es para seleccionar la dim del tiempo)
    for k in range(darray_stack.shape[1]):
        y = darray_stack[:, k] #y es todos los puntos de grilla del tiempo k
        [trendarray[k], b[k], r_va, pval[k], z] = stats.linregress(np.arange(len(darray_stack['time'])), y)    

    # Le regreso la forma de latxlon
    tend = np.reshape(trendarray, (len(darray['lat']), len(darray['lon']))) 
    pv_tend = np.reshape(pval, (len(darray['lat']), len(darray['lon'])))
    interc = np.reshape(b, (len(darray['lat']), len(darray['lon'])))
    return tend,pv_tend,interc

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

#Periodo largo
ANIO_MIN1=1979
ANIO_MAX1=2016

#Periodo corto
ANIO_MIN2=1981
ANIO_MAX2=2010

#Periodo en formato
TIMEMIN1 = str(ANIO_MIN1)+'-01-01'
TIMEMAX1 = str(ANIO_MAX1)+'-12-31'

#Periodo en formato
TIMEMIN2 = str(ANIO_MIN2)+'-01-01'
TIMEMAX2 = str(ANIO_MAX2)+'-12-31'

#Periodo en el cual se calculará la media para las anomalías
TIMEMIN_MEDIA = '1951-01-01'
TIMEMAX_MEDIA = '1980-12-31'

#Path a los datos
DATOS_SST='/home/dalia/Proyecto/BasesDatos/sst.mnmean.nc'
DATOS_PP_GPCC='/home/dalia/Proyecto/BasesDatos/precip.mon.total.v2018.nc'
DATOS_PP_CRU='/home/dalia/Proyecto/BasesDatos/cru_ts4.05.1901.2020.pre.dat.nc'
DATOS_GEOP='/home/dalia/Proyecto/BasesDatos/ERA5/geopotential_monthly_averaged.nc'
DATOS_MSLP='/home/dalia/Proyecto/BasesDatos/ERA5/mean_sea_level_pressure_monthly_averaged.nc'
#Path para guardar figuras
SALIDAS='/home/dalia/Proyecto/Salidas/'
#Vector con nombres de las estaciones 
KEYS=['Verano (DJF)', 'Otoño (MAM)', 'Invierno (JJA)', 'Primavera (SON)']
#%%

#Selección de periodo

per=1

if (per==1):
    
    TIMEMIN=TIMEMIN1
    TIMEMAX=TIMEMAX1
    PER='(1930-2016)'
    
else:
    
    TIMEMIN=TIMEMIN2
    TIMEMAX=TIMEMAX2
    PER='(1981-2010)'

#%%

############## ANOMALÍAS DE SEA SURFACE TEMPERATURE ######################

# Abro el archivo usando xarray
ds_tsm=xr.open_dataset(DATOS_SST)
print(len(ds_tsm.time))
print(dir(ds_tsm))
print(ds_tsm.info)
#Recorte temporal
tsm_recorte = ds_tsm.sel(time=slice(TIMEMIN,TIMEMAX))['sst']
print(len(tsm_recorte.time))
#Agrupando por meses, le resto la media climatologica mensual a los datos. Asi
#obtengo anomalías de tsm. 
tsm_anom=tsm_recorte.groupby('time.month')-tsm_recorte.groupby('time.month').mean('time')


#Remuestreo de los datos por estación (DJF-MAM-JJA-SON)
tsm_anom_est=tsm_anom.resample(time='QS-DEC').mean()

#Guardo los datos separados por estación en un diccionario
tsm_anom_est={KEYS[0]:tsm_anom.sel(time=tsm_anom["time.month"]==12),
              KEYS[1]:tsm_anom.sel(time=tsm_anom["time.month"]==3),
              KEYS[2]:tsm_anom.sel(time=tsm_anom["time.month"]==6),
              KEYS[3]:tsm_anom.sel(time=tsm_anom["time.month"]==9)}
 

#Excluyo el primer y ultimo datos de verano (Tienen meses creados que no estan en los datos)
tsm_anom_est[KEYS[0]]=tsm_anom_est[KEYS[0]].sel(time=slice(tsm_anom_est[KEYS[0]]['time'][0],
                                                           tsm_anom_est[KEYS[0]]['time'][-2]))

#%%

############ ANOMALÍAS DE PP EN EL CENTRO DE ARGENTINA #####################

# Abro el archivo usando xarray
ds_pp=xr.open_dataset(DATOS_PP_GPCC)
#Recorte temporal y espacial
pp_recorte = ds_pp.sel(lat=slice(LAT_NOR, LAT_SUR), lon=slice(LON_OESTE,LON_ESTE))['precip']

#Remuestreo los datos dividiendo por estaciones (DEF-MAM-JJA-SON)
pp_est = pp_recorte.resample(time="QS-DEC").sum()

pp_est={KEYS[0]:pp_est.sel(time=pp_est["time.month"]==12),
        KEYS[1]:pp_est.sel(time=pp_est["time.month"]==3),
        KEYS[2]:pp_est.sel(time=pp_est["time.month"]==6),
        KEYS[3]:pp_est.sel(time=pp_est["time.month"]==9)}

#Excluyo el primer y ultimo datos de verano (Tienen meses creados que no estan en los datos)
pp_est[KEYS[0]]=pp_est[KEYS[0]].sel(time=slice(pp_est[KEYS[0]]['time'][1],
                                                pp_est[KEYS[0]]['time'][-2]))


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
    pp_est_anom.append(pp_est[key].sel(time=slice(TIMEMIN,TIMEMAX))-pp_media)
    
    del pp_media

#Paso de list a dict 
pp_est_anom= dict(zip(KEYS, pp_est_anom))


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

############################# ANOMALÍAS DE GEOPOTENCIAL ######################

# Abro el archivo usando xarray
ds_geo=xr.open_dataset(DATOS_GEOP)
print(len(ds_tsm.time))
print(ds_geo.var)

alt=850
# #Recorte temporal
geo_recorte = ds_geo.sel(time=slice(TIMEMIN,TIMEMAX))['z'].sel(level=alt)
print(len(geo_recorte.time))
# #Agrupando por meses, le resto la media climatologica mensual a los datos. Asi
# #obtengo anomalías de tsm. 
geo_anom=geo_recorte.groupby('time.month')-geo_recorte.groupby('time.month').mean('time')

#Remuestreo de los datos por estación (DJF-MAM-JJA-SON)
geo_anom_est=geo_anom.resample(time='QS-DEC').mean()

#Guardo los datos separados por estación en un diccionario
geo_anom_est={KEYS[0]:geo_anom.sel(time=geo_anom["time.month"]==12),
              KEYS[1]:geo_anom.sel(time=geo_anom["time.month"]==3),
              KEYS[2]:geo_anom.sel(time=geo_anom["time.month"]==6),
              KEYS[3]:geo_anom.sel(time=geo_anom["time.month"]==9)}
 

#Excluyo el primer y ultimo datos de verano (Tienen meses creados que no estan en los datos)
geo_anom_est[KEYS[0]]=geo_anom_est[KEYS[0]].sel(time=slice(geo_anom_est[KEYS[0]]['time'][0],
                                                           geo_anom_est[KEYS[0]]['time'][-1]))


#%%

########################## CORRELACIÓN  ######################################

#pp_anom_est_PromReg --> 1 dimensión, tiempo
#tsm_anom_est --> 3 dimensiones, tiempo, lon y lat

#Con la función definida al principio

#Inicializo lista vacia
pp_tsm_corr=[]

#Con xarray.corr()

for key in KEYS:
    
    tsm=tsm_anom_est[key]
    pp=pp_anom_est_PromReg[key]
    
    corr=xr.corr(tsm,pp,'time') #No se si entiende que lo tiene que hacer
                                #para cada pto de grilla
    
    pp_tsm_corr.append(corr)
    
    del corr,tsm,pp


pp_tsm_corr=dict(zip(KEYS,pp_tsm_corr))
# #%%

# #Tamaño para mapas del centro de Argentina
# fig, ax = plt.subplots()

# #Define grid for subplots
# gs = gridspec.GridSpec(2,2)     

# #latitudes and longitudes to plot
# lons, lats = np.meshgrid(pp_tsm_corr[key]['lon'], pp_tsm_corr[key]['lat'])

# #Definimos los niveles para los contornos (inferior, superior, longitud de los intervalos)
# #Si utilizamos las funciones ds.min ds.max podemos darnos una idea de cómo definirlos
# lev_inf = -0.6
# lev_sup = 0.6
# lev_int =0.1
# clevs = np.arange(lev_inf, lev_sup+lev_int, lev_int)

# i=0
# for key in KEYS:
    
#     ax = plt.subplot(gs[i],projection=ccrs.PlateCarree(central_longitude=180))
#     crs_latlon = ccrs.PlateCarree()
#     #ax.set_extent([LON_OESTE_G, LON_ESTE_G, LAT_SUR_G, LAT_NOR_G], crs=crs_latlon)
#     #ax.add_feature(cartopy.feature.OCEAN, zorder=100, edgecolor='k')
#     ax.add_feature(cartopy.feature.COASTLINE, alpha=.5)
#     #ax.add_feature(cartopy.feature.BORDERS, linestyle='-', alpha=.1)
#     # Create a feature for States/Admin 1 regions at 1:50m from Natural Earth
#    # states_provinces = cartopy.feature.NaturalEarthFeature(
#       #  category='cultural',
#       #  name='admin_1_states_provinces_lines',
#       #  scale='10m',
#       #  facecolor='none')
    
#     #ax.add_feature(states_provinces, edgecolor='darkslategrey', linewidths=0.3)
    
#     im=ax.contourf(lons, lats, pp_tsm_corr[key],cmap="RdYlGn", levels=clevs,extend='both',transform=crs_latlon)
#     #ax.contour(lons, lats, pval[key],levels=[0.1],colors='k',linewidths=0.5 , transform=crs_latlon) 
    
#     #ax.set_xticks(np.arange(0,358,100), crs=crs_latlon)
#     #ax.set_yticks(np.arange(-88,88,100), crs=crs_latlon)
#     ax.grid(which='both', linewidth=0.3, linestyle='-')
#     ax.tick_params(axis='both', which='major', labelsize=8)
#     lon_formatter = LongitudeFormatter(zero_direction_label=True)
#     lat_formatter = LatitudeFormatter()
#     ax.xaxis.set_major_formatter(lon_formatter)
#     ax.yaxis.set_major_formatter(lat_formatter)    
    
#     ax.set_title(str(key),fontsize=10)
    
#     COLOR='black'

#     i=i+1

# cbar_ax = fig.add_axes([0.98, 0.1, 0.03, 0.8])

# fig.colorbar(im, cax=cbar_ax,orientation='vertical') 
# fig.suptitle("Correlación TSM y PP (Centro Argentina) \n"+PER ,fontsize=11)

# fig.savefig(SALIDAS+"Correlacion_TSM_PP_CRU_"+PER+".png", dpi=300, bbox_inches='tight')

#%%

#################### ELIMINACIÓN DE LA TENDENCIA LINEAL ######################

#ANOMALÍAS DE PRECIPITACIÓN

#(1) Calculo la pendiente y la ordenada de la regresión lineal del promedio 
#    regional para cada estación y se guardan las pendientes en "pend" y las 
#    ordenadas al origen en "interc"

pend=[]
interc=[]

for key in KEYS:
    m,b,pv,rv,c1=stats.linregress(pp_anom_est_PromReg[key]['time.year'],pp_anom_est_PromReg[key])
    pend.append(m)
    interc.append(b)
    
#Remuevo las variables temporales usadas en el ciclo
del m,b,pv,rv,c1

interc = dict(zip(KEYS,interc))
pend = dict(zip(KEYS,pend)) 

#(2) Calculo los valores de anomalía de precipitación según la regresión lineal
#    pp'=m*año+b para cada estación, y los guardo en "arreglo_lineal"

#Inicialización de lista vacía
arreglo_lineal = []

for key in KEYS:
    a=pend[key]*pp_anom_est_PromReg[key]['time.year']+interc[key]
    arreglo_lineal.append(a)

del a

arreglo_lineal=dict(zip(KEYS,arreglo_lineal))

#(3) Ahora a los valores originales de anomalía se le resta la tendencia lineal y 
#    se almacena en "anom_st"

#Inicalización de lista vacia
pp_anom_est_ST_PromReg=[]

#Ciclo para iterar diccionarios y realizar los calculos
for key in KEYS:
    pp_anom_est_ST_PromReg.append(pp_anom_est_PromReg[key]-arreglo_lineal[key])
    
pp_anom_est_ST_PromReg=dict(zip(KEYS,pp_anom_est_ST_PromReg))

#ANOMALÍAS DE TSM

#(1) Calculo con la función 

for key in KEYS:
    m,pv,b=trend_xarray(tsm_anom_est[key])
    pend[key]=m
    interc[key]=b

for key in KEYS:    
    a=[]
    for i in range(len(tsm_anom_est[key]['time.year'])-1):
        year=int(tsm_anom_est[key]['time.year'][i])
        a.append(pend[key]*year+interc[key])
    arreglo_lineal[key]=a

tsm_anom_est_ST=[]
for key in KEYS:
    tsm_anom_est_ST.append(tsm_anom_est[key]-arreglo_lineal[key][0])
    
tsm_anom_est_ST=dict(zip(KEYS,tsm_anom_est_ST))

#ANOMALÍAS DE GEOPOTENCIAL

for key in KEYS:
    m,pv,b=trend_xarray(tsm_anom_est[key])
    pend[key]=m
    interc[key]=b

for key in KEYS:    
    a=[]
    for i in range(len(tsm_anom_est[key]['time.year'])-1):
        year=int(tsm_anom_est[key]['time.year'][i])
        a.append(pend[key]*year+interc[key])
    arreglo_lineal[key]=a

tsm_anom_est_ST=[]
for key in KEYS:
    tsm_anom_est_ST.append(tsm_anom_est[key]-arreglo_lineal[key][0])
    
tsm_anom_est_ST=dict(zip(KEYS,tsm_anom_est_ST))

#ANOMALÍAS GEOPOTENCIAL

for key in KEYS:
    m,pv,b=trend_xarray(geo_anom_est[key])
    pend[key]=m
    interc[key]=b

for key in KEYS:    
    a=[]
    for i in range(len(geo_anom_est[key]['time.year'])-1):
        year=int(geo_anom_est[key]['time.year'][i])
        a.append(pend[key]*year+interc[key])
    arreglo_lineal[key]=a

geo_anom_est_ST=[]
for key in KEYS:
    geo_anom_est_ST.append(geo_anom_est[key]-arreglo_lineal[key][0])
    
geo_anom_est_ST=dict(zip(KEYS,geo_anom_est_ST))

#%% 

´###################### CORRELACIÓN SIN TENDENCIA ############################

# PRECIPITACIÓN + TSM

#Inicializo lista vacia
pp_tsm_corr_ST=[]

#Con xarray.corr()

for key in KEYS:
    
    tsm=tsm_anom_est_ST[key]
    pp=pp_anom_est_ST_PromReg[key]
    
    corr=xr.corr(tsm,pp,'time') #No se si entiende que lo tiene que hacer
                                #para cada pto de grilla
    
    pp_tsm_corr_ST.append(corr)
    
    del corr,tsm,pp

pp_tsm_corr_ST=dict(zip(KEYS,pp_tsm_corr_ST))

# PRECIPITACIÓN + GEOPOTENCIAL

#Inicializo lista vacia
pp_geo_corr_ST=[]

#Con xarray.corr()

for key in KEYS:
    
    geo=geo_anom_est_ST[key]
    pp=pp_anom_est_ST_PromReg[key]
    
    corr=xr.corr(geo,pp,'time') #No se si entiende que lo tiene que hacer
                                #para cada pto de grilla
    
    pp_geo_corr_ST.append(corr)
    
    del corr,geo,pp

pp_geo_corr_ST=dict(zip(KEYS,pp_geo_corr_ST))

#%%

var = pp_geo_corr_ST
VAR = 'GEOP(ERA5)'

#var = pp_tsm_corr_ST
#VAR = 'TSM()'




#Tamaño para mapas del centro de Argentina
fig, ax = plt.subplots()

#Define grid for subplots
gs = gridspec.GridSpec(2,2)     

#latitudes and longitudes to plot
lons, lats = np.meshgrid(var[key]['lon'], var[key]['lat'])

#Definimos los niveles para los contornos (inferior, superior, longitud de los intervalos)
#Si utilizamos las funciones ds.min ds.max podemos darnos una idea de cómo definirlos
lev_inf = -0.6
lev_sup = 0.6
lev_int =0.1
clevs = np.arange(lev_inf, lev_sup+lev_int, lev_int)

i=0
for key in KEYS:
    
    ax = plt.subplot(gs[i],projection=ccrs.PlateCarree(central_longitude=180))
    crs_latlon = ccrs.PlateCarree()
    #ax.set_extent([LON_OESTE_G, LON_ESTE_G, LAT_SUR_G, LAT_NOR_G], crs=crs_latlon)
    #ax.add_feature(cartopy.feature.OCEAN, zorder=100, edgecolor='k')
    ax.add_feature(cartopy.feature.COASTLINE, alpha=.5)
    #ax.add_feature(cartopy.feature.BORDERS, linestyle='-', alpha=.05)
    # Create a feature for States/Admin 1 regions at 1:50m from Natural Earth
    # states_provinces = cartopy.feature.NaturalEarthFeature(
      #  category='cultural',
      #  name='admin_1_states_provinces_lines',
      #  scale='10m',
      #  facecolor='none')
    
    #ax.add_feature(states_provinces, edgecolor='darkslategrey', linewidths=0.3)
    
    im=ax.contourf(lons, lats, var[key],cmap="RdYlGn", levels=clevs,extend='both',transform=crs_latlon)
    #ax.contour(lons, lats, pval[key],levels=[0.1],colors='k',linewidths=0.5 , transform=crs_latlon) 
    
    ax.set_xticks(np.arange(0,358,100), crs=crs_latlon)
    ax.set_yticks(np.arange(-88,88,100), crs=crs_latlon)
    ax.grid(which='both', linewidth=0.3, linestyle='-')
    ax.tick_params(axis='both', which='major', labelsize=8)
    lon_formatter = LongitudeFormatter(zero_direction_label=True)
    lat_formatter = LatitudeFormatter()
    ax.xaxis.set_major_formatter(lon_formatter)
    ax.yaxis.set_major_formatter(lat_formatter)    
    
    ax.set_title(str(key), fontsize=10)
    
    COLOR='black'

    i=i+1

cbar_ax = fig.add_axes([0.98, 0.1, 0.03, 0.8])

fig.colorbar(im, cax=cbar_ax,orientation='vertical') 
fig.suptitle("Correlación TSM y PP (Centro Argentina) sin tendencia \n"+PER,fontsize=11)

fig.savefig(SALIDAS+'Correlacion_'+VAR+'_PP(CRU)_SinTendencia_'+PER+".png", dpi=300, bbox_inches='tight')


#%%

######################### PROMEDIO MÓVIL ####################################

# Hago promedio movil para filtrar la variabilidad de mayor frecuencia y
# y calculo la correlación

tsm_anom_est_ST_PM7=[]
pp_anom_est_ST_PromReg_PM7=[]

for key in KEYS:
    pm_tsm = tsm_anom_est_ST[key].rolling(time=7).mean()
    pm_pp = pp_anom_est_ST_PromReg[key].rolling(time=7).mean()
    
    tsm_anom_est_ST_PM7.append(pm_tsm)
    pp_anom_est_ST_PromReg_PM7.append(pm_pp)
    
    del pm_tsm, pm_pp
    
tsm_anom_est_ST_PM7 = dict(zip(KEYS,tsm_anom_est_ST_PM7))
pp_anom_est_ST_PromReg_PM7 = dict(zip(KEYS,pp_anom_est_ST_PromReg_PM7))

############## CORRELACIÓN CON PROMEDIO MOVIL Y S/TENDENCIA ##################
#Correlacion con promedio movil 7 años

#Con la función definida al principio

#Inicializo lista vacia
pp_tsm_corr_ST_PM7=[]

#Con xarray.corr()

for key in KEYS:
    
    tsm=tsm_anom_est_ST_PM7[key]
    pp=pp_anom_est_ST_PromReg_PM7[key]
    
    corr=xr.corr(tsm,pp,'time') #No se si entiende que lo tiene que hacer
                                #para cada pto de grilla
    
    pp_tsm_corr_ST_PM7.append(corr)
    
    del corr,tsm,pp


pp_tsm_corr_ST_PM7=dict(zip(KEYS,pp_tsm_corr_ST_PM7))



#%%

#Tamaño para mapas del centro de Argentina
fig, ax = plt.subplots()

#Define grid for subplots
gs = gridspec.GridSpec(2,2)     

#latitudes and longitudes to plot
lons, lats = np.meshgrid(pp_tsm_corr_ST_PM7[key]['lon'], pp_tsm_corr_ST_PM7[key]['lat'])

#Definimos los niveles para los contornos (inferior, superior, longitud de los intervalos)
#Si utilizamos las funciones ds.min ds.max podemos darnos una idea de cómo definirlos
lev_inf = -0.6
lev_sup = 0.6
lev_int =0.1
clevs = np.arange(lev_inf, lev_sup+lev_int, lev_int)

i=0
for key in KEYS:
    
    ax = plt.subplot(gs[i],projection=ccrs.PlateCarree(central_longitude=180))
    crs_latlon = ccrs.PlateCarree()
    #ax.set_extent([LON_OESTE_G, LON_ESTE_G, LAT_SUR_G, LAT_NOR_G], crs=crs_latlon)
    #ax.add_feature(cartopy.feature.OCEAN, zorder=100, edgecolor='k')
    ax.add_feature(cartopy.feature.COASTLINE, alpha=.5)
    #ax.add_feature(cartopy.feature.BORDERS, linestyle='-', alpha=.05)
    # Create a feature for States/Admin 1 regions at 1:50m from Natural Earth
   # states_provinces = cartopy.feature.NaturalEarthFeature(
      #  category='cultural',
      #  name='admin_1_states_provinces_lines',
      #  scale='10m',
      #  facecolor='none')
    
    #ax.add_feature(states_provinces, edgecolor='darkslategrey', linewidths=0.3)
    
    im=ax.contourf(lons, lats, pp_tsm_corr_ST_PM7[key],cmap="RdYlGn", levels=clevs,extend='both',transform=crs_latlon)
    #ax.contour(lons, lats, pval[key],levels=[0.1],colors='k',linewidths=0.5 , transform=crs_latlon) 
    
    ax.set_xticks(np.arange(0,358,100), crs=crs_latlon)
    ax.set_yticks(np.arange(-88,88,100), crs=crs_latlon)
    ax.grid(which='both', linewidth=0.3, linestyle='-')
    ax.tick_params(axis='both', which='major', labelsize=8)
    lon_formatter = LongitudeFormatter(zero_direction_label=True)
    lat_formatter = LatitudeFormatter()
    ax.xaxis.set_major_formatter(lon_formatter)
    ax.yaxis.set_major_formatter(lat_formatter)    
    
    ax.set_title(str(key), fontsize=10)
    
    COLOR='black'

    i=i+1

cbar_ax = fig.add_axes([0.98, 0.1, 0.03, 0.8])

fig.colorbar(im, cax=cbar_ax,orientation='vertical') 
fig.suptitle("Correlación TSM y PP (Centro Argentina) con prom. movil 7 años (s/ tendencia)\n"+PER,fontsize=11)

fig.savefig(SALIDAS+"Correlacion_TSM_PP_CRU_PromedioMovil7años"+PER+".png", dpi=300, bbox_inches='tight')

#%%
