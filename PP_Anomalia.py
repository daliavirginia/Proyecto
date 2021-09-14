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

def trend_xarray(ds,VAR) :

    #Los parametros de la función son:ds instancia de Dataset y VAR la variable
    #que quiero del dataset
    
    #Apilo latitudes y longitudes
    ds_stack=ds.stack(points=['lat', 'lon']) 
    
    #Arrays vacios para todos los tiempos y el primer punto de grilla
    trends = np.empty_like(ds_stack[VAR][0,:]) # Arrays vacios
    pval = np.empty_like(ds_stack[VAR][0,:])
    b = np.empty_like(ds_stack[VAR][0,:])
    
    #Ciclo para completar los array vacios. Indice k recorre el rango de los
    #tiempos de ds_stack (shape[1] es para seleccionar la dim del tiempo)
    for k in range(ds_stack[VAR].shape[1]):
        y = ds_stack[VAR][:, k] #y es todos los puntos de grilla del tiempo k
        [trends[k], b[k], r_va, pval[k], z] = stats.linregress(np.arange(len(ds_stack['time'])), y)    

    # Le regreso la forma de latxlon
    tend = np.reshape(trends, (len(ds['lat']), len(ds['lon']))) 
    pv_tend = np.reshape(pval, (len(ds['lat']), len(ds['lon'])))
    interc = np.reshape(b, (len(ds['lat']), len(ds['lon'])))
    return tend,pv_tend,interc

#%%
# Constantes (en mayuscula)

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

ANIO_MIN=1930
ANIO_MAX=2016

TIMEMIN = str(ANIO_MIN)+'-01-01'
TIMEMAX = str(ANIO_MAX)+'-12-31'

TIMEMIN_MEDIA = '1951-01-01'
TIMEMAX_MEDIA = '1980-12-31'

# Directorio al archivo
DATOS_GPCC='/home/dalia/Proyecto/BasesDatos/precip.mon.total.v2018.nc'
DATOS_CMAP='/home/dalia/Proyecto/BasesDatos/CMAP_precip.mon.mean.nc'
DATOS_CRU='/home/dalia/Proyecto/BasesDatos/cru_ts4.05.1901.2020.pre.dat.nc'
# Directorio a salidas
SALIDAS='/home/dalia/Proyecto/Salidas/'

#%%

# Abro el archivo usando xarray
ds=xr.open_dataset(DATOS_GPCC)

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
ds_PromReg= dict(zip(ds_est.keys(), ds_PromReg))  

#Al calcular la media en cada punto de grilla debo considerar de que el area 
#en el cual se aplica el promedio es considerablemente mas grande en el Ecuador
#que en los Polos. Por eso, puedo aplicar una corrección, multiplicando la 
#media por el coseno de la latitud. Así se tendrá en cuenta que a menores latitudes
#el espacio es mayor y por ende el valor de pp pesa mas.

#%%

#Calculo la pendiente y la ordenada de la regresión lineal del promedio 
#regional para cada estación y se guardan las pendientes en "pend" y las 
#ordenadas al origen en "interc"

pend=[]
interc=[]

for key in ds_PromReg.keys():
    m,b,pv,rv,c1=stats.linregress(ds_PromReg[key]['time.year'],ds_PromReg[key]['precip'])
    pend.append(m)
    interc.append(b)

interc = dict(zip(ds_est.keys(),interc))
pend = dict(zip(ds_est.keys(),pend))

#%%

#Calculo los valores de anomalía de precipitación según la regresión lineal
# pp'=m*año+b para cada estación, y los guardo en "arreglo_lineal"

#Inicialización de lista vacía
arreglo_lineal = []

for key in interc.keys():
    a=pend[key]*ds_PromReg[key]['time.year']+interc[key]
    arreglo_lineal.append(a)

arreglo_lineal=dict(zip(ds_est.keys(),arreglo_lineal))

#%%

#Ahora a los valores originales de anomalía se le resta la tendencia lineal y 
#se almacena en "anom_st"

#Inicalización de lista vacia
anom_st=[]

#Ciclo para iterar diccionarios y realizar los calculos
for key in ds_PromReg.keys():
    anom_st.append(ds_PromReg[key]['precip']-arreglo_lineal[key])
    
anom_st=dict(zip(ds_est.keys(),anom_st))
    
#%%

#Grafico las anomalías

fig, ax = plt.subplots(figsize=(2*5.5,2*4.8))

gs = gridspec.GridSpec(2,2)

i=0
for key in ds_PromReg.keys():
    
    ax = plt.subplot(gs[i])
    bars = ax.bar(ds_PromReg[key]["precip"]["time.year"],
                  ds_PromReg[key]["precip"])
    
    for j in range(len(bars)):
    
        if bars[j].get_height()>0:
            bars[j].set_fc("darkturquoise")
        else:
            bars[j].set_fc("tomato")
            
    plt.ylim(-160,160)
    ax.set_ylabel("mm")
    ax.set_xlabel("Año")
    
    ax.set_title(str(key))
    
    i=i+1

fig.suptitle("Anomalía de precipitación", fontsize=18)

fig.savefig(SALIDAS+"AnomaliaPP_CentroArgentina_GPCC.png", dpi=300,  bbox_inches='tight')
plt.show()


#%%
#Grafico de las anomalías sin tendencia lineal

#Se crea una instancia de pyplot
fig, ax = plt.subplots(figsize=(2*5.5,2*4.8))

#Grilla de 2x2 para graficar cada estación
gs = gridspec.GridSpec(2,2)

#Ciclo para graficar cada estación
i=0
for key in ds_PromReg.keys():
    
    ax = plt.subplot(gs[i])
    bars = ax.bar(anom_st[key]["time.year"],
                  anom_st[key])
    
    #Ciclo para que las barras negativas sean rojas, y las positivas azules
    for j in range(len(bars)):
    
        if bars[j].get_height()>0:
            bars[j].set_fc("darkturquoise")
        else:
            bars[j].set_fc("tomato")
            
    plt.ylim(-160,160)
    ax.set_ylabel("mm")
    ax.set_xlabel("Año")
    
    ax.set_title(str(key))
    
    i=i+1

fig.suptitle("Anomalía de precipitación sin tendencia", fontsize=18)

fig.savefig(SALIDAS+"AnomaliaPP_SinTendencia_CentroArgentina_GPCC.png", dpi=300,  bbox_inches='tight')
plt.show()

#%%

#Ahora se le restará la tendencia a cada punto de grilla 

#ds_est_anom almacena las anomalías para cada punto de grilla

#Cálculo de la pendiente y ordenada de regresión lineal, sobreescribiendo los
#diccionarios "pend" e "interc" definidos antes

for key in ds_est_anom.keys():
    m,b,pv=trend_xarray(ds_est_anom[key],'precip')
    pend[key]=m
    interc[key]=b

#%%


for key in ds_est_anom.keys():    
    a=[]
    
    for i in range(len(ds_est_anom[key]['time.year'])-1):
        year=int(ds_est_anom[key]['time.year'][i])
        a.append(pend[key]*year+interc[key])
        
    arreglo_lineal[key]=a

#%%
anom_st_2D=[]

for key in ds_est_anom.keys():
    anom_st_2D.append(ds_est_anom[key]['precip']-arreglo_lineal[key][0])
    
anom_st_2D=dict(zip(ds_est.keys(),anom_st_2D))

#%%

#CORRELACION para cada punto de reticula con el promedio regional


#%%
corr=[]
pval=[]

for key in ds_PromReg.keys():
    cor,pv=corr_xarray(anom_st_2D[key], anom_st[key])
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
    
    im=ax.contourf(lons, lats, corr[key],cmap="RdYlGn", levels=clevs,extend='both',transform=crs_latlon)
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
    
    COLOR='black'

    ax.plot([LON_OESTE_M, LON_OESTE_M], [LAT_SUR_M,
    LAT_NOR_M],color=COLOR, linestyle='--',alpha=0.9,
    transform=ccrs.PlateCarree())
    ax.plot([LON_ESTE_M, LON_ESTE_M], [LAT_SUR_M,
    LAT_NOR_M],color=COLOR,
    linestyle='--',alpha=0.9,transform=ccrs.PlateCarree())
    ax.plot([LON_OESTE_M, LON_ESTE_M], [LAT_SUR_M,
    LAT_SUR_M],color=COLOR,
    linestyle='--',alpha=0.9,transform=ccrs.PlateCarree())
    ax.plot([LON_OESTE_M, LON_ESTE_M], [LAT_NOR_M,
    LAT_NOR_M],color=COLOR,
    linestyle='--',alpha=0.9,transform=ccrs.PlateCarree())
    
    i=i+1

cbar_ax = fig.add_axes([0.98, 0.1, 0.03, 0.8])

fig.colorbar(im, cax=cbar_ax,orientation='vertical') 
fig.suptitle("Correlación con promedio regional \n(sin tendencia)",fontsize=14)

fig.savefig(SALIDAS+"Correlacion_PP_GPCC.png", dpi=300, bbox_inches='tight')

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


# Hago recorte
ds_recorte = ds.sel(lat=slice(LAT_SUR, LAT_NOR), lon=slice(LON_OESTE,LON_ESTE),
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
    ds_media = ds_est[key].sel(lat=slice(LAT_SUR_M,LAT_NOR_M),
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
    ds_weighted = ds_est_anom[key].sel(lat=slice(LAT_SUR_M,LAT_NOR_M),
                                      lon=slice(LON_OESTE_M,LON_ESTE_M)).weighted(weights)
    ds_PromReg.append(ds_weighted.mean(("lon", "lat")))

                      
#Paso de list a dict
ds_PromReg= dict(zip(ds_est.keys(), ds_PromReg))  

#Al calcular la media en cada punto de grilla debo considerar de que el area 
#en el cual se aplica el promedio es considerablemente mas grande en el Ecuador
#que en los Polos. Por eso, puedo aplicar una corrección, multiplicando la 
#media por el coseno de la latitud. Así se tendrá en cuenta que a menores latitudes
#el espacio es mayor y por ende el valor de pp pesa mas.


#%%

#Grafico las anomalías

fig, ax = plt.subplots(figsize=(2*5.5,2*4.8))

gs = gridspec.GridSpec(2,2)

i=0
for key in ds_PromReg.keys():
    
    ax = plt.subplot(gs[i])
    bars = ax.bar(ds_PromReg[key]["pre"]["time.year"],
                  ds_PromReg[key]["pre"])
    
    for j in range(len(bars)):
    
        if bars[j].get_height()>0:
            bars[j].set_fc("darkturquoise")
        else:
            bars[j].set_fc("tomato")
            
    plt.ylim(-160,160)
    ax.set_ylabel("mm")
    ax.set_xlabel("Año")
    
    ax.set_title(str(key))
    
    i=i+1

fig.suptitle("Anomalía de precipitación", fontsize=18)

fig.savefig(SALIDAS+"AnomaliaPP_CentroArgentina_CRU.png", dpi=300,  bbox_inches='tight')
plt.show()


#%%
#Grafico de las anomalías sin tendencia lineal

#Se crea una instancia de pyplot
fig, ax = plt.subplots(figsize=(2*5.5,2*4.8))

#Grilla de 2x2 para graficar cada estación
gs = gridspec.GridSpec(2,2)

#Ciclo para graficar cada estación
i=0
for key in ds_PromReg.keys():
    
    ax = plt.subplot(gs[i])
    bars = ax.bar(anom_st[key]["time.year"],
                  anom_st[key])
    
    #Ciclo para que las barras negativas sean rojas, y las positivas azules
    for j in range(len(bars)):
    
        if bars[j].get_height()>0:
            bars[j].set_fc("darkturquoise")
        else:
            bars[j].set_fc("tomato")
            
    plt.ylim(-160,160)
    ax.set_ylabel("mm")
    ax.set_xlabel("Año")
    
    ax.set_title(str(key))
    
    i=i+1

fig.suptitle("Anomalía de precipitación sin tendencia", fontsize=18)

fig.savefig(SALIDAS+"AnomaliaPP_SinTendencia_CentroArgentina_CRU.png", dpi=300,  bbox_inches='tight')
plt.show()

#%%

#Ahora se le restará la tendencia a cada punto de grilla 

#ds_est_anom almacena las anomalías para cada punto de grilla

#Cálculo de la pendiente y ordenada de regresión lineal, sobreescribiendo los
#diccionarios "pend" e "interc" definidos antes

for key in ds_est_anom.keys():
    m,b,pv=trend_xarray(ds_est_anom[key],'pre')
    pend[key]=m
    interc[key]=b

#%%


for key in ds_est_anom.keys():    
    a=[]
    
    for i in range(len(ds_est_anom[key]['time.year'])-1):
        year=int(ds_est_anom[key]['time.year'][i])
        a.append(pend[key]*year+interc[key])
        
    arreglo_lineal[key]=a

#%%
anom_st_2D=[]

for key in ds_est_anom.keys():
    anom_st_2D.append(ds_est_anom[key]['pre']-arreglo_lineal[key][0])
    
anom_st_2D=dict(zip(ds_est.keys(),anom_st_2D))

#%%

#CORRELACION para cada punto de reticula con el promedio regional


#%%
corr=[]
pval=[]

for key in ds_PromReg.keys():
    cor,pv=corr_xarray(anom_st_2D[key], anom_st[key])
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
    
    im=ax.contourf(lons, lats, corr[key],cmap="RdYlGn", levels=clevs,extend='both',transform=crs_latlon)
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
    
    COLOR='black'

    ax.plot([LON_OESTE_M, LON_OESTE_M], [LAT_SUR_M,
    LAT_NOR_M],color=COLOR, linestyle='--',alpha=0.9,
    transform=ccrs.PlateCarree())
    ax.plot([LON_ESTE_M, LON_ESTE_M], [LAT_SUR_M,
    LAT_NOR_M],color=COLOR,
    linestyle='--',alpha=0.9,transform=ccrs.PlateCarree())
    ax.plot([LON_OESTE_M, LON_ESTE_M], [LAT_SUR_M,
    LAT_SUR_M],color=COLOR,
    linestyle='--',alpha=0.9,transform=ccrs.PlateCarree())
    ax.plot([LON_OESTE_M, LON_ESTE_M], [LAT_NOR_M,
    LAT_NOR_M],color=COLOR,
    linestyle='--',alpha=0.9,transform=ccrs.PlateCarree())
    
    i=i+1

cbar_ax = fig.add_axes([0.98, 0.1, 0.03, 0.8])

fig.colorbar(im, cax=cbar_ax,orientation='vertical') 
fig.suptitle("Correlación con promedio regional \n(sin tendencia)",fontsize=14)

fig.savefig(SALIDAS+"CorrelacionPP_SinTendencia_CRU.png", dpi=300, bbox_inches='tight')