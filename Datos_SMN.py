#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Sep 17 22:25:13 2021

@author: dalia
"""
#%%
#Load required libraries

import pandas as pd
import xarray as xr
import matplotlib.pyplot as plt
from matplotlib import gridspec
from scipy import stats
import numpy as np
import matplotlib.dates as mdates
import datetime
from scipy.fftpack import fft
import matplotlib.patches as mpatches

import cartopy.crs as ccrs
import cartopy.feature
from cartopy.mpl.ticker import LongitudeFormatter, LatitudeFormatter
#%%

#Función para contar faltantes 
def missing_per(df):
    
    f = df.isnull().sum() #Sumo la cant de nans para cada columna
    
    anio1 = str((df.date.values[0]))
    anio2 = str((df.date.values[-1]))
    time = pd.date_range(anio1, anio2, freq='D')
    
    f = f + len(time)-len(df) 
    f.code=df.code.notnull().sum()
    
    p = round((f*100)/len(time),2)
        
    return(p)
    


#%%

#Inicialización de valores constantes

#Directorio a los datos
DATOS='/home/dalia/Proyecto/BasesDatos/SMN/Exp186136.xlsx'

#Directorio a la carpeta con salidas
SALIDAS='/home/dalia/Proyecto/Salidas/SMN'

#%%

#Lectura de datos 
df = pd.read_excel(DATOS,1, usecols=[0,1,2,3,4,5])

#Renombro las columnas
df.columns=['station','date','max','min','pp','code']

#Cambio la clase de las fechas a datetime
df['date'] = pd.to_datetime(df['date']).dt.date

#La estacion 10034 cierra el 8/9/1998 y abre la estación 17970 en su lugar. 
#Entonces le cambio el numero de estación 10034 --> 17970 para facilitar el 
#analisis (como si fueran la misma estación)

df.loc[df.station==10034,"station"]= 17970

#%%

#Separo los datos por estación almacenandolos en un diccionario
estaciones = {"Rivadavia":df[df.station.eq(10006)],
              "Las Lomitas":df[df.station.eq(10011)],
              "Santiago del Estero":df[df.station.eq(10062)],
              "Roque Saenz Penia":df[df.station.eq(17970)],
              "Chamical":df[df.station.eq(10476)],
              "Chepes":df[df.station.eq(10102)],
              "Villa Dolores":df[df.station.eq(10117)],
              "Cordoba":df[df.station.eq(10100)],
              "Santa Rosa":df[df.station.eq(18030)],
              }

#No se usará la estación Comodoro Rivadavia, fuera del area de estudio
#%%

#Veo el periodo cubierto en cada estación
for key in estaciones.keys():
    i=estaciones[key].date.values[0]
    f=estaciones[key].date.values[-1]
    print("\n"+str(key)+": \n"+ 
          str(i)+ "\n"+
          str(f))

#%%

#Veo los faltantes por estación
for key in estaciones.keys():
    p=missing_per(estaciones[key])
    print("\n"+str(key))
    print(p)

#%%

#Ploteo los datos faltantes para ver si hay un periodo que no se pueda usar

for key in estaciones.keys():
    
    #Inicializo figura y ejes
    fig, ax = plt.subplots()
    
    #Inicializo una grilla para los subplots
    gs = gridspec.GridSpec(2,1)   
    
    #Grafico los datos faltantes de tmax, 0=NaN, 1!=NaN. 
    
    #Me posiciono en el primer subplot
    ax = plt.subplot(gs[0])
    #Inicializo un grafico de barras, eje x fecha, eje y datos de tmax como
    #booleans (True=1, False=0)
    bar = ax.bar(estaciones[key].date.values,
                 estaciones[key]['max'].astype(bool).values,
                 width=1)
    #Titulo
    ax.set_title( 'Max' , fontsize=12)
    
    #Grafico los datos faltantes de tmin
    
    #Me posiciono en el segundo subplot
    ax = plt.subplot(gs[1])
    bar = ax.bar(estaciones[key].date.values,
                 estaciones[key]['min'].astype(bool).values,
                 width=1)
    ax.set_title( 'Min' , fontsize=12)
    fig.suptitle("Faltantes "+str(key), fontsize=18)
    fig.subplots_adjust(left=0.05, bottom=0.05, right=0.98,top=0.85, wspace=0.5 ,hspace=0.35)
    
    fig.savefig(SALIDAS+"/Faltantes_"+str(key)+".png")

    plt.show()

#%%

#Mapa con la unbicación de las estaciones

#Lectura de la primer hoja del excel que contiene info de las estaciones
df_detalles = pd.read_excel(DATOS,0, usecols=[0,1,2,3,4])

#Reemplazo el nombre de las columnas
df_detalles.columns = ['station','name','lat','lon','height']

#Remuevo la fila de Comodoro Rivadavia porque no la voy a utilizar
df_detalles=df_detalles.drop(7,axis=0)

print(df_detalles.lat.min())
print(df_detalles.lat.max())
print(df_detalles.lon.min())
print(df_detalles.lon.max())


#Creo un diccionario, separando por estación el dataframe  
detalles = {"Rivadavia":df_detalles[df_detalles.station.eq(10006)],
            "Las Lomitas":df_detalles[df_detalles.station.eq(10011)],
            "Sgo.del Estero":df_detalles[df_detalles.station.eq(10062)],
            "Saenz Peña":df_detalles[df_detalles.station.eq(17970)],
            "Chamical":df_detalles[df_detalles.station.eq(10476)],
            "Chepes":df_detalles[df_detalles.station.eq(10102)],
            "Villa Dolores":df_detalles[df_detalles.station.eq(10117)],
            "Cordoba":df_detalles[df_detalles.station.eq(10100)],
            "Santa Rosa":df_detalles[df_detalles.station.eq(18030)],
    }

#%%
#Creamos la figura y definimos su tamaño
fig = plt.figure(figsize=(3.3, 5.3))

ax=plt.subplot(projection=ccrs.PlateCarree(central_longitude=180))

crs_latlon = ccrs.PlateCarree()
ax.set_extent([360-67,360-59, -33, -23], crs=crs_latlon)
ax.add_feature(cartopy.feature.COASTLINE)
ax.add_feature(cartopy.feature.BORDERS, linestyle='-', alpha=.5)

# Create a feature for States/Admin 1 regions at 1:50m from Natural Earth
states_provinces = cartopy.feature.NaturalEarthFeature(
     category='cultural',
     name='admin_1_states_provinces_lines',
     scale='10m',
     facecolor='none')

ax.add_feature(states_provinces, edgecolor='gray')

ax.set_xticks([295, 300], crs=crs_latlon)
ax.set_yticks([-35,-30, -25], crs=crs_latlon)
ax.grid(which='both', linewidth=0.3, linestyle='-')
ax.tick_params(axis='both', which='major', labelsize=6)
lon_formatter = LongitudeFormatter(zero_direction_label=True)
lat_formatter = LatitudeFormatter()
ax.xaxis.set_major_formatter(lon_formatter)
ax.yaxis.set_major_formatter(lat_formatter)

colors=['blue','green','red','yellow','k','orange','pink','cyan','magenta']

i=0
for key in detalles.keys():
    lat=detalles[key].lat.values*(-1)
    lon=360-detalles[key].lon.values
    
    ax.plot(lon,lat,color=colors[i],
            marker='o',markersize=6,
            transform=crs_latlon,label=str(key))
    i=i+1

ax.legend(bbox_to_anchor=(1.05, 1), loc='upper left', borderaxespad=0.)

fig.savefig(SALIDAS + '/Mapa.png', dpi=300, bbox_inches='tight')

#%%

#Análisis de precipitación

#Relleno los NaN con ceros.
for key in estaciones.keys():
   estaciones[key]['pp']=estaciones[key]['pp'].fillna(0)

#%%
#Códigos
#0 pp observada pero no medida --> NaN
#X no se midió -> acumulada dia posterior --> NaN
#A cantidad acumulada
#F cantidad diaria incompleta 
#SF dato faltante --> NaN
#- no se realizó observación --> NaN

for key in estaciones.keys():
    index=estaciones[key][(estaciones[key]['code']=='0')  | 
                          (estaciones[key]['code']=='X')  | 
                          (estaciones[key]['code']=='SF') |
                          (estaciones[key]['code']=='-')].index.to_list()
    estaciones[key]['pp'][index]=np.NaN

#%%
for key in estaciones.keys():
    print(key)
    f=estaciones[key]['pp'].isnull().sum()/len(estaciones[key])
    print(f)
    
