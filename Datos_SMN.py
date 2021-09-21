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
#%%

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

DATOS='/home/dalia/Proyecto/BasesDatos/SMN/Exp186136.xlsx'
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

#%%
#Veo el periodo cubierto en cada estación

for key in estaciones.keys():
    i=estaciones[key].date.values[0]
    f=estaciones[key].date.values[-1]
    print("\n"+str(key)+": \n"+ 
          str(i)+ "\n"+
          str(f))

#%%

for key in estaciones.keys():
    p=missing_per(estaciones[key])
    print("\n"+str(key))
    print(p)

#%%

#Ploteo los datos faltantes para ver si hay un periodo que no se pueda usar

for key in estaciones.keys():
    
    
    fig, ax = plt.subplots()
    #Define grid for subplots
    gs = gridspec.GridSpec(2,1)   
    
    
    ax = plt.subplot(gs[0])
    bar = ax.bar(estaciones[key].date.values,
                 estaciones[key]['max'].astype(bool).values,
                 width=0)
    ax.set_title( 'Max' , fontsize=12)
    
     
    ax = plt.subplot(gs[1])
    bar = ax.bar(estaciones[key].date.values,
                 estaciones[key]['min'].astype(bool).values
                 width=0)
    ax.set_title( 'Min' , fontsize=12)
    fig.suptitle("Faltantes "+str(key), fontsize=18)
    fig.subplots_adjust(left=0.03, bottom=0.03, right=0.93,top=0.96, wspace=0.5 ,hspace=0.35)
    
    fig.savefig(SALIDAS+"/Faltantes_"+str(key)+".png")
    
    
     
    plt.show()

