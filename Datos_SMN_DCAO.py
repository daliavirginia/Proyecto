#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Aug  8 22:03:43 2021

@author: dalia
"""

#%%
import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
import datetime as dt
import xarray as xr
import netCDF4 as nc
import xlrd 

#%%

def missing_per(df):
    
    f = df.isnull().sum() #Sumo la cant de nans para cada columna
    
    anio1 = str(int(df.year[0]))
    anio2 = str(int(df.year[len(df)-1]))
    time = pd.date_range(anio1+"-01-01", anio2+"-12-31")
    
    f = f + len(time)-len(df) 
    f.code=df.code.notnull().sum()
    
    p = round((f*100)/len(time),2)
        
    return(p)
    

#%%

# Directorio al archivo
DATOS='/home/dalia/Proyecto/BasesDatos/SMN_DCAO/Leandro_dd_SMN.dat'
DATOS_BRENAS='/home/dalia/Proyecto/BasesDatos/Estaciones/NH0416.csv'
DATOS_VDOLORES='/home/dalia/Proyecto/BasesDatos/Estaciones/Exp.178656.xlsx'

# Directorio a salidas
SALIDAS='/home/dalia/Proyecto/Salidas/'

#%%
#Leo el archivo .dat con pandas
data = pd.read_fwf(DATOS, header=None, widths=[6,4,4,6,6,6,6,6])


vd = pd.read_excel(DATOS_VDOLORES,1)
#Le pongo nombre a las columnas
data.columns = ["station","day","month","year","tmax","tmin","pp","code"]
#%%

vd.drop(vd.columns[4:9],axis=1,inplace=True)
vd.columns=["date","tmax","tmin","pp"]
vd['date']=pd.to_datetime(vd['date'])
vd['code']=np.nan
vd['year']=np.nan
vd['pp'][vd.pp.isnull()]=0

for i in range(len(vd)):
    vd['year'][i]=int(vd.date[i].year)

missing_per(vd)


#%%

# Uso la función eq() para filtrar los datos por estacion y guardarlos
#en un diccionario. 
estaciones = {"Rivadavia":data[data.station.eq(87065)],
              "LasLomitas":data[data.station.eq(87078)],
              "Santiago_aero":data[data.station.eq(87129)],
              "SaenzPenia":data[data.station.eq(87148)],
              "Chamical":data[data.station.eq(87322)],
              "Chepes":data[data.station.eq(87322)],
              "VillaDolores":data[data.station.eq(87328)],
              "VillaDolores2":vd,
              "Cordoba_aero":data[data.station.eq(87344)],
              "StaRosa_aero":data[data.station.eq(87444)]
              }
#%%

# Vuelvo a indexar cada dataframe para que empiece en 0 ya que arrastró
#los indices del dataframe original
for key in estaciones.keys():
    estaciones[key].index=range(len(estaciones[key]))
    
#%%
for key in estaciones.keys():
   print(str(key)+"\n"+str(missing_per(estaciones[key]))+"\n")

#%%
#87065:Rivadavia, Salta
#87078:Las Lomitas, Formosa
#87129:Santiago del Estero aero
#87148:Roque Saenz Peña, Chaco
#87320:Chamical aero, La Rioja
#87322:Chepes, La Rioja
#87328:Villa Dolores, Cordoba
#87344:Córdoba aero
#87444:Sta Rosa de Conlara aero, San Luis

#%%
#CODIGOS
#0 pp observada pero no medida
#X no se midio
#A cantidad acumulada
#F cantidad diaria incompleta
#SF synop faltante
#T dato telegrafico
#PF proveniente de faja
#- no se realizó observacion

#%%
   


