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

#%%

def missing_per(df):
    
    f = df.isnull().sum()
    anio1 = str(df.index[0])
    anio2 = str(df.index[len(df["day"])-1])
    time = pd.date_range(anio1+"-01-01", anio2+"-12-31")
    f = f + len(time)-len(df["day"]) 
    p = (f*100)/len(time)
    
    return(p)
    

#%%

 # Directorio al archivo
DATOS='/home/dalia/Proyecto/BasesDatos/SMN_DCAO/Leandro_dd_SMN.dat'

# Directorio a salidas
SALIDAS='/home/dalia/Proyecto/Salidas/'

#%%
#Leo el archivo .dat con pandas

data = pd.read_table(DATOS, sep='\s+', header=None, usecols=[0,1,2,3,4,5,6])

data.columns = ["station","day","month","year","tmax","tmin","pp"]


data2 = data.set_index(["station","year"]).sort_index()

print(data2.index.name)
print(data2.index.values)
print(data2.index.levels)


#87065:Rivadavia, Salta
#87078:Las Lomitas, Formosa
#87129:Santiago del Estero aero
#87148:Roque Saenz Peña, Chaco
#87320:Chamical aero, La Rioja
#87322:Chepes, La Rioja
#87328:Villa Dolores, Cordoba
#87344:Córdoba aero
#87444:Sta Rosa de Conlara aero, San Luis


rivadavia=data2.loc[(87065),:]
lomitas=data2.loc[(87078),:]
santiago=data2.loc[(87129),:]
saenz_penia=data2.loc[(87148),:]
chamical=data2.loc[(87320),:]
chepes=data2.loc[(87322),:]
vdolores=data2.loc[(87328),:]
cordoba_aero=data2.loc[(87344),:]
sta_rosa=data2.loc[(87444),:]

#Para hacer el netcdf:
    #vars: tmax,tmin,pp
    #dims: time
#Pruebo hacer uno con los datos de Rivadavia

 
time = pd.date_range("1959-01-01", "2020-12-31") 

### OJO FALTA LEER LA COLUMNA DE LAS LETRITAS #############

    
    

#%%

missing_per(lomitas)
missing_per(rivadavia) #chequear
missing_per(santiago)
missing_per(saenz_penia) # chequear
missing_per(chamical)
missing_per(chepes) #chequear
missing_per(vdolores) #chequear
missing_per(cordoba_aero)
missing_per(sta_rosa) #chequear 

# Usar Lomitas, Santiago, Chamical, 
#Cordoba aero , Villa Dolroes