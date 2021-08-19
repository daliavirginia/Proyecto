#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Aug 19 14:33:36 2021

@author: dalia
"""

#%%

#Load requested libraries
import numpy as np   # Manejar arrays
import xarray as xr  # Manejar arrays (especializado nc)
from scipy import stats # LLamo a la fun stats de scipy 
import netCDF4
import pandas as pd

from matplotlib import pyplot as plt #Graficar (ploteos basicos)
from matplotlib import gridspec 
import matplotlib

import cartopy.crs as ccrs	# Graficar mapas 
import cartopy.feature 	# Notacion de punto (escribir todo)
from cartopy.mpl.ticker import LongitudeFormatter, LatitudeFormatter

#%%

# Declaro variables

# Constantes (en mayuscula)
LAT_SUR=-60
LAT_NOR=-20
LON_ESTE=-40
LON_OESTE=-80

# Periodo 1
PER1_ANIO_MIN=1940
PER1_ANIO_MAX=1970

TIEMPO1_MIN = str(PER1_ANIO_MIN)+'-01-01'
TIEMPO1_MAX= str(PER1_ANIO_MAX)+'-12-31'

# Periodo 2
PER2_ANIO_MIN=1980
PER2_ANIO_MAX=2010

TIEMPO2_MIN=str(PER1_ANIO_MIN)+'-01-01'
TIEMPO2_MAX= str(PER2_ANIO_MAX)+'-12-31'

# Directorio al archivo
DATOS='/home/dalia/Proyecto/BasesDatos/inLand_gld_LU_CO2_SA_50_10_E.nc'

# Directorio a salidas
SALIDAS='/home/dalia/Proyecto/Salidas/'

#%%

# Abro el archivo usando xarray
data=xr.open_dataset(DATOS, decode_times=False)

# months y years no es una unidad que pueda decodificar

#extraigo de "time" el atributo unidad y que lo separe en dos en "since"
units, reference_date = data.time.attrs['units'].split('since')

#Con la funcion date_range de pandas armo el vector de fechas, usando
#la info que extraje antes
data['time'] = pd.date_range(start=reference_date, periods=data.sizes['time'], freq='MS')

# Hago recorte para el periodo 1
data_per1 = data.sel(latitude=slice(LAT_NOR, LAT_SUR), longitude=slice(LON_OESTE,LON_ESTE),
                        time=slice(TIEMPO1_MIN,TIEMPO1_MAX))
# Hago recorte para el período 2
data_per2 = data.sel(latitude=slice(LAT_NOR, LAT_SUR), longitude=slice(LON_OESTE,LON_ESTE),
                        time=slice(TIEMPO2_MIN,TIEMPO2_MAX))


