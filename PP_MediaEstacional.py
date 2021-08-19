#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Aug  3 11:31:30 2021

@author: Dalia
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

# Declaro variables

# Dominio
lat_sur=-60
lat_nor=-20
lon_este=320
lon_oeste=280

# 1950-2014 Available period
anio_min=1980
anio_max=2010

tiempo_min = str(anio_min)+'-01-01'
tiempo_max = str(anio_max)+'-12-31'


# Directorio al archivo
file='/home/dalia/Proyecto/BasesDatos/precip.mon.total.v2018.nc'

# Directorio a salidas
salidas='/home/dalia/Proyecto/Salidas/'

#%%

# Abro el archivo usando xarray
data=xr.open_dataset(file)

print(data["precip"]["time"])
# Imprimo datos de la variable
print(data["precip"]) # Precipitación mensual total

# Hago recorte 
data_recorte = data.sel(lat=slice(lat_nor, lat_sur), lon=slice(lon_oeste,lon_este),
                        time=slice(tiempo_min,tiempo_max))

# Redimensiono el array sumando trimestralmente
data_recorte_est = data_recorte.resample(time='QS-DEC').sum(skipna=False)
# Ahora cada trimestre es representado en la dimensión temporal con el primer 
# día del primer mes
# El string "QS_DEC" significa, Q: cuartos (dividir el año en cuartos),
# S: start, -DEC: el año empieza en diciembre (sin el S de start sería que el
# año termina en diciembre).

# Para poder completar los trimestres DEF, resample agrega un diciembre delante
# de la serie temporal, y un enero y febrero por detrás. Estos meses agregados
# tienen pp=0. Saco el EF del primer año y el

#Select summer accumulated (DJF season)
ds_djf = data_recorte_est.sel(time=data_recorte_est['time.month']==12)
#Exclude JF of first year and december of the last
ds_djf=ds_djf.sel(time=slice(ds_djf['time'][1],ds_djf['time'][-2]))
# Desde el segundo DEF hasta el  anteultimo DEF
#Select autumn accumulated (MAM season)
ds_mam = data_recorte_est.sel(time=data_recorte_est['time.month']==3)   
#Select winter accumulated (JJA season)
ds_jja = data_recorte_est.sel(time=data_recorte_est['time.month']==6)   
#Select spring accumulated (SON season)
ds_son = data_recorte_est.sel(time=data_recorte_est['time.month']==9) 

#Promedios

ds_djf_mean=ds_djf.mean("time")
ds_mam_mean=ds_mam.mean("time")
ds_jja_mean=ds_jja.mean("time")
ds_son_mean=ds_son.mean("time")

#%%

#Grafico

fig, ax = plt.subplots(figsize=(2*3.3,2*3.5))

#Define grid for subplots
gs = gridspec.GridSpec(2,2) 

#latitudes and longitudes to plot
lons, lats = np.meshgrid(data_recorte['lon'], data_recorte['lat'])

#Definimos los niveles para los contornos (inferior, superior, longitud de los intervalos)
#Si utilizamos las funciones ds.min ds.max podemos darnos una idea de cómo definirlos
lev_inf = 0
lev_sup = 1200
lev_int = 100
clevs = np.arange(lev_inf, lev_sup+lev_int, lev_int)

#Plot DJF mean

ax = plt.subplot(gs[0],projection=ccrs.PlateCarree(central_longitude=180))

crs_latlon = ccrs.PlateCarree()
ax.set_extent([lon_oeste, lon_este, lat_nor, lat_sur], crs=crs_latlon)
ax.add_feature(cartopy.feature.COASTLINE)
ax.add_feature(cartopy.feature.BORDERS, linestyle='-', alpha=.5)

# Create a feature for States/Admin 1 regions at 1:50m from Natural Earth
states_provinces = cartopy.feature.NaturalEarthFeature(
    category='cultural',
    name='admin_1_states_provinces_lines',
    scale='10m',
    facecolor='none')

ax.add_feature(states_provinces, edgecolor='gray')

#Las medias estacionales son Datasets y tengo que introducir un DataArray
#Por eso selecciono la variable ["precip"]. 
im=ax.contourf(lons, lats, ds_djf_mean["precip"],cmap="winter_r",extend='both',transform=crs_latlon)

ax.set_xticks(np.arange(lon_oeste,lon_este,10), crs=crs_latlon)
ax.set_yticks(np.arange(lat_sur,lat_nor,10), crs=crs_latlon)
ax.grid(which='both', linewidth=0.3, linestyle='-')
ax.tick_params(axis='both', which='major', labelsize=10)
lon_formatter = LongitudeFormatter(zero_direction_label=True)
lat_formatter = LatitudeFormatter()
ax.xaxis.set_major_formatter(lon_formatter)
ax.yaxis.set_major_formatter(lat_formatter)    

ax.set_title('Verano (DEF)')

#MAM

ax = plt.subplot(gs[1],projection=ccrs.PlateCarree(central_longitude=180))

crs_latlon = ccrs.PlateCarree()
ax.set_extent([lon_oeste, lon_este, lat_nor, lat_sur], crs=crs_latlon)
ax.add_feature(cartopy.feature.COASTLINE)
ax.add_feature(cartopy.feature.BORDERS, linestyle='-', alpha=.5)

# Create a feature for States/Admin 1 regions at 1:50m from Natural Earth
states_provinces = cartopy.feature.NaturalEarthFeature(
    category='cultural',
    name='admin_1_states_provinces_lines',
    scale='10m',
    facecolor='none')

ax.add_feature(states_provinces, edgecolor='gray')


im=ax.contourf(lons, lats, ds_mam_mean["precip"],cmap="winter_r",extend='both',transform=crs_latlon)

ax.set_xticks(np.arange(lon_oeste,lon_este,10), crs=crs_latlon)
ax.set_yticks(np.arange(lat_sur,lat_nor,10), crs=crs_latlon)
ax.grid(which='both', linewidth=0.3, linestyle='-')
ax.tick_params(axis='both', which='major', labelsize=10)
lon_formatter = LongitudeFormatter(zero_direction_label=True)
lat_formatter = LatitudeFormatter()
ax.xaxis.set_major_formatter(lon_formatter)
ax.yaxis.set_major_formatter(lat_formatter)    

ax.set_title('Otoño (MAM)')

#JJA

ax = plt.subplot(gs[2],projection=ccrs.PlateCarree(central_longitude=180))

crs_latlon = ccrs.PlateCarree()
ax.set_extent([lon_oeste, lon_este, lat_nor, lat_sur], crs=crs_latlon)
ax.add_feature(cartopy.feature.COASTLINE)
ax.add_feature(cartopy.feature.BORDERS, linestyle='-', alpha=.5)

# Create a feature for States/Admin 1 regions at 1:50m from Natural Earth
states_provinces = cartopy.feature.NaturalEarthFeature(
    category='cultural',
    name='admin_1_states_provinces_lines',
    scale='10m',
    facecolor='none')

ax.add_feature(states_provinces, edgecolor='gray')


im=ax.contourf(lons, lats, ds_jja_mean["precip"],cmap="winter_r",extend='both',transform=crs_latlon)

ax.set_xticks(np.arange(lon_oeste,lon_este,10), crs=crs_latlon)
ax.set_yticks(np.arange(lat_sur,lat_nor,10), crs=crs_latlon)
ax.grid(which='both', linewidth=0.3, linestyle='-')
ax.tick_params(axis='both', which='major', labelsize=10)
lon_formatter = LongitudeFormatter(zero_direction_label=True)
lat_formatter = LatitudeFormatter()
ax.xaxis.set_major_formatter(lon_formatter)
ax.yaxis.set_major_formatter(lat_formatter)    

ax.set_title('Invierno (JJA)')

#SON

ax = plt.subplot(gs[3],projection=ccrs.PlateCarree(central_longitude=180))

crs_latlon = ccrs.PlateCarree()
ax.set_extent([lon_oeste, lon_este, lat_nor, lat_sur], crs=crs_latlon)
ax.add_feature(cartopy.feature.COASTLINE)
ax.add_feature(cartopy.feature.BORDERS, linestyle='-', alpha=.5)

# Create a feature for States/Admin 1 regions at 1:50m from Natural Earth
states_provinces = cartopy.feature.NaturalEarthFeature(
    category='cultural',
    name='admin_1_states_provinces_lines',
    scale='10m',
    facecolor='none')

ax.add_feature(states_provinces, edgecolor='gray')

im=ax.contourf(lons, lats, ds_son_mean["precip"],cmap="winter_r",extend='both',transform=crs_latlon)

ax.set_xticks(np.arange(lon_oeste,lon_este,10), crs=crs_latlon)
ax.set_yticks(np.arange(lat_sur,lat_nor,10), crs=crs_latlon)
ax.grid(which='both', linewidth=0.3, linestyle='-')
ax.tick_params(axis='both', which='major', labelsize=10)
lon_formatter = LongitudeFormatter(zero_direction_label=True)
lat_formatter = LatitudeFormatter()
ax.xaxis.set_major_formatter(lon_formatter)
ax.yaxis.set_major_formatter(lat_formatter)    
ax.set_title('Primavera (SON)')

fig.suptitle("Precipitación acumulada estacional media (" + str(anio_min)+"-"+ 
             str(anio_max) + ")\n GPCC",fontsize=16)


#fig.subplots_adjust(left=0.02, bottom=0.03, right=0.8,top=0.9, wspace=0.04 ,hspace=0.3)

#Plot colorbar
cbar_ax = fig.add_axes([0.98, 0.1, 0.03, 0.8])
cbar_ax.set_title('mm')
fig.colorbar(im, cax=cbar_ax,orientation='vertical')

fig.savefig(salidas + 'MediaEstacionalGPCC_.png', dpi=300, bbox_inches='tight')

fig.show()


