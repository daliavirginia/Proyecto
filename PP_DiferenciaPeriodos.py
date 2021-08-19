# -*- coding: utf-8 -*-
"""
Created on Fri Aug  6 23:23:10 2021

@author: Dalia
"""
#%%
#DIFERENCIA DOS PERÍODOS:(1991-2020) - (1940-1970)

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

# Constantes (en mayuscula)
LAT_SUR=-60
LAT_NOR=-20
LON_ESTE=320
LON_OESTE=280

# Periodo 1
PER1_ANIO_MIN=1940
PER1_ANIO_MAX=1970

TIEMPO1_MIN = str(PER1_ANIO_MIN)+'-01-01'
TIEMPO1_MAX= str(PER1_ANIO_MAX)+'-12-31'

# Periodo 2
PER2_ANIO_MIN=1986
PER2_ANIO_MAX=2016

TIEMPO2_MIN=str(PER1_ANIO_MIN)+'-01-01'
TIEMPO2_MAX= str(PER2_ANIO_MAX)+'-12-31'

# Directorio al archivo
DATOS='/home/dalia/Proyecto/BasesDatos/precip.mon.total.v2018.nc'

# Directorio a salidas
SALIDAS='/home/dalia/Proyecto/Salidas/'


#%%

# Abro el archivo usando xarray
data=xr.open_dataset(DATOS)

# Hago recorte para el periodo 1
ds_per1 = data.sel(lat=slice(LAT_NOR, LAT_SUR), lon=slice(LON_OESTE,LON_ESTE),
                        time=slice(TIEMPO1_MIN,TIEMPO1_MAX))
# Hago recorte para el período 2
ds_per2 = data.sel(lat=slice(LAT_NOR, LAT_SUR), lon=slice(LON_OESTE,LON_ESTE),
                        time=slice(TIEMPO2_MIN,TIEMPO2_MAX))

#Remuestreo los datos dividiendo por estaciones (DEF-MAM-JJA-SON)
ds_per1_est = ds_per1.resample(time="QS-DEC").sum()
ds_per2_est = ds_per2.resample(time="QS-DEC").sum()

#Divido en cuartos (no me gusta esta forma de separar, muchas variables
#capaz se puede mejorar)
ds_djf_per1=ds_per1_est.sel(time=ds_per1_est["time.month"]==12)
ds_mam_per1=ds_per1_est.sel(time=ds_per1_est["time.month"]==3)
ds_jja_per1=ds_per1_est.sel(time=ds_per1_est["time.month"]==6)
ds_son_per1=ds_per1_est.sel(time=ds_per1_est["time.month"]==9)

ds_djf_per2=ds_per2_est.sel(time=ds_per2_est["time.month"]==12)
ds_mam_per2=ds_per2_est.sel(time=ds_per2_est["time.month"]==3)
ds_jja_per2=ds_per2_est.sel(time=ds_per2_est["time.month"]==6)
ds_son_per2=ds_per2_est.sel(time=ds_per2_est["time.month"]==9)

#De los datos de verano me quedo desde el segundo hasta el anteultimo ya que
#se agregaron meses vacios para poder completar el primer y ultimo cuarto.
ds_djf_per1=ds_djf_per1.sel(time=slice(ds_djf_per1["time"][1],ds_djf_per1["time"][-2]))
ds_djf_per2=ds_djf_per2.sel(time=slice(ds_djf_per2["time"][1],ds_djf_per2["time"][-2]))



#%%

#Resto las medias estacionales

ds_djf_dif=ds_djf_per2.mean("time")-ds_djf_per1.mean("time")
ds_mam_dif=ds_mam_per2.mean("time")-ds_mam_per1.mean("time")
ds_jja_dif=ds_jja_per2.mean("time")-ds_jja_per1.mean("time")
ds_son_dif=ds_son_per2.mean("time")-ds_son_per1.mean("time")

#%%

#Grafico

fig, ax = plt.subplots(figsize=(2*3.3,2*4))

#Define grid for subplots
gs = gridspec.GridSpec(2,2) 

#latitudes and longitudes to plot
lons, lats = np.meshgrid(ds_per1['lon'], ds_per1['lat'])

#Definimos los niveles para los contornos (inferior, superior, longitud de los intervalos)
#Si utilizamos las funciones ds.min ds.max podemos darnos una idea de cómo definirlos
lev_inf = 0
lev_sup = 1200
lev_int = 100
clevs = np.arange(lev_inf, lev_sup+lev_int, lev_int)

#Plot DJF mean

ax = plt.subplot(gs[0],projection=ccrs.PlateCarree(central_longitude=180))

crs_latlon = ccrs.PlateCarree()
ax.set_extent([LON_OESTE, LON_ESTE, LAT_NOR, LAT_SUR], crs=crs_latlon)
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
im=ax.contourf(lons, lats, ds_djf_dif["precip"],cmap="RdYlGn",levels=np.arange(-70,70,10),extend='both',transform=crs_latlon)

ax.set_xticks(np.arange(LON_OESTE,LON_ESTE,10), crs=crs_latlon)
ax.set_yticks(np.arange(LAT_SUR,LAT_NOR,10), crs=crs_latlon)
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
ax.set_extent([LON_OESTE, LON_ESTE, LAT_NOR, LAT_SUR], crs=crs_latlon)
ax.add_feature(cartopy.feature.COASTLINE)
ax.add_feature(cartopy.feature.BORDERS, linestyle='-', alpha=.5)

# Create a feature for States/Admin 1 regions at 1:50m from Natural Earth
states_provinces = cartopy.feature.NaturalEarthFeature(
    category='cultural',
    name='admin_1_states_provinces_lines',
    scale='10m',
    facecolor='none')

ax.add_feature(states_provinces, edgecolor='gray')

im=ax.contourf(lons, lats, ds_mam_dif["precip"],cmap="RdYlGn",levels=np.arange(-70,70,10),extend='both',transform=crs_latlon)

ax.set_xticks(np.arange(LON_OESTE,LON_ESTE,10), crs=crs_latlon)
ax.set_yticks(np.arange(LAT_SUR,LAT_NOR,10), crs=crs_latlon)
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
ax.set_extent([LON_OESTE, LON_ESTE, LAT_NOR, LAT_SUR], crs=crs_latlon)
ax.add_feature(cartopy.feature.COASTLINE)
ax.add_feature(cartopy.feature.BORDERS, linestyle='-', alpha=.5)

# Create a feature for States/Admin 1 regions at 1:50m from Natural Earth
states_provinces = cartopy.feature.NaturalEarthFeature(
    category='cultural',
    name='admin_1_states_provinces_lines',
    scale='10m',
    facecolor='none')

ax.add_feature(states_provinces, edgecolor='gray')

im=ax.contourf(lons, lats, ds_jja_dif["precip"],cmap="RdYlGn",levels=np.arange(-70,70,10),extend='both',transform=crs_latlon)

ax.set_xticks(np.arange(LON_OESTE,LON_ESTE,10), crs=crs_latlon)
ax.set_yticks(np.arange(LAT_SUR,LAT_NOR,10), crs=crs_latlon)
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
ax.set_extent([LON_OESTE, LON_ESTE, LAT_NOR, LAT_SUR], crs=crs_latlon)
ax.add_feature(cartopy.feature.COASTLINE)
ax.add_feature(cartopy.feature.BORDERS, linestyle='-', alpha=.5)
# Create a feature for States/Admin 1 regions at 1:50m from Natural Earth
states_provinces = cartopy.feature.NaturalEarthFeature(
    category='cultural',
    name='admin_1_states_provinces_lines',
    scale='10m',
    facecolor='none')

ax.add_feature(states_provinces, edgecolor='gray')

im=ax.contourf(lons, lats, ds_son_dif["precip"],cmap="RdYlGn",levels=np.arange(-70,70,10),extend='both',transform=crs_latlon)

ax.set_xticks(np.arange(LON_OESTE,LON_ESTE,10), crs=crs_latlon)
ax.set_yticks(np.arange(LAT_SUR,LAT_NOR,10), crs=crs_latlon)
ax.grid(which='both', linewidth=0.3, linestyle='-')
ax.tick_params(axis='both', which='major', labelsize=10)
lon_formatter = LongitudeFormatter(zero_direction_label=True)
lat_formatter = LatitudeFormatter()
ax.xaxis.set_major_formatter(lon_formatter)
ax.yaxis.set_major_formatter(lat_formatter)    
ax.set_title('Primavera (SON)')

fig.suptitle("Diferencia de precipitación acumulada\n entre periodos (" + str(PER2_ANIO_MIN)+"-"+ 
             str(PER2_ANIO_MAX) + ") - (" + str(PER1_ANIO_MIN)+
             "-"+str(PER1_ANIO_MAX)+")",fontsize=16)


#fig.subplots_adjust(left=0.02, bottom=0.03, right=0.8,top=0.9, wspace=0.04 ,hspace=0.3)

#Plot colorbar
cbar_ax = fig.add_axes([0.98, 0.1, 0.03, 0.8])
cbar_ax.set_title('mm')
fig.colorbar(im, cax=cbar_ax,orientation='vertical')


fig.savefig(SALIDAS + 'DiferenciaPeriodos.png', dpi=300, bbox_inches='tight')

fig.show()


