## Librerías
require(ncdf4)
require(lubridate)
require(here)
require(RColorBrewer)
require(mapdata)
require(maps)
require(scales)
require(ggplot2)


## Carga de datos al environment

archivo = here("data","CMAP_precip.mon.mean.nc")
nc = nc_open(archivo)

## Extraigo info
# Nombre de las variables
nombre_var = names(nc$var)
# Nombre de las dimensiones
nombres_dim = names(nc$dim)
# Extraigo datos de las dimensiones (lat,long,tiempo) 
datos_dim = 0
for (i in 1:length(nombres_dim)) {
  datos_dim[i] = list(ncvar_get(nc, nombres_dim[i]))
}
names(datos_dim) = nombres_dim

## Extraigo datos de pp
# Hago un recorte del dominio . Me centro en el sur de sudamerica.

lats = which((datos_dim$lat < (-20)) & (datos_dim$lat >(-60)))
lons = which((datos_dim$lon < 320) & (datos_dim$lon > 280))

## Ahora si extraigo datos de la pp
variable = ncvar_get(nc, "precip", start = c(lons[1],lats[1], 1), count = c(length(lons),length(lats), -1))

## Array de 3 dimensiones latxlonxtiempo

# Ahora quiero hacer promedios estacionales para el (1981-2010). Tengo que promediar , mar+abr+may,
# jun+jul+ago y sep+oct+nov. 

# Tiempo

ncatt_get(nc, "time", "units")

fechas = as.Date(datos_dim$time/24, origin = "1800-01-01")

## Extraigo los datos de 1981 - 2010

anios = as.integer(substr(fechas,1,4)) #Vector con los años como enteros
variable_recorte = variable[,, which((anios>1980) & (anios<2011))] #Uso el vector de anios para hacer el recorte
fechas = fechas[which((anios>1980) & (anios<2011))]

## Ahora hago variables por estacion 

pp_verano = variable_recorte[,,(months(fechas)== "enero"
                                | months(fechas)== "diciembre"
                                | months(fechas)== "febrero")]
pp_otonio = variable_recorte[,,(months(fechas)== "marzo"
                                | months(fechas)== "abril"
                                | months(fechas)== "mayo")]
pp_invierno = variable_recorte[,,(months(fechas)== "junio"
                                  | months(fechas)== "julio"
                                  | months(fechas)== "agosto")]
pp_primavera = variable_recorte[,,(months(fechas)== "septiembre"
                                  | months(fechas)== "octubre"
                                  | months(fechas)== "noviembre")]
# Calculo las medias
media_pp_DEF = apply(pp_verano, c(1,2), mean)
media_pp_MAM = apply(pp_otonio, c(1,2), mean)
media_pp_JJA = apply(pp_invierno, c(1,2), mean)
media_pp_SON = apply(pp_primavera, c(1,2), mean)

# Grafico 

## Para usar ggplot2 necesito generar data frames con los datos

df_DEF = data.frame(x=rep(datos_dim$lon[lons], length(lats)),
                    y=rep(datos_dim$lat[lats], each=length(lons)),
                    z=array(media_pp_DEF, length(lats)*length(lons)))

df_MAM = data.frame(x=rep(datos_dim$lon[lons], length(lats)),
                    y=rep(datos_dim$lat[lats], each=length(lons)),
                    z=array(media_pp_MAM, length(lats)*length(lons)))

df_JJA = data.frame(x=rep(datos_dim$lon[lons], length(lats)),
                    y=rep(datos_dim$lat[lats], each=length(lons)),
                    z=array(media_pp_JJA, length(lats)*length(lons)))

df_SON = data.frame(x=rep(datos_dim$lon[lons], length(lats)),
                    y=rep(datos_dim$lat[lats], each=length(lons)),
                    z=array(media_pp_SON, length(lats)*length(lons)))

## Ajustes estéticos 
map.world <- map_data("world2")
my_fill <- scale_fill_distiller(palette='GnBu', direction = 1,limits = c(0,8),breaks=pretty_breaks(10),labs(fill = "mm/día"))
my_theme <- theme_bw() + theme(panel.ontop=TRUE,
                               panel.background=element_blank(),
                               plot.title = element_text(hjust = 0.5),panel.grid.major = element_blank(),
                               panel.grid.minor = element_blank())
## DEF

p = ggplot(df_DEF, aes(x=x, y=y)) + geom_tile(aes(fill=z)) + my_fill +
  ggtitle("Climatología precipitación DEF (1981-2010)")

p = p + geom_map(dat=map.world, map = map.world, aes(map_id=region), xlim=range(df_DEF$x),
                 ylim=range(df_DEF$y), fill="NA", color="black", inherit.aes = F) +
  coord_quickmap(xlim = range(df_DEF$x), ylim = range(df_DEF$y), expand = FALSE) +
  scale_x_continuous(limits = range(df_DEF$x)) + scale_y_continuous(limits = range(df_DEF$y)) +
  xlab("longitud") +
  ylab("latitud") + my_theme

ggsave(here("salidas","pp DEF.png"))

## MAM

p = ggplot(df_MAM, aes(x=x, y=y)) + geom_tile(aes(fill=z)) + my_fill +
  ggtitle("Climatología precipitación MAM (1981-2010)")

p = p + geom_map(dat=map.world, map = map.world, aes(map_id=region), xlim=range(df_MAM$x),
                 ylim=range(df_MAM$y), fill="NA", color="black", inherit.aes = F) +
  coord_quickmap(xlim = range(df_MAM$x), ylim = range(df_MAM$y), expand = FALSE) +
  scale_x_continuous(limits = range(df_MAM$x)) + scale_y_continuous(limits = range(df_DEF$y)) +
  xlab("longitud") +
  ylab("latitud") + my_theme

ggsave(here("salidas","pp MAM.png"))

## JJA

p = ggplot(df_JJA, aes(x=x, y=y)) + geom_tile(aes(fill=z)) + my_fill +
  ggtitle("Climatología precipitación JJA (1981-2010)")

p = p + geom_map(dat=map.world, map = map.world, aes(map_id=region), xlim=range(df_JJA$x),
                 ylim=range(df_JJA$y), fill="NA", color="black", inherit.aes = F) +
  coord_quickmap(xlim = range(df_JJA$x), ylim = range(df_JJA$y), expand = FALSE) +
  scale_x_continuous(limits = range(df_JJA$x)) + scale_y_continuous(limits = range(df_DEF$y)) +
  xlab("longitud") +
  ylab("latitud") + my_theme

ggsave(here("salidas","pp JJA.png"))


## SON

p = ggplot(df_SON, aes(x=x, y=y)) + geom_tile(aes(fill=z)) + my_fill +
  ggtitle("Climatología precipitación SON (1981-2010)")

p = p + geom_map(dat=map.world, map = map.world, aes(map_id=region), xlim=range(df_SON$x),
                 ylim=range(df_SON$y), fill="NA", color="black", inherit.aes = F) +
  coord_quickmap(xlim = range(df_SON$x), ylim = range(df_SON$y), expand = FALSE) +
  scale_x_continuous(limits = range(df_SON$x)) + scale_y_continuous(limits = range(df_DEF$y)) +
  xlab("longitud") +
  ylab("latitud") + my_theme

ggsave(here("salidas","pp SON.png"))

rm(media_pp_DEF,media_pp_JJA,media_pp_MAM,media_pp_SON,variable_recorte,z)

################################################################################

# Cambios entre 2 períodos

## Recorte de datos

var_per1 = variable[,, which((anios>1990) & (anios<2021))]
var_per2 = variable[,, which((anios>1939) & (anios<1971))]

## promedios y resta

media_per1 = apply(var_per1, c(1,2), mean)
media_per2 = apply(var_per2, c(1,2), mean)

resta = media_per1 - media_per2

################################################################################

# Centro del pais

lats = which((datos_dim$lat < (-25)) & (datos_dim$lat >(-37)))
lons = which((datos_dim$lon < 330) & (datos_dim$lon > 300))

## Extraigo datos variable

variable = ncvar_get(nc, "precip", start = c(lons[1],lats[1], 1), count = c(length(lons),length(lats), -1))

## Promedio espacial

prom_espacial = apply(variable, 3, mean)

df1=data.frame()


## Dimension tiempo

fechas = as.Date(datos_dim$time/24, origin = "1800-01-01") # Objeto clase Date
anios = substr(fechas,1,4) # Vector con los años y meses como integers

## Quiero hacer el acumulado anual, quito los ultimos 10 datos

prom_espacial = prom_espacial[which(anios<2019)]

## Redimensiono el array

prom_espacial= array(prom_espacial, dim = c(12,40))
anios=1979:2018

## Hago el acumulado anual

acumulado = apply(prom_espacial, 2, sum)

## Armo el data frame para ggplot

df = data.frame(Fecha=anios, Precipitación=acumulado)

## Ploteo 

p = ggplot(df, aes(x=Fecha,y=Precipitación)) +
  geom_point(aes(group=1)) +
  geom_line(aes(group=1), size=1, color="navy")+
  geom_smooth(method = "lm")+
  ylab("mm")+xlab("Año")+
  ggtitle("Precipitación acumulada media")+
  theme_bw()+
  theme(plot.title = element_text(hjust = 0.5),
        panel.background=element_blank())

p

lm(Precipitación ~ Fecha, df)
