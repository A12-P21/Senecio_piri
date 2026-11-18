#######################################################################################
#  Modelo de distribución del senecio inaequidens en el Pirineo Catalán
#  Autómatas celulares en tiempo discreto
#  Bloques: 1. datos  2. idoneidad (MaxEnt)  3. AC estocástico  4. calibración  5. sensibilidad
#######################################################################################
# Versión corregida. Los cambios respecto al original están marcados con "# CORREGIDO".
# Las capas de carreteras y usos del suelo se extraen de data/cobertes_sol_2024.tif.

# IMPORTANT: WORK ALWAYS IN THE SAME PROJECT "senecio_piri"

# 0. Instalación de paquetes  ---------------------------------------------

pkgs <- c("terra","sf","mapview", "ggplot2", "stars","usdm",
          "rgbif", "geodata", "maxnet", "osmdata", "readr", "lattice")   # packages needed

to_install <- setdiff(pkgs, rownames(installed.packages()))      # evaluate which packages are missing
if (length(to_install)) install.packages(to_install)             # install the missing packages
invisible(lapply(pkgs, library, character.only = TRUE))          # load pkgs in the script


set.seed(123)                                                    # fix the random process
dir.create("data", showWarnings = FALSE)                         # create a data folder
dir.create("out",  showWarnings = FALSE)                         # create an output folder

# 1. Ámbito de estudio  ---------------------------------------------------

# Rectangle coordinates (in degrees, WGS84)
xmin <- 1.572456        # western long
xmax <- 2.154267        # eastern long
ymin <- 42.246714       # southern latit
ymax <- 42.450201       # northern latit

# Create the bounding box with its coordinate reference system (CRS)
bbox <- st_bbox(c(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
                crs = 4326)

# Convert it into a polygon (sf object)
study_area <- st_as_sfc(bbox) |> st_sf(name = "study_area")

# visualize rectangle in a basemap
map_study_area <- mapview(study_area)
# save the study area
st_write(study_area, "data/study_area.gpkg", delete_dsn = TRUE, quiet = TRUE)   # GeoPackage


# 1.0 Raster --------------------------------------------------------------
size_cel <- 1500                              # cell size in metres (ÚNICO valor para todo el script)
study_area <- st_transform(study_area, 25831) # transform area to UTM 31N

# CORREGIDO: plantilla raster de terra ÚNICA (antes había "r_500m", que no existía,
# y una segunda resolución de 1000 m más abajo, de modo que las matrices no coincidían)
r_base <- rast(vect(study_area), resolution = size_cel)

# Malla sf solo para visualizar (antes se llamaba "raster")
grid_sf <- st_make_grid(study_area, cellsize = size_cel, square = TRUE) |>
  st_sf()
grid_sf <- grid_sf[lengths(st_intersects(grid_sf, study_area)) > 0, ]

map_rast <- mapview(grid_sf, alpha.regions = 0, color = "red", lwd = 0.5) +
  mapview(study_area, alpha.regions = 0.2)
nrow(grid_sf)


# 1. Data  ---------------------------------------------------------------

# 1.1 Distribución de la especie ------------------------------------------
# download data of occurrence in GBIF (.csv): https://www.gbif.org/

occ_sene_piri <- read_delim("data/occ_sene_piri.csv")          # import data
coord_occ <- occ_sene_piri[, c("decimalLatitude",              # select coord columns
                               "decimalLongitude")]
coord_occ <- coord_occ[complete.cases(coord_occ), ]            # CORREGIDO: quitar filas sin coordenadas

# convertimos las coordenadas en puntos
points <- st_as_sf(coord_occ,
                   coords = c("decimalLongitude", "decimalLatitude"),
                   crs = 4326)
# generamos un mapa de puntos con la distribución
map_points <- mapview(points, color = "darkgreen", col.regions = "darkgreen",
                      cex = 1, alpha.regions = 0.8, cluster = TRUE)


# 1.2 Carreteras: distancia a la carretera más cercana (km) ----------------
# Se extraen directamente del raster cobertes_sol_2024.tif.
# Según la leyenda: 354 = Xarxa viària.
# La distancia se calcula a partir de los píxeles originales de carreteras,
# antes de adaptar los usos del suelo a la resolución del modelo. Así evitamos
# perder carreteras estrechas durante el remuestreo a 1500 m.

f_cobertes <- "data/cobertes_sol_2024.tif"
if (!file.exists(f_cobertes)) {
  stop("No se encuentra data/cobertes_sol_2024.tif. Guarda el raster en la carpeta data.")
}

lc_raw <- rast(f_cobertes)
if (nlyr(lc_raw) > 1) {
  warning("El TIFF tiene varias bandas; se utilizará la primera banda.")
  lc_raw <- lc_raw[[1]]
}
if (is.na(crs(lc_raw)) || crs(lc_raw) == "") {
  stop("El raster cobertes_sol_2024.tif no tiene CRS definido. Asígnale su sistema de coordenadas antes de continuar.")
}

# Recortar en el CRS original para reducir el volumen de datos.
study_area_src <- project(vect(study_area), crs(lc_raw))
lc_raw <- crop(lc_raw, study_area_src, snap = "out")
# 1.2 Carreteras: distancia a la carretera más cercana
# Extraemos las carreteras directamente del raster original.

# Proyectar el área de estudio al CRS del raster
study_area_src <- project(vect(study_area), crs(lc_raw))

# Recortar el raster al área de estudio
lc_raw <- crop(lc_raw, study_area_src, snap = "out")

# Extraer carreteras (código 351 = Xarxa viària)
road_src <- ifel(lc_raw == 351, 1, NA)

# Limitar las carreteras al área de estudio
road_src <- mask(road_src, study_area_src)

# Guardar raster de carreteras
writeRaster(
  road_src,
  "data/carreteras_original.tif",
  overwrite = TRUE
)

# Calcular distancia a la carretera más cercana
# No es necesario convertir los píxeles en puntos
d_road_src <- terra::distance(
  road_src,
  filename = "data/distancia_carreteras_original.tif",
  overwrite = TRUE
)

# Proyectar las distancias a la resolución del modelo (1500 m)
d_road <- project(
  d_road_src,
  r_base,
  method = "bilinear",
  filename = "data/distancia_carreteras_1500m.tif",
  overwrite = TRUE
)

# Convertir a matriz y expresar en kilómetros
d <- as.matrix(d_road, wide = TRUE) / 1000

# Distancia desde cada celda de la rejilla del modelo a los píxeles de carretera
# originales; terra devuelve metros, convertimos a kilómetros.
d_r <- distance(r_base, road_pts)
d <- as.matrix(d_r, wide = TRUE) / 1000


# 1.3 Usos del suelo: idoneidad S (0-1) por celda ---------------------------
# Se agrupan las 41 categorías de la leyenda en cinco grandes grupos ecológicos.
# Los valores S son hipótesis iniciales y deberán calibrarse/someterse a sensibilidad.

# Adaptar el raster categórico a la rejilla común del modelo.
lc <- project(lc_raw, r_base, method = "near")
lc <- mask(lc, area_r)
names(lc) <- "cobertes_sol"

# Tabla de reclasificación: código original -> idoneidad S.
# Agrícola y pastos: 111-116 y 228.
# Bosques densos: 221-223; bosques claros: 225-227; ribera: 229.
# Matorral y espacios abiertos: 224, 230-234.
# Zonas urbanas/transformadas: 341-353 y 355 (con valores diferenciados).
# La clase 354 (carreteras) se excluye de S porque ya se incorpora mediante d.
# Agua: 461-466.
rcl_S <- matrix(c(
  111, 0.8,
  112, 0.8,
  113, 0.8,
  114, 0.8,
  115, 0.8,
  116, 0.8,
  221, 0.2,
  222, 0.2,
  223, 0.2,
  224, 0.7,
  225, 0.4,
  226, 0.4,
  227, 0.4,
  228, 0.9,
  229, 0.3,
  230, 0.7,
  231, 0.8,
  232, 0.5,
  233, 0.3,
  234, 0.1,
  341, 0.3,
  342, 0.3,
  343, 0.5,
  344, 0.6,
  345, 0.5,
  346, 0.7,
  347, 0.5,
  348, 0.6,
  349, 0.8,
  350, 0.5,
  351, 0.7,
  352, 0.7,
  353, 0.6,
  354, 0.0,
  355, 0.3,
  461, 0.0,
  462, 0.0,
  463, 0.0,
  464, 0.0,
  465, 0.0,
  466, 0.0
), ncol = 2, byrow = TRUE)

S_r <- classify(lc, rcl = rcl_S, others = NA)
S_r <- mask(S_r, area_r)
names(S_r) <- "idoneidad_S"
writeRaster(S_r, "data/idoneidad_usos_suelo_1500m.tif", overwrite = TRUE)

# Raster categórico simplificado: 1 agrícola/pastos, 2 forestal,
# 3 espacios abiertos/naturales, 4 urbanizado/transformado, 5 agua.
# La clase 354 (carreteras) queda como NA porque se guarda por separado.
rcl_grupo <- matrix(c(
  111, 1, 112, 1, 113, 1, 114, 1, 115, 1, 116, 1, 228, 1,
  221, 2, 222, 2, 223, 2, 225, 2, 226, 2, 227, 2, 229, 2,
  224, 3, 230, 3, 231, 3, 232, 3, 233, 3, 234, 3,
  341, 4, 342, 4, 343, 4, 344, 4, 345, 4, 346, 4, 347, 4,
  348, 4, 349, 4, 350, 4, 351, 4, 352, 4, 353, 4, 355, 4,
  461, 5, 462, 5, 463, 5, 464, 5, 465, 5, 466, 5
), ncol = 2, byrow = TRUE)
lc_grupo <- classify(lc, rcl = rcl_grupo, others = NA)
lc_grupo <- mask(lc_grupo, area_r)
names(lc_grupo) <- "grupo_usos_suelo"
writeRaster(lc_grupo, "data/usos_suelo_agrupados_1500m.tif", overwrite = TRUE)

# Convertir a matriz para el autómata celular.
S <- as.matrix(S_r, wide = TRUE)
# Las celdas sin dato/clase no reconocida quedan con idoneidad 0.
S[is.na(S)] <- 0

# Guardar también el raster de usos del suelo original, adaptado a la rejilla.
writeRaster(lc, "data/cobertes_sol_1500m.tif", overwrite = TRUE)

# 4. Modelo AC  -----------------------------------------------------------

# 4.1 convertir el raster en una matriz de presencia (0/1) -----------------

# CORREGIDO: "puntos" no existía (el objeto se llama "points") y
# crs(raster) apuntaba a la malla sf; ahora se usa el CRS de study_area (UTM)
pres_vect <- vect(st_transform(points, st_crs(study_area)))

# Marcar con 1 las celdas que contienen al menos un punto
# - field = 1: valor que se asigna a las celdas con presencia
# - background = 0: valor del resto de celdas (en vez de NA)
pres_rast <- rasterize(pres_vect, r_base, field = 1, background = 0)   # CORREGIDO: r_base en vez de r_500m

# Dejar NA fuera del área de estudio
area_r    <- rasterize(vect(study_area), r_base, field = 1)            # 1 dentro, NA fuera
pres_rast <- mask(pres_rast, area_r)
names(pres_rast) <- "presencia"

mat <- as.matrix(pres_rast, wide = TRUE)     # filas x columnas

dim(mat)
mat[is.na(mat)] <- 0
table(mat, useNA = "ifany")                  # nº de celdas con 0 y 1
global(pres_rast, "sum", na.rm = TRUE)       # nº de celdas con presencia

# CORREGIDO: eliminado el segundo bloque con tam_celda <- 1000 que sobrescribía pres_rast

# Celdas dentro del área de estudio y comprobación de dimensiones
dentro <- !is.na(as.matrix(area_r, wide = TRUE))
S[!dentro] <- 0                              # fuera del área no se puede colonizar
stopifnot(identical(dim(S), dim(mat)), identical(dim(d), dim(mat)))


# 4.2 Función base modelo -------------------------------------------------

# CORREGIDO: "vecinos.R" no estaba disponible; esta función cuenta los vecinos de
# Moore (8) con presencia. Si tu vecinos.R ya define neighbours(), puedes borrar este bloque
neighbours <- function(x, sp1 = 1) {
  nr  <- nrow(x); nc <- ncol(x)
  pad <- matrix(0L, nr + 2, nc + 2)
  pad[2:(nr + 1), 2:(nc + 1)] <- (x == sp1) * 1L
  out <- matrix(0L, nr, nc)
  for (di in -1:1) for (dj in -1:1) {
    if (di == 0 && dj == 0) next
    out <- out + pad[(2:(nr + 1)) + di, (2:(nc + 1)) + dj]
  }
  out
}

# Ecuaciones:
#   Colonización (celda vacía):  Pc  = S · [ 1 − (1 − a)^Y · (1 − b · exp(−d/c)) ]
#   Extinción (celda ocupada):   Pex = e · (1 − S)
# S = matriz de idoneidad por uso del suelo
# Y = nº de vecinos de Moore ocupados
# a = prob. de propagación local (por vecino)
# b = prob. de propagación por carreteras
# d = distancia a carreteras (km);  c = escala de decaimiento (km)
# e = prob. máxima de extinción

# CORREGIDO: S, y, d se pasan como argumentos (antes eran variables globales inexistentes, y
# se usaba "Y" en vez de "y"); e^(...) pasa a exp(...) porque "e" es la extinción, no el
# número de Euler; los errores usan stop(); el doble bucle se vectoriza (mismo resultado, mucho más rápido)
mod_sene_piri <- function(x, y, S, d, a, b, c, e) {
  if (!all(dim(x) == dim(y))) stop("x e y deben tener las mismas dimensiones")
  if (!all(dim(x) == dim(S))) stop("x y S deben tener las mismas dimensiones")
  if (!all(dim(x) == dim(d))) stop("x y d deben tener las mismas dimensiones")

  pc  <- S * (1 - (1 - a)^y * (1 - b * exp(-d / c)))   # probability of colonization
  pc  <- pmin(pc, 1)
  ext <- pmin(e * (1 - S), 1)                           # probability of extinction

  u     <- matrix(runif(length(x)), nrow = nrow(x), ncol = ncol(x))
  vacia <- x == 0
  ocup  <- x == 1

  x_new <- x
  x_new[vacia] <- as.integer(u[vacia] <  pc[vacia])     # vacía -> ocupada con prob. pc
  x_new[ocup]  <- as.integer(u[ocup]  >= ext[ocup])     # ocupada -> vacía con prob. ext
  x_new
} # end model

# 5. Managment ------------------------------------------------------------

M <- matrix(data = NA, nrow = nrow(mat), ncol = ncol(mat))


# 5.1 by roads ------------------------------------------------------------


# 5.2 by land cover -------------------------------------------------------


# 5.3 Random  -------------------------------------------------------------



# 5. Simulaciones ---------------------------------------------------------


# Parámetros de ejemplo (valores provisionales, a calibrar)
par_mod <- list(a = 0.15,   # propagación local por vecino
                b = 0.05,   # propagación por carreteras
                c = 2,      # escala de decaimiento (km); d está en km
                e = 0.2)    # extinción máxima

# CORREGIDO: eliminada la llamada "U <- mod_sene_piri()" sin argumentos;
# en el bucle faltaban x, S, d, a, b y se pasaba "ext =" en lugar de "e ="
U <- mat
maxIter <- 10
ocupacion <- numeric(maxIter + 1)
ocupacion[1] <- mean(U[dentro])

for (j in 1:maxIter) {
  V <- neighbours(x = U, sp1 = 1)
  U <- mod_sene_piri(x = U, y = V, S = S, d = d,
                     a = par_mod$a, b = par_mod$b, c = par_mod$c, e = par_mod$e)
  ocupacion[j + 1] <- mean(U[dentro])
}

plot(0:maxIter, ocupacion, type = "b",
     xlab = "Paso de tiempo", ylab = "Proporción de celdas ocupadas")

r_U <- rast(U, crs = crs(pres_rast), extent = ext(pres_rast))
r_U <- mask(r_U, area_r)
writeRaster(r_U, "out/sene_final.tif", overwrite = TRUE)
mapview(r_U)

