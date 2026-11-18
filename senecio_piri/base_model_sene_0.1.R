#######################################################################################
#  Modelo de distribución del senecio inaequidens en el Pirineo Catalán
#  Autómatas celulares en tiempo discreto
#######################################################################################
# Versión corregida. Los cambios respecto al original están marcados con "# CORREGIDO".

# IMPORTANT: WORK ALWAYS IN THE SAME PROJECT "senecio_piri"

# 0. Instalación de paquetes  ---------------------------------------------

pkgs <- c("terra","sf","mapview", "ggplot2", "stars","usdm",
          "rgbif", "geodata", "maxnet", "readr", "lattice")   # packages needed

to_install <- setdiff(pkgs, rownames(installed.packages()))      # evaluate which packages are missing
if (length(to_install)) install.packages(to_install)             # install the missing packages
invisible(lapply(pkgs, library, character.only = TRUE))          # load pkgs in the script


set.seed(123)                                                    # fix the random process
dir.create("data", showWarnings = FALSE)                         # create a data folder
dir.create("out",  showWarnings = FALSE)                         # create an output folder

# 1. Ámbito de estudio  ---------------------------------------------------

# Rectangle coordinates (in degrees, WGS84)
xmin <- 1.30       # western long
xmax <- 2.00        # eastern long
ymin <- 42.03      # southern latit
ymax <- 42.30       # northern latit

# Create the bounding box with its coordinate reference system (CRS)
bbox <- st_bbox(c(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
                crs = 4326)

# Convert it into a polygon (sf object)
study_area <- st_as_sfc(bbox) |> st_sf(name = "study_area")

# visualize rectangle in a basemap
map_study_area <- mapview(study_area)
# save the study area
st_write(study_area, "data/study_area.gpkg", delete_dsn = TRUE, quiet = TRUE)   # GeoPackage

map_study_area
# 1.0 Raster --------------------------------------------------------------
size_cel <- 500                              # cell size in metres (ÚNICO valor para todo el script)
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

map_rast

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


# 1.2 / 1.3 Carreteras y usos del suelo desde cobertes_sol_2024.tif (ICGC) ------------------
# Un único archivo da las dos capas del modelo:
#   D : distancia (km) a la infraestructura más cercana (carreteras, ferrocarril y, opcional, núcleos urbanos)
#   S : idoneidad del uso del suelo (0-1), con las 41 clases del ICGC agrupadas en 10
# Las matrices se guardan en data/S_d_<tamaño>m.rds: el proceso a 1 m es lento y solo hay que hacerlo una vez.

# Agrupación de clases (se usa el número entre paréntesis de la leyenda del ICGC, 1-41)
grupos <- list(
  cultivos         = c(1, 2, 3, 4, 5, 6),       # herbáceos, huerta, viñas, olivos, otros leñosos, en transformación
  prados           = 14,                         # prados y herbazales
  matorral         = 10,
  bosque_denso     = c(7, 8, 9, 15),             # densos (aciculifolios, caducifolios, esclerófilos) y de ribera
  bosque_claro     = c(11, 12, 13),
  suelo_perturbado = c(16, 17, 29, 30, 32),      # suelo desnudo forestal, quemado, extracción, en transformación, desnudo urbano
  roca_nieve       = 18,                         # roquedos y neveros
  urbano           = 21:28,                      # casco urbano ... zonas deportivas
  agua_humedal     = c(20, 36, 37, 38, 39, 40),  # humedales, embalses, lagos, cursos, balsas, canales
  infraestructura  = c(31, 34)                   # carreteras y ferrocarril
)
# Clases eliminadas (S = 0; no aparecen en el Pirineo): 19 playas, 33 aeropuertos, 35 puertos, 41 mar, 0 sin datos

# Idoneidad S por grupo: VALORES PROVISIONALES, a calibrar 
S_grupo <- c(cultivos = 0.8, prados = 0.7, matorral = 0.5, bosque_denso = 0.15, bosque_claro = 0.4,
             suelo_perturbado = 0.9, roca_nieve = 0.1, urbano = 0.1, agua_humedal = 0, infraestructura = 0.5)

incluir_nucleos <- TRUE   # TRUE: d = distancia a carreteras/ferrocarril O núcleos urbanos (clases 21, 22, 23, 27)

# Código de 3 cifras de la leyenda (111-466) a partir del número 1-41
a_cod3 <- function(i) i + ifelse(i <= 6, 110, ifelse(i <= 20, 214, ifelse(i <= 35, 320, 425)))

f_cache <- paste0("data/S_d_", size_cel, "m.rds")

if (file.exists(f_cache)) {
  tmp <- readRDS(f_cache); S <- tmp$S; D <- tmp$D
} else {
  lc <- rast("data/cobertes_sol_2024.tif")
  if (!isTRUE(crs(lc, describe = TRUE)$code == "25831"))
    warning("El CRS del ráster no es EPSG:25831: comprueba que coincide con r_base")
  lc <- crop(lc, ext(r_base))                      # solo la zona de estudio (resolución original: 1 m)

  # tabla clase -> S, con los dos códigos posibles (1-41 o 111-466)
  tab <- data.frame(idx = unname(unlist(grupos)),
                    S   = rep(S_grupo[names(grupos)], lengths(grupos)))
  tab$cod3 <- a_cod3(tab$idx)

  # clase dominante a 25 m (aligera el cálculo) y detección del tipo de código del ráster
  lc25 <- aggregate(lc, fact = 25, fun = "modal")
  vals <- na.omit(unique(lc25)[, 1])
  usa_cod3 <- max(vals) > 41                       # TRUE si el ráster guarda 111-466 en vez de 1-41
  cod <- if (usa_cod3) tab$cod3 else tab$idx
  message("Clases presentes sin asignar (S = 0): ", paste(setdiff(vals, cod), collapse = ", "))

  # --- Distancia (km) a carreteras, ferrocarril y, opcionalmente, núcleos urbanos
  inf_idx <- c(grupos$infraestructura, if (incluir_nucleos) c(21, 22, 23, 27))
  inf_cod <- if (usa_cod3) a_cod3(inf_idx) else inf_idx
  inf25   <- aggregate(lc %in% inf_cod, fact = 25, fun = "max")  # "max": no se pierden vías estrechas
  inf_pts <- as.points(classify(inf25, cbind(0, NA)))            # solo los píxeles con infraestructura
  D <- as.matrix(distance(r_base, inf_pts), wide = TRUE) / 1000

  # --- Idoneidad S: reclasificar a 25 m y promediar en cada celda del modelo
  S25 <- classify(lc25, cbind(cod, tab$S), others = 0)           # eliminadas / sin datos -> 0
  S_r <- project(S25, r_base, method = "average")
  S   <- as.matrix(S_r, wide = TRUE)
  S[is.na(S)] <- 0

  saveRDS(list(S = S, D = D), f_cache)
}

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
stopifnot(identical(dim(S), dim(mat)), identical(dim(D), dim(mat)))


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
mod_sene_piri <- function(x, y, S, D, a, b, c, e) {
  if (!all(dim(x) == dim(y))) stop("x e y deben tener las mismas dimensiones")
  if (!all(dim(x) == dim(S))) stop("x y S deben tener las mismas dimensiones")
  if (!all(dim(x) == dim(D))) stop("x y d deben tener las mismas dimensiones")

  pc  <- S * (1 - (1 - a)^y * (1 - b * exp(-D / c)))   # probability of colonization
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

# 5. Management ------------------------------------------------------------

M0 <- matrix(0, nrow(mat), ncol(mat)) #Management matrix empty

# 5.1 by roads ------------------------------------------------------------

M_by_roads_matrix <- function(M0, D, d) {
  
  d <- matrix(d, nrow = nrow(M0), ncol = ncol(M0))
  
  M_roads <- matrix(as.integer(D < d),
                  nrow = nrow(M0),
                  ncol = ncol(M0))
  M_roads
}


# 5.2 by land cover -------------------------------------------------------


M_by_land_matrix <- function(M0, S, s) {        #value suitanibility to manage 
  
  s <- matrix(s, nrow = nrow(M0), ncol = ncol(M0))
  
  M_land <- matrix(as.integer(S > s),
                  nrow = nrow(M0),
                  ncol = ncol(M0))
 M_land
}

# 5.3 Random  -------------------------------------------------------------
M_rand <- matrix(
  sample(c(0, 1), nrow(mat)*ncol(mat),
         replace = TRUE,
         prob = c(0.5, 0.5)),
  nrow = nrow(mat),
  ncol = ncol(mat)
)

# 5.4 management funtion --------------------------------------------------
mod_sene_piri_mana <- function(x, y, S, D, M, a, b, c, e, n) {
  if (!all(dim(x) == dim(y))) stop("x e y deben tener las mismas dimensiones")
  if (!all(dim(x) == dim(S))) stop("x y S deben tener las mismas dimensiones")
  if (!all(dim(x) == dim(D))) stop("x y d deben tener las mismas dimensiones") 
  if (!all(dim(x) == dim(M))) stop("x y d deben tener las mismas dimensiones")
  
  pc  <- S * (1 - (1 - a)^y * (1 - b * exp(-D / c)))   # probability of colonization
  pc  <- pmin(pc, 1)
  ext <- e * (1 - S)
  resist <- (1-ext)*(1-n*M)
  resist <- pmax(pmin(resist, 1), 0)                         # probability of survive 
  
  u     <- matrix(runif(length(x)), nrow = nrow(x), ncol = ncol(x))
  vacia <- x == 0
  ocup  <- x == 1
  
  x_new <- x
  x_new[vacia] <- as.integer(u[vacia] <  pc[vacia])        # vacía -> ocupada con prob. pc
  x_new[ocup]  <- as.integer(u[ocup]  <  resist[ocup])     # ocupada -> ocupada con prob. resist
  x_new
} # end model



# 5. Simulaciones ---------------------------------------------------------

# Parámetros de ejemplo (valores provisionales, a calibrar)
par_mod <- list(a = 0.15,   # propagación local por vecino
                b = 0.05,   # propagación por carreteras
                c = 2,      # escala de decaimiento (km); d está en km
                e = 0.2)    # extinción máxima

# 5.0 Sin tratamiento -----------------------------------------------------

U <- mat
maxIter <- 10
ocupacion <- numeric(maxIter + 1)
ocupacion[1] <- mean(U[dentro])

for (j in 1:maxIter) {
  V <- neighbours(x = U, sp1 = 1)
  U <- mod_sene_piri(x = U, y = V, S = S, D = D,
                     a = par_mod$a, b = par_mod$b, c = par_mod$c, e = par_mod$e)
  ocupacion[j + 1] <- mean(U[dentro])
}

plot(0:maxIter, ocupacion, type = "b",
     xlab = "Paso de tiempo", ylab = "Proporción de celdas ocupadas")

r_U <- rast(U, crs = crs(pres_rast), extent = ext(pres_rast))
r_U <- mask(r_U, area_r)
writeRaster(r_U, "out/sene_final.tif", overwrite = TRUE)
mapview(r_U)


# 5.1  --------------------------------------------------------------------



simular <- function(M = NULL,        # matriz 0/1, o función(U, j) que devuelva una matriz 0/1
                    cada = 1,        # tratar cada 'cada' años (1, 2, 3, 4, 5...)
                    inicio = 1,      # primer año de tratamiento
                    maxIter = 10,
                    eta = 1,         # eficacia de eliminación
                    par = par_mod,
                    x0 = mat) {
  
  U <- x0
  ocupacion <- numeric(maxIter + 1)
  ocupacion[1] <- mean(U[dentro])
  
  for (j in 1:maxIter) {
    toca <- !is.null(M) && j >= inicio && ((j - inicio) %% cada == 0)
    Mj <- if (!toca) M0 else if (is.function(M)) M(U, j) else M
    
    V <- neighbours(x = U, sp1 = 1)
    U <- mod_sene_piri_mana(x = U, y = V, S = S, D = D, M = Mj,
                            a = par$a, b = par$b, c = par$c, e = par$e, n = eta)
    ocupacion[j + 1] <- mean(U[dentro])
  }
  list(U = U, ocupacion = ocupacion)
}

# 5.2 Tratamiento aleatorio -----------------------------------------------
nrep <- 30
frecuencias <- c(1, 2, 3, 4, 5)

res <- sapply(frecuencias, function(f) {
  rowMeans(replicate(nrep, simular(M = M_rand, cada = f)$ocupacion))
})
colnames(res) <- paste0("cada_", frecuencias)

matplot(0:10, res, type = "l", lty = 1, lwd = 2,
        xlab = "Paso de tiempo", ylab = "Proporción de celdas ocupadas")
legend("topright", colnames(res), col = 1:5, lty = 1, bty = "n")
# 5.3 tratamiento by roads ------------------------------------------------

# 5.4 Tratamiento por land cover ------------------------------------------


