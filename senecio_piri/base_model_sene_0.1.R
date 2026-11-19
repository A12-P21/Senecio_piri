#######################################################################################
#  Senecio inaequidens en el Pirineo Catalán
#  Autómata celular con gestión: tratamientos periódicos
#######################################################################################
#  Qué hace este script (en orden):
#    1. Parámetros            -> todo lo que se puede cambiar está aquí
#    2. Datos                 -> área, presencia inicial (mat), idoneidad (S) y distancia (D)
#    3. Modelo                -> funciones que pasan de un año al siguiente
#    4. Gestión               -> qué celdas se tratan según la estrategia
#    5. Simulación            -> repetir el modelo año a año, muchas veces
#    6. Escenarios            -> sin gestión + 3 estrategias x 5 frecuencias de tratamiento
#    7. Resultados            -> tablas y gráficos (carpetas out/ y out/fig/)
#
#  Trabajar siempre dentro del proyecto "senecio_piri".
#######################################################################################


# 0. Paquetes -------------------------------------------------------------

pkgs  <- c("terra", "sf", "ggplot2", "readr")
falta <- setdiff(pkgs, rownames(installed.packages()))   # paquetes que faltan
if (length(falta)) install.packages(falta)               # se instalan
invisible(lapply(pkgs, library, character.only = TRUE))  # se cargan

set.seed(123)                                            # para repetir los mismos resultados
dir.create("data",    showWarnings = FALSE)
dir.create("out",     showWarnings = FALSE)
dir.create("out/fig", showWarnings = FALSE)


# 1. Parámetros -----------------------------------------------------------

size_cel <- 500        # tamaño de celda (m)
maxIter  <- 10         # años simulados
nrep     <- 30         # repeticiones de cada escenario (el modelo es aleatorio)

# Parámetros del modelo (provisionales, a calibrar)
par_mod <- list(a = 0.15,   # propagación local (por cada vecino ocupado)
                b = 0.05,   # propagación por carreteras
                c = 2,      # cuánto decae la propagación con la distancia (km)
                e = 0.2)    # extinción máxima

# Gestión
eta         <- 1           # eficacia del tratamiento (1 = elimina todo; <1 = parcial)
inicio      <- 1           # primer año en que se trata
frecuencias <- 1:5         # se trata cada 1, 2, 3, 4 y 5 años
cada_ref    <- 1           # frecuencia que se usa en los mapas comparativos

# Umbrales: se trata TODA celda que cumple el criterio
d_trat  <- 1               # carreteras: celdas a menos de d_trat km
s_trat  <- 0.5             # uso del suelo: celdas con idoneidad S mayor que s_trat
p_aleat <- 0.5             # azar: probabilidad de tratar cada celda

# Nombres de las estrategias (para tablas y gráficos)
estrategias <- c(aleatoria  = "A. Aleatoria",
                 carreteras = "B. Carreteras",
                 suelo      = "C. Uso del suelo")

colores <- c("A. Aleatoria"     = "#7f7f7f",
             "B. Carreteras"    = "#1f77b4",
             "C. Uso del suelo" = "#2ca02c",
             "Sin gestión"      = "black")


# 2. Datos ----------------------------------------------------------------

# 2.1 Área de estudio y malla de celdas -------------------------------------
bbox <- st_bbox(c(xmin = 1.30, xmax = 2.00, ymin = 42.03, ymax = 42.30), crs = 4326)
study_area <- st_as_sfc(bbox) |> st_sf(name = "study_area") |> st_transform(25831)  # UTM 31N

r_base <- rast(vect(study_area), resolution = size_cel)       # malla de celdas
area_r <- rasterize(vect(study_area), r_base, field = 1)      # 1 dentro del área, NA fuera
dentro <- !is.na(as.matrix(area_r, wide = TRUE))              # matriz TRUE/FALSE: celda dentro del área

# 2.2 Presencia inicial: celdas con al menos un registro de GBIF ------------
occ   <- read_delim("data/occ_sene_piri.csv", show_col_types = FALSE)
coord <- occ[, c("decimalLatitude", "decimalLongitude")]
coord <- coord[complete.cases(coord), ]                       # quitar filas sin coordenadas
points <- st_as_sf(coord, coords = c("decimalLongitude", "decimalLatitude"), crs = 4326)

pres_vect <- vect(st_transform(points, st_crs(study_area)))   # puntos en UTM
pres_rast <- rasterize(pres_vect, r_base, field = 1, background = 0)  # 1 = hay registro, 0 = no
pres_rast <- mask(pres_rast, area_r)                          # NA fuera del área

mat <- as.matrix(pres_rast, wide = TRUE)                      # matriz de presencia x(0)
mat[is.na(mat)] <- 0

# 2.3 Idoneidad S (uso del suelo) y distancia D (km a carreteras) -----------
# Se calculan a partir de cobertes_sol_2024.tif (ICGC, 41 clases).
# El cálculo es lento, así que se hace una vez y se guarda en data/S_d_500m.rds.

f_cache <- paste0("data/S_d_", size_cel, "m.rds")

if (file.exists(f_cache)) {
  tmp <- readRDS(f_cache); S <- tmp$S; D <- tmp$D
} else {
  # Clases del ICGC agrupadas (número de la leyenda, 1-41)
  grupos <- list(
    cultivos         = c(1, 2, 3, 4, 5, 6),
    prados           = 14,
    matorral         = 10,
    bosque_denso     = c(7, 8, 9, 15),
    bosque_claro     = c(11, 12, 13),
    suelo_perturbado = c(16, 17, 29, 30, 32),
    roca_nieve       = 18,
    urbano           = 21:28,
    agua_humedal     = c(20, 36, 37, 38, 39, 40),
    infraestructura  = c(31, 34)          # carreteras (31) y ferrocarril (34): solo afecta a S
  )
  # Idoneidad de cada grupo (0-1): VALORES PROVISIONALES, a calibrar
  S_grupo <- c(cultivos = 0.8, prados = 0.7, matorral = 0.5, bosque_denso = 0.15,
               bosque_claro = 0.4, suelo_perturbado = 0.9, roca_nieve = 0.1,
               urbano = 0.1, agua_humedal = 0, infraestructura = 0.5)

  incluir_nucleos <- TRUE   # TRUE: D mide también la distancia a núcleos urbanos

  # Si el ráster usa códigos de 3 cifras (111-466) en vez de 1-41, se convierten
  a_cod3 <- function(i) i + ifelse(i <= 6, 110, ifelse(i <= 20, 214, ifelse(i <= 35, 320, 425)))

  lc <- rast("data/cobertes_sol_2024.tif")
  if (!isTRUE(crs(lc, describe = TRUE)$code == "25831"))
    warning("El CRS del ráster no es EPSG:25831: comprueba que coincide con r_base")
  lc <- crop(lc, ext(r_base))                       # solo la zona de estudio

  tab <- data.frame(idx = unname(unlist(grupos)),   # tabla clase -> idoneidad
                    S   = rep(S_grupo[names(grupos)], lengths(grupos)))
  tab$cod3 <- a_cod3(tab$idx)

  lc25     <- aggregate(lc, fact = 25, fun = "modal")   # clase dominante a 25 m (más rápido)
  vals     <- na.omit(unique(lc25)[, 1])
  usa_cod3 <- max(vals) > 41                            # ¿el ráster usa códigos de 3 cifras?
  cod      <- if (usa_cod3) tab$cod3 else tab$idx

  # D: distancia (km) a la carretera más cercana (clase 31; el ferrocarril, clase 34, no cuenta)
  inf_idx <- c(31, if (incluir_nucleos) c(21, 22, 23, 27))
  inf_cod <- if (usa_cod3) a_cod3(inf_idx) else inf_idx
  inf25   <- aggregate(lc %in% inf_cod, fact = 25, fun = "max")
  inf_pts <- as.points(classify(inf25, cbind(0, NA)))
  D <- as.matrix(distance(r_base, inf_pts), wide = TRUE) / 1000

  # S: cambiar cada clase por su idoneidad y promediar en cada celda del modelo
  S25 <- classify(lc25, cbind(cod, tab$S), others = 0)
  S   <- as.matrix(project(S25, r_base, method = "average"), wide = TRUE)
  S[is.na(S)] <- 0

  saveRDS(list(S = S, D = D), f_cache)
}

S[!dentro] <- 0     # fuera del área no se puede colonizar
mat[!dentro] <- 0
stopifnot(identical(dim(S), dim(mat)), identical(dim(D), dim(mat)))   # las 3 matrices miden lo mismo

n_dentro <- sum(dentro)
message("Celdas dentro del área: ", n_dentro, " | ocupadas al inicio: ", sum(mat))


# 3. Modelo ---------------------------------------------------------------
# Ecuaciones:
#   Celda vacía (colonización):   Pc  = S * [ 1 - (1 - a)^Y * (1 - b * exp(-D/c)) ]
#   Celda ocupada (extinción):    Pex = e * (1 - S)
#   Con gestión, la celda ocupada sobrevive con prob.:  (1 - Pex) * (1 - n * M)
# S = idoneidad, Y = nº de vecinos ocupados, D = distancia a carreteras (km)
# n = eficacia del tratamiento, M = 1 si la celda se trata (0 si no)

# Cuenta los vecinos ocupados (los 8 que rodean a cada celda)
neighbours <- function(x, sp1 = 1) {
  nr  <- nrow(x); nc <- ncol(x)
  pad <- matrix(0L, nr + 2, nc + 2)                       # borde de ceros alrededor
  pad[2:(nr + 1), 2:(nc + 1)] <- (x == sp1) * 1L
  out <- matrix(0L, nr, nc)
  for (di in -1:1) for (dj in -1:1) {
    if (di == 0 && dj == 0) next                          # la propia celda no cuenta
    out <- out + pad[(2:(nr + 1)) + di, (2:(nc + 1)) + dj]
  }
  out
}

# Un año SIN gestión
mod_sene_piri <- function(x, y, S, D, a, b, c, e) {
  if (!all(dim(x) == dim(y))) stop("x e y deben tener las mismas dimensiones")
  if (!all(dim(x) == dim(S))) stop("x y S deben tener las mismas dimensiones")
  if (!all(dim(x) == dim(D))) stop("x y D deben tener las mismas dimensiones")

  pc  <- S * (1 - (1 - a)^y * (1 - b * exp(-D / c)))   # prob. de colonización
  pc  <- pmin(pc, 1)
  ext <- pmin(e * (1 - S), 1)                           # prob. de extinción

  u     <- matrix(runif(length(x)), nrow = nrow(x), ncol = ncol(x))  # un número al azar por celda
  vacia <- x == 0
  ocup  <- x == 1

  x_new <- x
  x_new[vacia] <- as.integer(u[vacia] <  pc[vacia])     # vacía -> ocupada con prob. pc
  x_new[ocup]  <- as.integer(u[ocup]  >= ext[ocup])     # ocupada -> vacía con prob. ext
  x_new
}

# Un año CON gestión (M = matriz de celdas tratadas, n = eficacia)
mod_sene_piri_mana <- function(x, y, S, D, M, a, b, c, e, n) {
  if (!all(dim(x) == dim(y))) stop("x e y deben tener las mismas dimensiones")
  if (!all(dim(x) == dim(S))) stop("x y S deben tener las mismas dimensiones")
  if (!all(dim(x) == dim(D))) stop("x y D deben tener las mismas dimensiones")
  if (!all(dim(x) == dim(M))) stop("x y M deben tener las mismas dimensiones")

  pc  <- S * (1 - (1 - a)^y * (1 - b * exp(-D / c)))   # prob. de colonización
  pc  <- pmin(pc, 1)
  ext <- e * (1 - S)                                    # prob. de extinción natural
  resist <- (1 - ext) * (1 - n * M)                     # prob. de sobrevivir (extinción + tratamiento)
  resist <- pmax(pmin(resist, 1), 0)

  u     <- matrix(runif(length(x)), nrow = nrow(x), ncol = ncol(x))
  vacia <- x == 0
  ocup  <- x == 1

  x_new <- x
  x_new[vacia] <- as.integer(u[vacia] <  pc[vacia])     # vacía -> ocupada con prob. pc
  x_new[ocup]  <- as.integer(u[ocup]  <  resist[ocup])  # ocupada -> sigue ocupada con prob. resist
  x_new
}


# 4. Gestión --------------------------------------------------------------
# Cada estrategia devuelve una matriz M con 1 en las celdas que se tratan.

M0 <- matrix(0, nrow(mat), ncol(mat))                   # matriz vacía: no se trata nada

M_by_roads_matrix <- function(M0, D, d) {               # carreteras: celdas con D < d
  matrix(as.integer(D < d), nrow = nrow(M0), ncol = ncol(M0))
}

M_by_land_matrix <- function(M0, S, s) {                # uso del suelo: celdas con S > s
  matrix(as.integer(S > s), nrow = nrow(M0), ncol = ncol(M0))
}

M_random_matrix <- function(M0, p) {                    # azar: cada celda con prob. p
  matrix(as.integer(runif(length(M0)) < p), nrow = nrow(M0), ncol = ncol(M0))
}

# Elige la matriz según la estrategia (y solo cuenta celdas dentro del área)
matriz_trat <- function(estrategia) {
  M <- switch(estrategia,
              aleatoria  = M_random_matrix(M0, p_aleat),
              carreteras = M_by_roads_matrix(M0, D, d_trat),
              suelo      = M_by_land_matrix(M0, S, s_trat),
              stop("Estrategia desconocida: ", estrategia))
  M * dentro
}


# 5. Simulación -----------------------------------------------------------

# Una simulación de maxIter años.
#   estrategia: "ninguna", "aleatoria", "carreteras" o "suelo"
#   cada:       se trata cada 'cada' años, a partir del año 'inicio'
simular <- function(estrategia = "ninguna", cada = 1) {
  U    <- mat                                    # estado inicial
  ocup <- numeric(maxIter + 1)                   # proporción de celdas ocupadas cada año
  ocup[1] <- mean(U[dentro])
  trat <- matrix(0, nrow(U), ncol(U))            # nº de veces que se trata cada celda

  for (j in 1:maxIter) {
    V <- neighbours(x = U, sp1 = 1)              # vecinos ocupados
    if (estrategia == "ninguna") {
      U <- mod_sene_piri(x = U, y = V, S = S, D = D,
                         a = par_mod$a, b = par_mod$b, c = par_mod$c, e = par_mod$e)
    } else {
      toca <- j >= inicio && ((j - inicio) %% cada == 0)     # ¿se trata este año?
      Mj   <- if (toca) matriz_trat(estrategia) else M0
      U <- mod_sene_piri_mana(x = U, y = V, S = S, D = D, M = Mj,
                              a = par_mod$a, b = par_mod$b, c = par_mod$c,
                              e = par_mod$e, n = eta)
      trat <- trat + Mj
    }
    ocup[j + 1] <- mean(U[dentro])
  }
  list(U = U, ocup = ocup, trat = trat)
}

# Repite la simulación nrep veces (el modelo es aleatorio, cada vez sale algo distinto)
repetir <- function(estrategia = "ninguna", cada = 1) {
  ocup  <- matrix(0, maxIter + 1, nrep)          # una columna por repetición
  p_fin <- 0                                     # suma de los mapas finales
  trat  <- 0                                     # suma de los mapas de tratamiento
  for (r in 1:nrep) {
    s <- simular(estrategia, cada)
    ocup[, r] <- s$ocup
    p_fin <- p_fin + s$U
    trat  <- trat + s$trat
  }
  list(ocup  = ocup,
       p_fin = p_fin / nrep,    # probabilidad de presencia al final (0-1) en cada celda
       trat  = trat / nrep)     # nº medio de tratamientos por celda
}


# 6. Escenarios -----------------------------------------------------------

# 6.1 Sin gestión
sin_gestion <- repetir("ninguna")

# 6.2 Tres estrategias x cinco frecuencias (15 escenarios)
escenarios <- expand.grid(estrategia = names(estrategias), cada = frecuencias,
                          stringsAsFactors = FALSE)

res <- vector("list", nrow(escenarios))
for (i in seq_len(nrow(escenarios))) {
  message("Escenario ", i, " de ", nrow(escenarios), ": ",
          estrategias[[escenarios$estrategia[i]]], ", cada ", escenarios$cada[i], " año(s)")
  res[[i]] <- repetir(escenarios$estrategia[i], escenarios$cada[i])
}


# 7. Resultados -----------------------------------------------------------

# 7.1 Tablas ------------------------------------------------------------------

# Evolución de la ocupación: media y banda 5-95 % de las repeticiones
resumen <- function(ocup, estrategia, cada) {
  data.frame(estrategia = estrategia, cada = cada, paso = 0:maxIter,
             media = rowMeans(ocup),
             q05   = apply(ocup, 1, quantile, 0.05),
             q95   = apply(ocup, 1, quantile, 0.95))
}

tray_sg <- resumen(sin_gestion$ocup, "Sin gestión", NA)
tray <- do.call(rbind, lapply(seq_along(res), function(i)
  resumen(res[[i]]$ocup, estrategias[[escenarios$estrategia[i]]], escenarios$cada[i])))
tray$estrategia <- factor(tray$estrategia, levels = unname(estrategias))

# Resultado final de cada escenario y reducción respecto a no gestionar
final_sg <- mean(sin_gestion$ocup[maxIter + 1, ])
tabla <- escenarios
tabla$ocup_final    <- sapply(res, function(r) mean(r$ocup[maxIter + 1, ]))
tabla$reduccion_pct <- 100 * (1 - tabla$ocup_final / final_sg)
tabla$estrategia    <- factor(unname(estrategias[tabla$estrategia]), levels = unname(estrategias))

write.csv(tabla, "out/resumen_escenarios.csv", row.names = FALSE)
write.csv(rbind(tray_sg, transform(tray, estrategia = as.character(estrategia))),
          "out/trayectorias.csv", row.names = FALSE)

# 7.2 Gráficos ------------------------------------------------------------------

tema      <- theme_minimal(base_size = 11) + theme(legend.position = "bottom")
tema_mapa <- theme_void(base_size = 11)    + theme(legend.position = "bottom")

guardar <- function(p, nombre, w, h)
  ggsave(file.path("out/fig", nombre), p, width = w, height = h, dpi = 200, bg = "white")

lab_cada <- labeller(cada = function(x) paste0("cada ", x, " año(s)"))

# Convierte una matriz en tabla (x, y, valor) para dibujarla como mapa
mat_a_df <- function(m) {
  m <- m * 1
  m[!dentro] <- NA
  g <- expand.grid(fila = seq_len(nrow(m)), col = seq_len(ncol(m)))
  data.frame(x = xFromCol(r_base, g$col), y = yFromRow(r_base, g$fila), valor = as.vector(m))
}

# Escala de color común para los mapas de probabilidad (0-1)
escala_prob <- scale_fill_viridis_c(option = "magma", direction = -1, limits = c(0, 1),
                                    na.value = "grey95", name = "Probabilidad de presencia")

# --- Gráfico 1. Sin gestión: evolución de la ocupación ---------------------------
rep_sg <- data.frame(paso = rep(0:maxIter, nrep),
                     rep  = rep(1:nrep, each = maxIter + 1),
                     ocup = as.vector(sin_gestion$ocup))

g1 <- ggplot() +
  geom_line(data = rep_sg, aes(paso, 100 * ocup, group = rep), colour = "grey60", alpha = .5) +
  geom_ribbon(data = tray_sg, aes(paso, ymin = 100 * q05, ymax = 100 * q95), alpha = .15) +
  geom_line(data = tray_sg, aes(paso, 100 * media), linewidth = 1.2) +
  labs(title = "Sin gestión: evolución de la ocupación",
       subtitle = paste0(nrep, " repeticiones (gris), media (negro) y banda 5-95 %"),
       x = "Año", y = "Celdas ocupadas (%)") + tema
guardar(g1, "1_sin_gestion_evolucion.png", 7, 4.5)

# --- Gráfico 2. Sin gestión: mapa inicial y final ---------------------------------
d2 <- rbind(cbind(mat_a_df(mat),                panel = "Año 0 (presencia inicial)"),
            cbind(mat_a_df(sin_gestion$p_fin),  panel = paste0("Año ", maxIter, " (prob. de presencia)")))
g2 <- ggplot(d2, aes(x, y, fill = valor)) + geom_raster() + coord_equal() +
  escala_prob + facet_wrap(~panel) +
  labs(title = "Sin gestión: situación inicial y final") + tema_mapa
guardar(g2, "2_sin_gestion_mapas.png", 10, 4.5)

# --- Gráfico 3. Evolución con tratamiento (filas = estrategia, columnas = frecuencia) ---
sg_rep <- do.call(rbind, lapply(frecuencias, function(f) transform(tray_sg, cada = f)))
sg_rep$estrategia <- NULL        # así la línea "sin gestión" se repite en todas las filas

g3 <- ggplot() +
  geom_ribbon(data = sg_rep, aes(paso, ymin = 100 * q05, ymax = 100 * q95), fill = "grey70", alpha = .4) +
  geom_line(data = sg_rep, aes(paso, 100 * media), colour = "black", linetype = "dashed") +
  geom_ribbon(data = tray, aes(paso, ymin = 100 * q05, ymax = 100 * q95, fill = estrategia), alpha = .3) +
  geom_line(data = tray, aes(paso, 100 * media, colour = estrategia), linewidth = 1) +
  scale_colour_manual(values = colores, guide = "none") +
  scale_fill_manual(values = colores, guide = "none") +
  facet_grid(estrategia ~ cada, labeller = labeller(cada = function(x) paste0("cada ", x, " año(s)"))) +
  labs(title = "Evolución de la ocupación con tratamientos periódicos",
       subtitle = paste0("Línea discontinua: sin gestión. Eficacia = ", eta,
                         "; carreteras D < ", d_trat, " km; suelo S > ", s_trat,
                         "; azar p = ", p_aleat),
       x = "Año", y = "Celdas ocupadas (%)") + tema
guardar(g3, "3_evolucion_estrategias.png", 14, 7)

# --- Gráfico 4. Reducción de la ocupación respecto a no gestionar -----------------
g4 <- ggplot(tabla, aes(factor(cada), estrategia, fill = reduccion_pct)) +
  geom_tile(colour = "white") +
  geom_text(aes(label = sprintf("%.0f%%", reduccion_pct)), size = 4) +
  scale_fill_gradient2(low = "#b2182b", mid = "white", high = "#2166ac",
                       midpoint = 0, name = "Reducción (%)") +
  scale_y_discrete(limits = rev) +
  labs(title = paste0("Reducción de la ocupación en el año ", maxIter, " respecto a no gestionar"),
       x = "Frecuencia de tratamiento (cada n años)", y = NULL) + tema
guardar(g4, "4_reduccion.png", 7.5, 4)

# --- Gráfico 5. Mapas de presencia final: 15 escenarios ---------------------------
d5 <- do.call(rbind, lapply(seq_along(res), function(i)
  cbind(mat_a_df(res[[i]]$p_fin),
        estrategia = estrategias[[escenarios$estrategia[i]]],
        cada = escenarios$cada[i])))
d5$estrategia <- factor(d5$estrategia, levels = unname(estrategias))

g5 <- ggplot(d5, aes(x, y, fill = valor)) + geom_raster() + coord_equal() +
  escala_prob +
  facet_grid(estrategia ~ cada, labeller = labeller(cada = function(x) paste0("cada ", x, " año(s)"))) +
  labs(title = paste0("Probabilidad de presencia en el año ", maxIter)) + tema_mapa
guardar(g5, "5_mapas_presencia_final.png", 14, 7)

# --- Gráfico 6. Dónde se trata cada estrategia (frecuencia de referencia) ---------
i_ref <- which(escenarios$cada == cada_ref)
d6 <- do.call(rbind, lapply(i_ref, function(i)
  cbind(mat_a_df(res[[i]]$trat), panel = estrategias[[escenarios$estrategia[i]]])))
d6$panel <- factor(d6$panel, levels = unname(estrategias))

g6 <- ggplot(d6, aes(x, y, fill = valor)) + geom_raster() + coord_equal() +
  scale_fill_viridis_c(na.value = "grey95", name = "Nº medio de tratamientos") +
  facet_wrap(~panel, nrow = 1) +
  labs(title = paste0("Dónde se trata (cada ", cada_ref, " año(s))")) + tema_mapa
guardar(g6, "6_mapas_tratamiento.png", 14, 4.5)

message("Listo. Tablas en out/, gráficos en out/fig/")
