library(tidyverse)
library(corrplot)
library(sf)
library(spdep)
library(tmap)
library(RColorBrewer)
library(cluster)
library(dendextend)
library(spatialreg)
library(stargazer)


### FUNCIONS CREADES ###
## Creació mapes
crearMapa <- function(df, variable, bar=TRUE){
  dadesMapa <- inner_join(mapa, df, by=c("CODIMUNI" = "codiGenINE")) 
  if(!bar){
    dadesMapa[dadesMapa$nom=="Barcelona",variable] <- NA
  }
  grafic <-
    ggplot(dadesMapa) +
    geom_sf(aes(fill = .data[[variable]]), color = "grey70", size = 0.1) +
    {
      # Els colors usats s'han anat canviant al llarg del treball,
      # en funció del tipus de mapa creat
      if (is.numeric(dadesMapa[[variable]])) {
        scale_fill_gradient(low = "white", high = "steelblue", na.value = "grey70")
      } else {
        scale_fill_brewer(palette = "Set3", na.value = "grey70")
      }
    } +
    theme_minimal() +
    labs(fill = NULL,
         title = paste0("mapa ", variable),
         subtitle = "Catalunya") +
    theme(axis.text = element_blank(),
          axis.ticks = element_blank(),
          axis.title = element_blank(),
          panel.grid = element_blank(),
          panel.background = element_blank(),
          axis.line = element_blank())
  return(grafic)
}


# Funció per a fer mapes amb tm_shape
mapaEstEsp <- function(df,var,k=5){
  tm_shape(st_as_sf(df)) +
    tm_polygons(
      fill = var, 
      fill.scale = tm_scale_intervals(
        style = "quantile",
        n=k,
        #breaks = c(-0.3, -0.1, 0, 0.1, 0.3, 0.5), # si es vol fer breaks manual cal treure style
        values = "PuOr"  # Per a local Moran
        #values = "blues"  # Per a local G
      ),
      title.fill = "Local Moran's I") +
    tm_title(paste0(var)) +
    tm_layout(legend.outside = TRUE)
}


# Dataframes obtinguts a partir de les dades de l'INE, IDESCAT i Open Data BCN 
# i les corresponents transformacions, tal i com es descriu al treball. 
df # Dataframe amb les dades originals
dfLog # Dataframe amb les transformacions logarítmiques, si escau
dfLogEsc # Dataframe amb les transformacions logarítmiques, lògit i normalitzacions
mapa  # Dataframe de tipus shapefile amb els límits territorials dels 
      # districtes de Barcelona i els municipis de Catalunya



### ANÀLISI UNIVARIANT ###
# Histograma IST
ggplot(df, aes(x=IST)) + 
  geom_histogram(color="black",fill="#0072B2")

# Boxplot variables
boxplot <- list()
it <- 1
for(i in names(dfLog)[9:20]){
    boxplot[[it]] <- ggplot(dfLog, aes(x = .data[[i]])) + 
    geom_boxplot(color = "black", fill = "#0072B2")+
    coord_flip()
    it <- it + 1
}

# Mapes variables
mapes <- list()
it <- 1
for(i in names(dfLog)[9:20]){
  mapes[[it]] <- crearMapa(dfLog,i)
  it <- it + 1
}

# Estadístics bàsics
# Mínim, 1r quartil, mediana, 3r quartil, màxim
apply(dfLog %>% select(
  edatMitjanaHabitatges,
  superficieMitjana,
  percentatgeHabitatgesNoPrincipals,
  placesTuristiques,
  placesPerHabTuristic,
  habitatgesTurPerHab,
  poblacio,
  placesHotels,
  placesCampings,
  placesTurismeRural,
  placesHotLuxe,
), 2, quantile, probs = c(0, 0.25, 0.5, 0.75, 1), na.rm = TRUE) %>% round(2)

# Mitjana
apply(dfLog %>% select(
  edatMitjanaHabitatges,
  superficieMitjana,
  percentatgeHabitatgesNoPrincipals,
  placesTuristiques,
  placesPerHabTuristic,
  habitatgesTurPerHab,
  poblacio,
  placesHotels,
  placesCampings,
  placesTurismeRural,
  placesHotLuxe,
),2,mean) %>% round(2)



### ANÀLISI BIVARIANT ###
# Matriu de correlacions
dfLogEsc %>% select(
  IST,
  edatMitjanaHabitatges,
  superficieMitjana,
  percentatgeHabitatgesNoPrincipals,
  placesTuristiques,
  placesPerHabTuristic,
  habitatgesTurPerHab,
  poblacio,
  placesHotels,
  placesCampings,
  placesTurismeRural,
  placesHotLuxe,
) %>% 
  cor() %>% corrplot(
    tl.col = "black",
    tl.cex = 0.7
  )



### MATRIU DE PESOS ###
mapaUnics <- mapa %>%
  group_by(CODIMUNI) %>%
  summarise(geometry = st_union(geometry)) %>%
  ungroup()

# Pesos basats en distància entre els centroides
coords <- st_coordinates(st_centroid(mapaUnics))
row.names(coords) <- mapaUnics$CODIMUNI

dist.mat<-as.matrix(dist(coords, method="euclidean"))

# Invertim la distància només si és <= dmax, altrament es deixa a 0
dmax <- 30000
dist.mat.inv <- ifelse(dist.mat <= dmax & dist.mat > 0, 1 / dist.mat, 0)

diag(dist.mat.inv)<-0

pesos<-mat2listw(dist.mat.inv, style="W")

# Ordenem els df de manera que coincideixi amb la matriu de pesos
df <- df %>% arrange(codiGenINE)
dfLog <- dfLog %>% arrange(codiGenINE)
dfLogEsc <- dfLogEsc %>% arrange(codiGenINE)
dfPercent <- dfPercent %>% arrange(codiGenINE)



### ESTADÍSTICS ESPACIALS ###
## Contrastos generals d'autocorrelació espacial
# I de Moran
idemoran <- matrix(nrow=2,ncol=ncol(df)-8)
colnames(idemoran) <- names(df)[9:ncol(df)]
rownames(idemoran) <- c("Estadístic","pvalor")
it <- 1
for(j in 9:ncol(df)){
  x <- moran.test(dfLogEsc[,j], pesos, alternative="two.sided")
  idemoran[1,it] <- x$statistic
  idemoran[2,it] <- x$p.value
  it <- it+1
}

# G de Getis i Ord
# Creem pesos binaris per a aquest estadístic (funciona millor)
nbBin <- dnearneigh(coords, d1 = 0, d2 = 30000)  # 30.000 metres
pesosBin <- nb2listw(nbBin, style = "B")      # B = binary weights

gdegetis <- matrix(nrow=2,ncol=ncol(df)-8)
colnames(gdegetis) <- names(df)[9:ncol(df)]
rownames(gdegetis) <- c("Estadístic","pvalor")
it <- 1
for(j in 9:ncol(df)){
  x <- globalG.test(df[,j], pesosBin, alternative = "two.sided")
  gdegetis[1,it] <- x$statistic
  gdegetis[2,it] <- x$p.value
  it <- it+1
}


# Moran Plots
for(j in 9:ncol(dfLogEsc)){
  moran.plot(dfLogEsc[,j], pesos,
             xlab=names(dfLogEsc)[j],
             ylab=paste0("Spatially lagged ",names(dfLogEsc)[j]))
}


## Contrastos locals d'autocorrelació espacial
# Local Moran
localMoranZEst <- matrix(nrow=nrow(df),ncol=ncol(df)-8)
colnames(localMoranZEst) <- names(df)[9:ncol(df)]
rownames(localMoranZEst) <- mapaUnics$CODIMUNI

it <- 1
for(j in 9:ncol(df)){
  x <- as.data.frame(localmoran(dfLogEsc[,j], pesos, alternative="two.sided"))
  localMoranZEst[,it] <- x$Z.Ii
  it <- it+1
}

# Es passa a df per a poder ajuntar mapa unic i poder fer mapes
localMoranZEst <- as.data.frame(localMoranZEst) %>%
  mutate(codiGenINE=mapaUnics$CODIMUNI) %>%
  inner_join(mapaUnics, by=c("codiGenINE"="CODIMUNI"))

# Mapes
mapesLocalMoran <- list()
it <- 1
for(i in names(dfLog)[9:20][1:3]){
  mapesLocalMoran[[it]] <- mapaEstEsp(localMoranZEst,i)
  it <- it + 1
}

# New G
localGZEst <- matrix(nrow=nrow(df),ncol=ncol(df)-8)
colnames(localGZEst) <- names(df)[9:ncol(df)]
rownames(localGZEst) <- mapaUnics$CODIMUNI

it <- 1
for(j in 9:ncol(df)){
  x <- localG(dfLogEsc[,j], pesos, alternative="two.sided")
  x <- as.data.frame(attr(x,"internals"))
  localGZEst[,it] <- x$`Z(Gi)`
  it <- it+1
}

# Es passa a df per a poder ajuntar mapa unic i poder fer mapes
localGZEst <- as.data.frame(localGZEst) %>%
  mutate(codiGenINE=mapaUnics$CODIMUNI) %>%
  inner_join(mapaUnics, by=c("codiGenINE"="CODIMUNI"))

# Mapes
mapesLocalG <- list()
it <- 1
for(i in names(dfLog)[9:20]){
  mapesLocalG[[it]] <- mapaEstEsp(localGZEst,i)
  it <- it + 1
}



### CLUSTERING ###
data_scaled <- localGZEst %>% select(
  edatMitjanaHabitatges,
  superficieMitjana,
  percentatgeHabitatgesNoPrincipals,
  placesTuristiques,
  placesPerHabTuristic,
  habitatgesTurPerHab,
  poblacio,
  placesHotels,
  placesCampings,
  placesTurismeRural,
  placesHotLuxe,
)  # Al ser valors Z ja està escalat

dist <- daisy(data_scaled, metric = "euclidean")^2
hc <- hclust(dist, method = "ward.D2")


plot(hc,hang=-1,cex=0.6,labels=FALSE)

# Mètode de la silhouette
sil <- sapply(2:15, function(k) {
  grups <- cutree(hc, k = k)
  mean(silhouette(grups, dist)[, 3])
})

# Visualitzar resultats
plot(2:15, sil, type = "b", pch = 19,
     xlab = "Nombre de clústers (k)",
     ylab = "Mitjana silhouette")

# Es fa el dendograma amb 11 clústers
dend <- as.dendrogram(hc)
colors <- brewer.pal(11, "Set3")
# Colors ordenats per a que mapa i dendograma coincideixin
colors_mapa <- colors[c(1,11,2,6,5,4,8,10,7,3,9)] 
colors_dend <- color_branches(dend, k = 11, col = colors_mapa)
plot(colors_dend)

# Dendograma
dfLogEsc$jerarquic <- as.factor(cutree(hc,k=11))
crearMapa(dfLogEsc,"jerarquic")

# Es guarden els resultats a la resta de df
dfLog$jerarquic <- dfLogEsc$jerarquic
df$jerarquic <- dfLogEsc$jerarquic



### PROFILING ###
# Centroides
mitjanaClusters <- df %>%
  select(
    jerarquic,
    IST,
    edatMitjanaHabitatges,
    superficieMitjana,
    percentatgeHabitatgesNoPrincipals,
    placesTuristiques,
    placesPerHabTuristic,
    habitatgesTurPerHab,
    poblacio,
    placesHotels,
    placesCampings,
    placesTurismeRural,
    placesHotLuxe,
  ) %>%
  group_by(jerarquic) %>%
  summarise(across(where(is.numeric), mean, na.rm = TRUE))

mitjanes_long <- mitjanaClusters %>%
  pivot_longer(
    cols = -jerarquic,
    names_to = "variable",
    values_to = "valor"
  )

mitjana_global <- mitjanes_long %>%
  group_by(variable) %>%
  summarise(mitjana = mean(valor))

# Barplot amb línia de mitjana
ggplot() +
  geom_col(
    data = mitjanes_long,
    aes(x = jerarquic, y = valor, fill = jerarquic)
  ) +
  geom_hline(
    data = mitjana_global,
    aes(yintercept = mitjana),
    linetype = "dashed",
    color = "red",
    linewidth = 0.7
  ) +
  scale_fill_brewer(palette = "Set3") +
  scale_x_discrete(drop = FALSE) +
  facet_wrap(~ variable, scales = "free_y") +
  labs(
    title = "Mitjana per variable i jerarquia\nAmb línia de mitjana global",
    x = "Jerarquia",
    y = "Valor",
    fill = "Jerarquia"
  ) +
  theme_minimal(base_size = 13) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    strip.text = element_text(face = "bold"),
    panel.grid = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  )

# Mapa mitjana IST per clúster
mitjanes <- dfLogEsc %>%
  select(
    IST,
    edatMitjanaHabitatges,
    superficieMitjana,
    percentatgeHabitatgesNoPrincipals,
    placesTuristiques,
    placesPerHabTuristic,
    habitatgesTurPerHab,
    poblacio,
    placesHotels,
    placesCampings,
    placesTurismeRural,
    placesHotLuxe,
    jerarquic
  ) %>%
  group_by(jerarquic) %>%
  summarise(across(where(is.numeric),
                   ~ mean(.),
                   .names = "mitjana_{.col}")) %>%
  ungroup()

mapaCluster <- mapa %>%
  inner_join(dfLogEsc %>% select(codiGenINE,jerarquic),
             by=c("CODIMUNI" = "codiGenINE")) %>%
  inner_join(mitjanes %>% select(mitjana_IST,jerarquic))

mapaCluster$IST <- round(mapaCluster$mitjana_IST,2)

clusters_sf <- mapaCluster %>%
  group_by(jerarquic) %>%
  summarise(geometry = st_union(geometry)) %>%
  st_make_valid()

centroides <- st_point_on_surface(clusters_sf)

centroides <- centroides %>%
  inner_join(mitjanes %>% select(jerarquic,mitjana_IST))

centroides$IST <- round(centroides$mitjana_IST,2)

ggplot(mapaCluster) +
  geom_sf(aes(fill = mitjana_IST), color = "grey70", size = 0.1) +
  scale_fill_gradient(low = "white", high = "steelblue", na.value = "grey70") +
  geom_sf_text(data = centroides, aes(label = jerarquic), size = 4, fontface="bold") +
  theme_minimal() +
  labs(fill = NULL,
       title = "Mapa IST",
       subtitle = "Catalunya") +
  theme(axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.title = element_blank(),
        panel.grid = element_blank(),
        panel.background = element_blank(),
        axis.line = element_blank())

# Correlacions
cor <- dfLogEsc %>%
  select(
    IST,
    edatMitjanaHabitatges,
    superficieMitjana,
    percentatgeHabitatgesNoPrincipals,
    placesTuristiques,
    placesPerHabTuristic,
    habitatgesTurPerHab,
    poblacio,
    placesHotels,
    placesCampings,
    placesTurismeRural,
    placesHotLuxe,
    jerarquic
  ) %>%
  group_by(jerarquic) %>%
  summarise(across(where(is.numeric),
                   ~ round(cor(., IST, use = "pairwise.complete.obs"),2),
                   .names = "cor_{.col}")) %>%
  ungroup()

mapaCluster <- mapa %>%
  inner_join(dfLogEsc %>% select(codiGenINE,jerarquic),
             by=c("CODIMUNI" = "codiGenINE")) %>%
  inner_join(cor %>% select(
    cor_IST,
    cor_edatMitjanaHabitatges,
    cor_superficieMitjana,
    cor_percentatgeHabitatgesNoPrincipals,
    cor_placesTuristiques,
    cor_placesPerHabTuristic,
    cor_habitatgesTurPerHab,
    cor_poblacio,
    cor_placesHotels,
    cor_placesCampings,
    cor_placesTurismeRural,
    cor_placesHotLuxe,
    jerarquic))

clusters_sf <- mapaCluster %>%
  group_by(jerarquic) %>%
  summarise(geometry = st_union(geometry)) %>%
  st_make_valid()

centroides <- st_point_on_surface(clusters_sf)

correlacions <- list()
it <- 1
for(i in names(mapaCluster)[5:15]){
  correlacions[[it]] <- ggplot(mapaCluster) +
    geom_sf(aes(fill = .data[[i]]), color = "grey70", size = 0.1) +
    scale_fill_gradient2(
      low = "firebrick",
      mid = "white",
      high = "steelblue",
      midpoint = 0,
      limits = c(-1, 1),       # Escala fixa de -1 a +1
      na.value = "grey70"
    ) +
    geom_sf_text(data = centroides, aes(label = jerarquic), size = 4, fontface="bold") +
    theme_minimal() +
    labs(fill = NULL,
         title = paste0("Mapa ",i),
         subtitle = "Catalunya") +
    theme(axis.text = element_blank(),
          axis.ticks = element_blank(),
          axis.title = element_blank(),
          panel.grid = element_blank(),
          panel.background = element_blank(),
          axis.line = element_blank())
  it <- it+1
}



### MODELS ###
# Model Habitatges Turístics
model1ols <- lm(IST ~ placesTuristiques + I(placesTuristiques^2),
                  data=dfLog)
lm.LMtests(model1ols, listw = pesos, test = "all")
model1 <- errorsarlm(
  IST ~ placesTuristiques + I(placesTuristiques^2),
  listw=pesos,
  data=dfLog,
  Durbin= TRUE
)
summary(model1, Nagelkerke=TRUE)
impacts(model1,listw=pesos)
as.numeric(-coef(model1)[3]/(2*coef(model1)[4])) # Màxim local
max <- -impacts(model1,listw=pesos)$impacts$total[1]/
  (2*impacts(model1,listw=pesos)$impacts$total[2])

# Model Hotels
model2ols <- lm(IST ~ placesHotels + I(placesHotels^2),
                  data=dfLog)
lm.LMtests(model2ols, listw = pesos, test = "all")
pesosHot <- mat2listw(listw2mat(pesos)[dfLog$placesHotels>min(dfLog$placesHotels),
                                       dfLog$placesHotels>min(dfLog$placesHotels)],
                      style = "W")
model2 <- errorsarlm(
  IST ~ placesHotels,
  listw=pesosHot,
  data=dfLog %>%
    filter(placesHotels>0)
)
summary(model2, Nagelkerke=TRUE)
moran.test(model2$residuals,pesosHot)

# Model Càmpings
model3ols <- lm(IST ~ placesCampings + I(placesCampings^2),
                  data=dfLog)
lm.LMtests(model3ols, listw = pesos, test = "all")
pesosCamp <- mat2listw(listw2mat(pesos)[dfLog$placesCampings>min(dfLog$placesCampings),
                                        dfLog$placesCampings>min(dfLog$placesCampings)],
                       style = "W")
# moran.test(dfCamp %>% filter(grup==1) %>% select(IST) %>% unlist(),pesosCamp) # Autocorrelació espacial
moran.test(dfLog %>% filter(placesCampings>0) %>% select(IST) %>% unlist(),pesosCamp) 
model3 <- errorsarlm(
  IST ~ placesCampings + I(placesCampings^2),
  listw=pesosCamp,
  data=dfLog %>%
    filter(placesCampings>0)
)
summary(model3, Nagelkerke=TRUE)
moran.test(model3$residuals,pesosCamp)
as.numeric(-coef(model3)[3]/(2*coef(model3)[4])) # Màxim local

# Model Turisme Rural
model4ols <- lm(IST ~ placesTurismeRural + I(placesTurismeRural^2),
                  data=dfLog)
lm.LMtests(model4ols, listw = pesos, test = "all")
pesosRur <- mat2listw(listw2mat(pesos)[dfLogEsc$placesTurismeRural>min(dfLogEsc$placesTurismeRural),
                                       dfLogEsc$placesTurismeRural>min(dfLogEsc$placesTurismeRural)],
                      style = "W")
model4 <- errorsarlm(
  IST ~ placesTurismeRural,
  listw=pesosRur,
  data=dfLog %>%
    filter(placesTurismeRural>0)
)
summary(model4, Nagelkerke=TRUE)
moran.test(model4$residuals,pesosRur)

# Presentació testos
fmt_pval <- function(p) {
  ifelse(p < 2.2e-16, "<2.2e-16",
         formatC(p, format = "e", digits = 2))
}

extract_LM <- function(lmtest_obj) {
  tests <- c("RSerr", "RSlag", "adjRSerr", "adjRSlag")
  
  sapply(tests, function(t) {
    c(
      stat  = lmtest_obj[[t]]$statistic,
      pval  = lmtest_obj[[t]]$p.value
    )
  })
}

LM1 <- extract_LM(lm.LMtests(model1ols, listw = pesos, test = "all"))
LM2 <- extract_LM(lm.LMtests(model2ols, listw = pesos, test = "all"))
LM3 <- extract_LM(lm.LMtests(model3ols, listw = pesos, test = "all"))
LM4 <- extract_LM(lm.LMtests(model4ols, listw = pesos, test = "all"))

make_cell <- function(stat, pval) {
  paste0(
    round(stat, 3),
    " (", fmt_pval(pval), ")"
  )
}

LM_table <- matrix(NA, nrow = 4, ncol = 4)
rownames(LM_table) <- c("LM error", "LM lag", 
                        "LM-EL", "LM-LE")
colnames(LM_table) <- c("Places Turístiques", "Places Hotels",
                        "Places Càmpings", "Places Turisme Rural")

for (i in 1:4) {
  LM_table[, i] <- c(
    make_cell(LM1[1,"RSerr"],  LM1["pval","RSerr"]),
    make_cell(LM1[1,"RSlag"],  LM1["pval","RSlag"]),
    make_cell(LM1[1,"adjRSerr"], LM1["pval","adjRSerr"]),
    make_cell(LM1[1,"adjRSlag"], LM1["pval","adjRSlag"])
  )
  
  LM_list <- list(LM1, LM2, LM3, LM4)[[i]]
  
  LM_table[, i] <- c(
    make_cell(LM_list[1,"RSerr"],  LM_list["pval","RSerr"]),
    make_cell(LM_list[1,"RSlag"],  LM_list["pval","RSlag"]),
    make_cell(LM_list[1,"adjRSerr"], LM_list["pval","adjRSerr"]),
    make_cell(LM_list[1,"adjRSlag"], LM_list["pval","adjRSlag"])
  )
}

stargazer(
  as.data.frame(LM_table),
  type = "latex",
  summary = FALSE,
  title = "Test basats en el Multiplicador de Lagrange per a dependència espacial",
  label = "tab:LMtests",
  rownames = TRUE,
  table.placement = "htbp",
  notes = "Es reporta l'estadístic LM amb els corresponents pvalors en parèntesis.",
  float.env = "table",
  header = FALSE
)

# Presentació models
# Places turístiques
lambda <- model1$lambda
lambda_se <- model1$lambda.se
r2 <- summary(model1,Nagelkerke=TRUE)$NK
stargazer(
  model1,
  type = "latex",
  title = "Spatial error model (SEM)",
  dep.var.labels = "IST",
  covariate.labels = c(
    "log(Places Turístiques)",
    "log(Places Turístiques)$^2$",
    "Lag espacial de log(Places Turístiques)",
    "Lag espacial de log(Places Turístiques)$^2$",
    "Constant"
  ),
  add.lines = list(
    c("Spatial error parameter ($\\lambda$)",
      paste0(
        round(lambda, 3),
        " (", round(lambda_se, 3), ")"
      )
    ),
    c("Nagelkerke $R^2$", round(r2, 4)),
    c("Z value Moran Test", paste0(round(as.numeric(moran.test(model1$residuals,pesos)$statistic),3),"***"))
  ),
  omit.stat = c("f", "ser"),
  digits = 2,
  notes = "Spatial error model estimated by maximum likelihood."
)

impactes <- as.data.frame(impacts(model1,pesos)$impacts)
rownames(impactes) <- c(
  "log(Places Turístiques)",
  "log(Places Turístiques)$^2$"
)
colnames(impactes) <- c("directe","indirecte","total")
stargazer(
  impactes,
  type = "latex",
  title = "Impactes espacials del model",
  label = "tab:impactes",
  summary = FALSE,
  digits = 3
)

# Pendent model per municipis
crearMapa(dfLog %>%
            mutate(pendent=impactes[1,3]+2*impactes[2,3]*placesTuristiques),
          "pendent") # Posar colors adequats

# Pendent model per clúster
pendents <- dfLog %>%
  mutate(pendent=impactes[1,3]+2*impactes[2,3]*placesTuristiques) %>%
  group_by(jerarquic) %>%
  summarise(pendent2=mean(pendent))

mapaCluster <- mapa %>%
  inner_join(dfLogEsc %>% select(codiGenINE,jerarquic),
             by=c("CODIMUNI" = "codiGenINE")) %>%
  inner_join(pendents)

clusters_sf <- mapaCluster %>%
  group_by(jerarquic) %>%
  summarise(geometry = st_union(geometry)) %>%
  st_make_valid()

centroides <- st_point_on_surface(clusters_sf)

ggplot(mapaCluster) +
  geom_sf(aes(fill = pendent2), color = "grey70", size = 0.1) +
  scale_fill_gradient2(
    #limits=c(-max_abs,max_abs)
    low = "#8B1A1A",
    mid = "white",
    high = "#1A9850",
    midpoint = 0
  ) +
  geom_sf_text(data = centroides, aes(label = jerarquic), size = 4, fontface="bold") +
  theme_minimal() +
  labs(fill = NULL,
       title = paste0("Mapa ",i),
       subtitle = "Catalunya") +
  theme(axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.title = element_blank(),
        panel.grid = element_blank(),
        panel.background = element_blank(),
        axis.line = element_blank())

# Places hoteleres
moran.test(dfLog %>% filter(placesHotels>0) %>% select(IST) %>% unlist(),
           pesosHot)
lambda <- model2$lambda
lambda_se <- model2$lambda.se
r2 <- summary(model2,Nagelkerke=TRUE)$NK
stargazer(
  model2,
  type = "latex",
  title = "Spatial error model (SEM)",
  dep.var.labels = "IST",
  covariate.labels = c(
    "log(Places Hoteles)",
    "Constant"
  ),
  add.lines = list(
    c("Spatial error parameter ($\\lambda$)",
      paste0(
        round(lambda, 3),
        " (", round(lambda_se, 3), ")"
      )
    ),
    c("Nagelkerke $R^2$", round(r2, 4)),
    c("Z value Moran Test", paste0(round(as.numeric(moran.test(model2$residuals,pesosHot)$statistic),3),"")) # Posar tantes *** com escaigui pel pvalor
  ),
  omit.stat = c("f", "ser"),
  digits = 2,
  notes = "Spatial error model estimated by maximum likelihood."
)
dfLog %>%
  group_by(placesHotels!=0) %>%
  summarise(mean(IST),sd(IST))

# Places càmpings
moran.test(dfLog %>% filter(placesCampings>0) %>% select(IST) %>% unlist(),
           pesosCamp)
lambda <- model3$lambda
lambda_se <- model3$lambda.se
r2 <- summary(model3,Nagelkerke=TRUE)$NK
stargazer(
  model3,
  type = "latex",
  title = "Spatial error model (SEM)",
  dep.var.labels = "IST",
  covariate.labels = c(
    "log(Places Càmpings)",
    "log(Places Càmpings)$^2$",
    "Constant"
  ),
  add.lines = list(
    c("Spatial error parameter ($\\lambda$)",
      paste0(
        round(lambda, 3),
        " (", round(lambda_se, 3), ")"
      )
    ),
    c("Nagelkerke $R^2$", round(r2, 4)),
    c("Z value Moran Test", paste0(round(as.numeric(moran.test(model3$residuals,pesosCamp)$statistic),3),"")) # Posar tantes *** com escaigui pel pvalor
  ),
  omit.stat = c("f", "ser"),
  digits = 2,
  notes = "Spatial error model estimated by maximum likelihood."
)
-coef(model3)[3]/(2*coef(model3)[4])
crearMapa(dfLog %>%
            mutate(pendent=ifelse(placesCampings==0,NA,
                                  coef(model3)[3]+2*coef(model3)[4]*placesCampings)),
          "pendent") # Posar colors adequats
dfLog %>%
  group_by(placesCampings!=0) %>%
  summarise(mean(IST),sd(IST))

# Places turisme rural
moran.test(dfLog %>% filter(placesTurismeRural>0) %>% select(IST) %>% unlist(),
           pesosRur)
lambda <- model4$lambda
lambda_se <- model4$lambda.se
r2 <- summary(model4,Nagelkerke=TRUE)$NK
stargazer(
  model4,
  type = "latex",
  title = "Spatial error model (SEM)",
  dep.var.labels = "IST",
  covariate.labels = c(
    "log(Places Turisme Rural)",
    "Constant"
  ),
  add.lines = list(
    c("Spatial error parameter ($\\lambda$)",
      paste0(
        round(lambda, 3),
        " (", round(lambda_se, 3), ")"
      )
    ),
    c("Nagelkerke $R^2$", round(r2, 4)),
    c("Z value Moran Test", paste0(round(as.numeric(moran.test(model4$residuals,pesosRur)$statistic),3),"*")) # Posar tantes *** com escaigui pel pvalor
  ),
  omit.stat = c("f", "ser"),
  digits = 2,
  notes = "Spatial error model estimated by maximum likelihood."
)
dfLog %>%
  group_by(placesTurismeRural!=0) %>%
  summarise(mean(IST),sd(IST))
