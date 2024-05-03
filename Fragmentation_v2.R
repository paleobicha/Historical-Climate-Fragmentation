# -----------------------------------------------------------------------
# Project: Historical Climatic Fragmentation
# File name: Fragmentation_v2.R
# Last updated: 2023-12-12
# Author: Sara Gamboa
# Email: sara.gamboa@uvigo.es
# Repository: https://github.com/paleobicha/Historical-Climate-Fragmentation
# -----------------------------------------------------------------------

library(raster)
library(sp)
library(RColorBrewer)
library(tictoc)
library(dplyr)
library(rgdal)
library(landscapemetrics)
library(beepr)
library(ggplot2)
library(plyr)
library(tidyverse)
library(terra)
library(jmv)
######


setwd("~/Dropbox/Mac/Documents/Sara Gamboa/Papers/Fragmentation/Reclassified_data")

moll_1to500 <- stack("Stack_reclassified_moll_1to500.tif")
moll_501to1000 <- stack("Stack_reclassified_moll_501to1000_V4.tif")
moll_1001to1500 <- stack("Stack_reclassified_moll_1001to1500.tif")
moll_1501to2000 <- stack("Stack_reclassified_moll_1501to2000.tif")
moll_2001to2500 <- stack("Stack_reclassified_moll_2001to2500.tif")
moll_2501to3000 <- stack("Stack_reclassified_moll_2501to3000.tif")
moll_3001to3500 <- stack("Stack_reclassified_moll_3001to3500.tif")
moll_3501to4000 <- stack("Stack_reclassified_moll_3501to4000.tif")
moll_4001to4500 <- stack("Stack_reclassified_moll_4001to4500_V2.tif")
moll_4501to5001 <- stack("Stack_reclassified_moll_4501to5001.tif")

stack_all <- stack(moll_1to500, moll_501to1000, moll_1001to1500, moll_1501to2000,
                   moll_2001to2500, moll_2501to3000, moll_3001to3500, moll_3501to4000,
                   moll_4001to4500, moll_4501to5001)

rm(moll_1to500, moll_501to1000,moll_1001to1500, moll_1501to2000, moll_2001to2500, 
   moll_2501to3000, moll_3001to3500, moll_3501to4000, moll_4001to4500, moll_4501to5001)

stack_all #Stack con 5001 rasters.
plot(stack_all[[20]])
plot(Eurasia,add=T)
Eurasia@data

library(landscapemetrics)
check_landscape(stack_all[[1]]) #OK

t <- list_lsm()

#### Estimate the time before running the hole code
tic("tiempo")
patch_metrics <- dplyr::bind_rows(
  lsm_p_area(stack_all[[1:10]]),
  lsm_p_enn(stack_all[[1:10]]),
  lsm_p_para(stack_all[[1:10]]),
  lsm_p_shape(stack_all[[1:10]]),
  lsm_c_ai(stack_all[[1:10]]),
)
toc()

###Run all functions at the same time
patch_metrics <- dplyr::bind_rows(
  lsm_p_area(stack_all),
  lsm_p_enn(stack_all),
  lsm_p_para(stack_all),
  lsm_p_shape(stack_all),
  lsm_c_ai(stack_all),
)

# look at the results
patch_metrics_full_names <- dplyr::left_join(x = patch_metrics,
                                             y = lsm_abbreviations_names, 
                                             by = "metric")
patch_metrics_full_names

setwd("~/Dropbox/Mac/Documents/Sara Gamboa/Papers/Fragmentation/Results")
save(patch_metrics, file = "patch_metrics.RData")
patch_metrics_full_names <- patch_metrics_full_names[patch_metrics_full_names$class>0,]
save(patch_metrics_full_names, file = "patch_metrics_full_names.RData")

load("patch_metrics_full_names.RData")

####Masking by continent

setwd("~/Dropbox/Mac/Documents/Sara Gamboa/Papers/Fragmentation/Reclassified_data/Continents_masks")
list.files()
Africa <- readOGR("Africa.shp")
America <- readOGR("America3.shp")
Eurasia <- readOGR("Eurasia_mod.shp")

America_p <- raster(ncol=302, nrow=732)
extent(America_p) <- extent(stack_all)
values(America_p) <- 1
crs.mol <- CRS("+proj=moll +lon_0=0 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs")
crs(America_p) <- crs.mol
result <- America_p %>% mask(America)
plot(result)
yyy <- resample(result,stack_all,method='bilinear')
yyy@data@values[yyy@data@values>0] <- 1
plot(yyy)
America_stack <- stack_all*yyy

writeRaster(America_stack, filename = "America_stack.tif", format = "GTiff")

stack_0 <- stack_all
for(i in c(1:5001)){
  stack_0[[i]][is.na(stack_0[[i]])] <- 0
  print(i)
}

Eurasia_p <- raster(ncol=302, nrow=732)
extent(Eurasia_p) <- extent(stack_all)
values(Eurasia_p) <- 1
crs.mol <- CRS("+proj=moll +lon_0=0 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs")
crs(Eurasia_p) <- crs.mol
result <- Eurasia_p %>% mask(Eurasia)
plot(result)
yyy <- resample(result,stack_all,method='bilinear')
yyy@data@values[yyy@data@values>0] <- 1
plot(yyy)
Eurasia_stack <- stack_0*yyy

for(i in c(1:5001)){
  Eurasia_stack[[i]][Eurasia_stack[[1]]<1]<- NA
  print(i)
}

writeRaster(Eurasia_stack, filename = "Eurasia_stack.tif", format = "GTiff")

Africa_p <- raster(ncol=302, nrow=732)
extent(Africa_p) <- extent(stack_all)
values(Africa_p) <- 1
crs.mol <- CRS("+proj=moll +lon_0=0 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs")
crs(Africa_p) <- crs.mol
result <- Africa_p %>% mask(Africa)
plot(result)
yyy <- resample(result,stack_all,method='bilinear')
yyy@data@values[yyy@data@values>0] <- 1
plot(yyy)
Africa_stack <- stack_0*yyy
plot(Africa_stack[[1]])

for(i in c(1:5001)){
  Africa_stack[[i]][Africa_stack[[1]]<1]<- NA
  print(i)
}

writeRaster(Africa_stack, filename = "Africa_stack.tif", format = "GTiff", overwrite=T)
setwd("~/Dropbox/Mac/Documents/Sara Gamboa/Papers/Fragmentation/Reclassified_data/Continents_masks")
list.files()
Africa_stack <- stack("Africa_stack.tif")
Eurasia_stack <- stack("Eurasia_stack.tif")
America_stack <- stack("America_stack.tif")
plot(Africa_stack[[2]])
quartz()
########Metrics for region
setwd("~/Dropbox/Mac/Documents/Sara Gamboa/Papers/Fragmentation/Metrics")

tic("tiempo")
America_metrics <- dplyr::bind_rows(
  lsm_p_area(America_stack[[1:10]]),
  lsm_p_enn(America_stack[[1:10]]),
  lsm_p_para(America_stack[[1:10]]),
  lsm_p_shape(America_stack[[1:10]]),
  lsm_c_ai(America_stack[[1:10]]),
)
toc()

###Run all functions at the sime time
America_metrics <- dplyr::bind_rows(
  lsm_p_area(America_stack),
  lsm_p_enn(America_stack),
  lsm_p_para(America_stack),
  lsm_p_shape(America_stack),
  lsm_c_ai(America_stack),
)

America_metrics_full_names <- dplyr::left_join(x = America_metrics,
                                               y = lsm_abbreviations_names, 
                                               by = "metric")
America_metrics_full_names

save(America_metrics, file = "America_metrics.RData")
save(America_metrics_full_names, file = "America_metrics_full_names.RData")

Africa_metrics <- dplyr::bind_rows(
  lsm_p_area(Africa_stack),
  lsm_p_enn(Africa_stack),
  lsm_p_para(Africa_stack),
  lsm_p_shape(Africa_stack),
  lsm_c_ai(Africa_stack),
)
beep()

Africa_metrics_full_names <- dplyr::left_join(x = Africa_metrics,
                                              y = lsm_abbreviations_names, 
                                              by = "metric")
Africa_metrics_full_names <- Africa_metrics_full_names[Africa_metrics_full_names$class>0,]
Africa_metrics_full_names

save(Africa_metrics, file = "Africa_metrics.RData")
save(Africa_metrics_full_names, file = "Africa_metrics_full_names.RData")

Eurasia_metrics <- dplyr::bind_rows(
  lsm_p_area(Eurasia_stack),
  lsm_p_enn(wwrasia_stack),
  lsm_p_para(wwrasia_stack),
  lsm_p_shape(wwrasia_stack),
  lsm_c_ai(wwrasia_stack),
)
beep()

Eurasia_metrics_full_names <- dplyr::left_join(x = Eurasia_metrics,
                                               y = lsm_abbreviations_names, 
                                               by = "metric")
Eurasia_metrics_full_names

save(Eurasia_metrics, file = "Eurasia_metrics.RData")
Eurasia_metrics_full_names <- Eurasia_metrics_full_names[Eurasia_metrics_full_names$class>0,]
save(Eurasia_metrics_full_names, file = "Eurasia_metrics_full_names.RData")

######
world_metrics <- patch_metrics_full_names[,c(1,3,5,6)]

unique(world_metrics$metric)
unique(patch_metrics_full_names$name)
world_metrics <- patch_metrics_full_names[,c(1,3,5,6)]
area_metrics <- world_metrics[world_metrics$metric=="area",]
enn_metrics <- world_metrics[world_metrics$metric=="enn",]
para_metrics <- world_metrics[world_metrics$metric=="para",]
shape_metrics <- world_metrics[world_metrics$metric=="shape",]
aggre_metrics <- world_metrics[world_metrics$metric=="ai",]

##Africa
african_metrics <- Africa_metrics_full_names[,c(1,3,5,6)]
african_metrics$continent <- "Africa"
area_Africa <- african_metrics[african_metrics$metric=="area",]
enn_Africa <- african_metrics[african_metrics$metric=="enn",]
para_Africa <- african_metrics[african_metrics$metric=="para",]
shape_Africa <- african_metrics[african_metrics$metric=="shape",]
aggre_Africa <- african_metrics[african_metrics$metric=="ai",]

##America
american_metrics <- America_metrics_full_names[,c(1,3,5,6)]
american_metrics$continent <- "America"
area_America <- american_metrics[american_metrics$metric=="area",]
enn_America <- american_metrics[american_metrics$metric=="enn",]
para_America <- american_metrics[american_metrics$metric=="para",]
shape_America <- american_metrics[american_metrics$metric=="shape",]
aggre_America <- american_metrics[american_metrics$metric=="ai",]

##Eurasia
eurasian_metrics <- Eurasia_metrics_full_names[,c(1,3,5,6)]
eurasian_metrics$continent <- "eurasia"
area_Eurasia <- eurasian_metrics[eurasian_metrics$metric=="area",]
enn_eurasia <- eurasian_metrics[eurasian_metrics$metric=="enn",]
para_eurasia <- eurasian_metrics[eurasian_metrics$metric=="para",]
shape_eurasia <- eurasian_metrics[eurasian_metrics$metric=="shape",]
aggre_eurasia <- eurasian_metrics[eurasian_metrics$metric=="ai",]

full_metrics <- rbind(american_metrics, african_metrics, eurasian_metrics)
area_metrics <- full_metrics[full_metrics$metric=="area",]

area_ww <- data.frame(matrix(data=NA, ncol = 5, nrow = 2000000))


max_ln <- max(c(length(x), length(y)))
area_names <- c("Tropical", "Arid", "Temperate", "Cold", "Polar")

i <- 1
for (i in 1:5) {
  area_ww[,i] <- area_metrics[area_metrics$class==i,4]
  colnames(area_ww)[i] <- area_names[i]
}

gfg_data<- data.frame(col1 = c(x,rep(NA, max_ln - length(x))))

area_trop_ww <- area_metrics[area_metrics$class==1,]
area_arid_ww <- area_metrics[area_metrics$class==2,]
area_temp_ww <- area_metrics[area_metrics$class==3,]
area_cold_ww <- area_metrics[area_metrics$class==4,]
area_pola_ww <- area_metrics[area_metrics$class==5,]

area_trop_am <- area_trop_ww[area_trop_ww$continent=="America",]
area_arid_am <- area_arid_ww[area_arid_ww$continent=="America",]
area_temp_am <- area_temp_ww[area_temp_ww$continent=="America",]
area_cold_am <- area_cold_ww[area_cold_ww$continent=="America",]
area_pola_am <- area_pola_ww[area_pola_ww$continent=="America",]

area_trop_af <- area_trop_ww[area_trop_ww$continent=="Africa",]
area_arid_af <- area_arid_ww[area_arid_ww$continent=="Africa",]
area_temp_af <- area_temp_ww[area_temp_ww$continent=="Africa",]
area_cold_af <- area_cold_ww[area_cold_ww$continent=="Africa",]
area_pola_af <- area_pola_ww[area_pola_ww$continent=="Africa",]

area_trop_eu <- area_trop_ww[area_trop_ww$continent=="eurasia",]
area_arid_eu <- area_arid_ww[area_arid_ww$continent=="eurasia",]
area_temp_eu <- area_temp_ww[area_temp_ww$continent=="eurasia",]
area_cold_eu <- area_cold_ww[area_cold_ww$continent=="eurasia",]
area_pola_eu <- area_pola_ww[area_pola_ww$continent=="eurasia",]

unique(area_trop_ww$continent)

####Categories by patch size
##Worldwide
area_trop_ww_s <- area_trop_ww[area_trop_ww$value<310000,] ###because value is in hectares
area_trop_ww_m <- area_trop_ww[area_trop_ww$value>310000 & area_trop_ww$value<3000000,]
area_trop_ww_l <- area_trop_ww[area_trop_ww$value>3000000 & area_trop_ww$value<60000000,]
area_trop_ww_xl <- area_trop_ww[area_trop_ww$value>60000000,]

area_arid_ww_s <- area_arid_ww[area_arid_ww$value<310000,] ###because value is in hectares
area_arid_ww_m <- area_arid_ww[area_arid_ww$value>310000 & area_arid_ww$value<3000000,]
area_arid_ww_l <- area_arid_ww[area_arid_ww$value>3000000 & area_arid_ww$value<60000000,]
area_arid_ww_xl <- area_arid_ww[area_arid_ww$value>60000000,]

area_temp_ww_s <- area_temp_ww[area_temp_ww$value<310000,] ###because value is in hectares
area_temp_ww_m <- area_temp_ww[area_temp_ww$value>300000 & area_temp_ww$value<3000000,]
area_temp_ww_l <- area_temp_ww[area_temp_ww$value>3000000 & area_temp_ww$value<60000000,]
area_temp_ww_xl <- area_temp_ww[area_temp_ww$value>60000000,]

area_cold_ww_s <- area_cold_ww[area_cold_ww$value<310000,] ###because value is in hectares
area_cold_ww_m <- area_cold_ww[area_cold_ww$value>300000 & area_cold_ww$value<3000000,]
area_cold_ww_l <- area_cold_ww[area_cold_ww$value>3000000 & area_cold_ww$value<60000000,]
area_cold_ww_xl <- area_cold_ww[area_cold_ww$value>60000000,]

area_pola_ww_s <- area_pola_ww[area_pola_ww$value<310000,] ###because value is in hectares
area_pola_ww_m <- area_pola_ww[area_pola_ww$value>300000 & area_pola_ww$value<3000000,]
area_pola_ww_l <- area_pola_ww[area_pola_ww$value>3000000 & area_pola_ww$value<60000000,]
area_pola_ww_xl <- area_pola_ww[area_pola_ww$value>60000000,]

##America
area_trop_am_s <- area_trop_am[area_trop_am$value<310000,] ###because value is in hectares
area_trop_am_m <- area_trop_am[area_trop_am$value>300000 & area_trop_am$value<3000000,]
area_trop_am_l <- area_trop_am[area_trop_am$value>3000000 & area_trop_am$value<60000000,]
area_trop_am_xl <- area_trop_am[area_trop_am$value>60000000,]

area_arid_am_s <- area_arid_am[area_arid_am$value<310000,] ###because value is in hectares
area_arid_am_m <- area_arid_am[area_arid_am$value>300000 & area_arid_am$value<3000000,]
area_arid_am_l <- area_arid_am[area_arid_am$value>3000000 & area_arid_am$value<60000000,]
area_arid_am_xl <- area_arid_am[area_arid_am$value>60000000,]

area_temp_am_s <- area_temp_am[area_temp_am$value<310000,] ###because value is in hectares
area_temp_am_m <- area_temp_am[area_temp_am$value>300000 & area_temp_am$value<3000000,]
area_temp_am_l <- area_temp_am[area_temp_am$value>3000000 & area_temp_am$value<60000000,]
area_temp_am_xl <- area_temp_am[area_temp_am$value>60000000,]

area_cold_am_s <- area_cold_am[area_cold_am$value<310000,] ###because value is in hectares
area_cold_am_m <- area_cold_am[area_cold_am$value>300000 & area_cold_am$value<3000000,]
area_cold_am_l <- area_cold_am[area_cold_am$value>3000000 & area_cold_am$value<60000000,]
area_cold_am_xl <- area_cold_am[area_cold_am$value>60000000,]

area_pola_am_s <- area_pola_am[area_pola_am$value<310000,] ###because value is in hectares
area_pola_am_m <- area_pola_am[area_pola_am$value>300000 & area_pola_am$value<3000000,]
area_pola_am_l <- area_pola_am[area_pola_am$value>3000000 & area_pola_am$value<60000000,]
area_pola_am_xl <- area_pola_am[area_pola_am$value>60000000,]

##Africa
area_trop_af_s <- area_trop_af[area_trop_af$value<310000,] ###because value is in hectares
area_trop_af_m <- area_trop_af[area_trop_af$value>300000 & area_trop_af$value<3000000,]
area_trop_af_l <- area_trop_af[area_trop_af$value>3000000 & area_trop_af$value<60000000,]
area_trop_af_xl <- area_trop_af[area_trop_af$value>60000000,]

area_arid_af_s <- area_arid_af[area_arid_af$value<310000,] ###because value is in hectares
area_arid_af_m <- area_arid_af[area_arid_af$value>300000 & area_arid_af$value<3000000,]
area_arid_af_l <- area_arid_af[area_arid_af$value>3000000 & area_arid_af$value<60000000,]
area_arid_af_xl <- area_arid_af[area_arid_af$value>60000000,]

area_temp_af_s <- area_temp_af[area_temp_af$value<310000,] ###because value is in hectares
area_temp_af_m <- area_temp_af[area_temp_af$value>300000 & area_temp_af$value<3000000,]
area_temp_af_l <- area_temp_af[area_temp_af$value>3000000 & area_temp_af$value<60000000,]
area_temp_af_xl <- area_temp_af[area_temp_af$value>60000000,]

area_cold_af_s <- area_cold_af[area_cold_af$value<310000,] ###because value is in hectares
area_cold_af_m <- area_cold_af[area_cold_af$value>300000 & area_cold_af$value<3000000,]
area_cold_af_l <- area_cold_af[area_cold_af$value>3000000 & area_cold_af$value<60000000,]
area_cold_af_xl <- area_cold_af[area_cold_af$value>60000000,]

area_pola_af_s <- area_pola_af[area_pola_af$value<310000,] ###because value is in hectares
area_pola_af_m <- area_pola_af[area_pola_af$value>300000 & area_pola_af$value<3000000,]
area_pola_af_l <- area_pola_af[area_pola_af$value>3000000 & area_pola_af$value<60000000,]
area_pola_af_xl <- area_pola_af[area_pola_af$value>60000000,]

##eurasia
area_trop_eu_s <- area_trop_eu[area_trop_eu$value<310000,] ###because value is in hectares
area_trop_eu_m <- area_trop_eu[area_trop_eu$value>300000 & area_trop_eu$value<3000000,]
area_trop_eu_l <- area_trop_eu[area_trop_eu$value>3000000 & area_trop_eu$value<60000000,]
area_trop_eu_xl <- area_trop_eu[area_trop_eu$value>60000000,]

area_arid_eu_s <- area_arid_eu[area_arid_eu$value<310000,] ###because value is in hectares
area_arid_eu_m <- area_arid_eu[area_arid_eu$value>300000 & area_arid_eu$value<3000000,]
area_arid_eu_l <- area_arid_eu[area_arid_eu$value>3000000 & area_arid_eu$value<60000000,]
area_arid_eu_xl <- area_arid_eu[area_arid_eu$value>60000000,]

area_temp_eu_s <- area_temp_eu[area_temp_eu$value<310000,] ###because value is in hectares
area_temp_eu_m <- area_temp_eu[area_temp_eu$value>300000 & area_temp_eu$value<3000000,]
area_temp_eu_l <- area_temp_eu[area_temp_eu$value>3000000 & area_temp_eu$value<60000000,]
area_temp_eu_xl <- area_temp_eu[area_temp_eu$value>60000000,]

area_cold_eu_s <- area_cold_eu[area_cold_eu$value<310000,] ###because value is in hectares
area_cold_eu_m <- area_cold_eu[area_cold_eu$value>300000 & area_cold_eu$value<3000000,]
area_cold_eu_l <- area_cold_eu[area_cold_eu$value>3000000 & area_cold_eu$value<60000000,]
area_cold_eu_xl <- area_cold_eu[area_cold_eu$value>60000000,]

area_pola_eu_s <- area_pola_eu[area_pola_eu$value<310000,] ###because value is in hectares
area_pola_eu_m <- area_pola_eu[area_pola_eu$value>300000 & area_pola_eu$value<3000000,]
area_pola_eu_l <- area_pola_eu[area_pola_eu$value>3000000 & area_pola_eu$value<60000000,]
area_pola_eu_xl <- area_pola_eu[area_pola_eu$value>60000000,]

####Number of patches

number_patches <- count(area_metrics, vars = c("layer", "class"))
number_patches$class <- as.factor(number_patches$class)
number_America <- count(area_America, vars = c("layer", "class"))
number_Africa <- count(area_Africa, vars = c("layer", "class"))
number_Eurasia <- count(area_Eurasia, vars = c("layer", "class"))

# beanplot(freq~class, data = number_patches, col=list(c("#f9d14a","#f9d14a","#f9d14a","#f9d14a"),
#                                                      c("#ab3329","#ab3329","#ab3329","#ab3329"),
#                                                      c("#ed968c","#ed968c","#ed968c","#ed968c"),
#                                                      c("#7c4b73","#7c4b73","#7c4b73","#7c4b73"),
#                                                      c("#88a0dc","#88a0dc","#88a0dc","#88a0dc")),
#          overallline = "median",
#          method = "overplot", beanlines = "median", what = c(0,1,1,0), 
#          horizontal = F,frame.plot = F, border="gray48", names = climnames, log="", ylim=c(0,300))

ggplot(number_patches, aes(y=freq, group=class, color=class, fill=class)) + 
  geom_boxplot( outlier.shape=8,
               outlier.size=1) + ylab("Number of patches") +
  scale_color_manual(values = c("#f99900","#4b0000", "#D6666F","#6B0082","#000070" ))+
  scale_fill_manual(values = c("#f9d14a","#ab3239", "#ed968c","#7c4b73","#88a0dc" )) +
  theme_classic() + theme(legend.position="top",panel.border = element_blank(), axis.text.x=element_blank()) +
  geom_hline(yintercept=200, linetype="dashed", color = "lightgrey") +
  geom_hline(yintercept=100, linetype="dashed", color = "lightgrey")
quartz()

kkk <- number_patches[number_patches$class==5,]
mean(np_trop_ww_m$freq)

np_trop_ww_s <- count(area_trop_ww_s, vars = c("layer", "class"))
np_trop_ww_m <- count(area_trop_ww_m, vars = c("layer", "class"))
np_trop_ww_l <- count(area_trop_ww_l, vars = c("layer", "class"))
np_trop_ww_xl <- count(area_trop_ww_xl, vars = c("layer", "class"))
NT_trop_ww_s <- sum(np_trop_ww_s$freq)
NT_trop_ww_m <- sum(np_trop_ww_m$freq)
NT_trop_ww_l <- sum(np_trop_ww_l$freq)
NT_trop_ww_xl <- sum(np_trop_ww_xl$freq)

np_arid_ww_s <- count(area_arid_ww_s, vars = c("layer", "class"))
np_arid_ww_m <- count(area_arid_ww_m, vars = c("layer", "class"))
np_arid_ww_l <- count(area_arid_ww_l, vars = c("layer", "class"))
np_arid_ww_xl <- count(area_arid_ww_xl, vars = c("layer", "class"))
NT_arid_ww_s <- sum(np_arid_ww_s$freq)
NT_arid_ww_m <- sum(np_arid_ww_m$freq)
NT_arid_ww_l <- sum(np_arid_ww_l$freq)
NT_arid_ww_xl <- sum(np_arid_ww_xl$freq)

np_temp_ww_s <- count(area_temp_ww_s, vars = c("layer", "class"))
np_temp_ww_m <- count(area_temp_ww_m, vars = c("layer", "class"))
np_temp_ww_l <- count(area_temp_ww_l, vars = c("layer", "class"))
np_temp_ww_xl <- count(area_temp_ww_xl, vars = c("layer", "class"))
NT_temp_ww_s <- sum(np_temp_ww_s$freq)
NT_temp_ww_m <- sum(np_temp_ww_m$freq)
NT_temp_ww_l <- sum(np_temp_ww_l$freq)
NT_temp_ww_xl <- sum(np_temp_ww_xl$freq)

np_cold_ww_s <- count(area_cold_ww_s, vars = c("layer", "class"))
np_cold_ww_m <- count(area_cold_ww_m, vars = c("layer", "class"))
np_cold_ww_l <- count(area_cold_ww_l, vars = c("layer", "class"))
np_cold_ww_xl <- count(area_cold_ww_xl, vars = c("layer", "class"))
NT_cold_ww_s <- sum(np_cold_ww_s$freq)
NT_cold_ww_m <- sum(np_cold_ww_m$freq)
NT_cold_ww_l <- sum(np_cold_ww_l$freq)
NT_cold_ww_xl <- sum(np_cold_ww_xl$freq)

np_pola_ww_s <- count(area_pola_ww_s, vars = c("layer", "class"))
np_pola_ww_m <- count(area_pola_ww_m, vars = c("layer", "class"))
np_pola_ww_l <- count(area_pola_ww_l, vars = c("layer", "class"))
np_pola_ww_xl <- count(area_pola_ww_xl, vars = c("layer", "class"))
NT_pola_ww_s <- sum(np_pola_ww_s$freq)
NT_pola_ww_m <- sum(np_pola_ww_m$freq)
NT_pola_ww_l <- sum(np_pola_ww_l$freq)
NT_pola_ww_xl <- sum(np_pola_ww_xl$freq)

####America
np_trop_am_s <- count(area_trop_am_s, vars = c("layer", "class"))
np_trop_am_m <- count(area_trop_am_m, vars = c("layer", "class"))
np_trop_am_l <- count(area_trop_am_l, vars = c("layer", "class"))
np_trop_am_xl <- count(area_trop_am_xl, vars = c("layer", "class"))
NT_trop_am_s <- sum(np_trop_am_s$freq)
NT_trop_am_m <- sum(np_trop_am_m$freq)
NT_trop_am_l <- sum(np_trop_am_l$freq)
NT_trop_am_xl <- sum(np_trop_am_xl$freq)

np_arid_am_s <- count(area_arid_am_s, vars = c("layer", "class"))
np_arid_am_m <- count(area_arid_am_m, vars = c("layer", "class"))
np_arid_am_l <- count(area_arid_am_l, vars = c("layer", "class"))
np_arid_am_xl <- count(area_arid_am_xl, vars = c("layer", "class"))
NT_arid_am_s <- sum(np_arid_am_s$freq)
NT_arid_am_m <- sum(np_arid_am_m$freq)
NT_arid_am_l <- sum(np_arid_am_l$freq)
NT_arid_am_xl <- sum(np_arid_am_xl$freq)

np_temp_am_s <- count(area_temp_am_s, vars = c("layer", "class"))
np_temp_am_m <- count(area_temp_am_m, vars = c("layer", "class"))
np_temp_am_l <- count(area_temp_am_l, vars = c("layer", "class"))
np_temp_am_xl <- count(area_temp_am_xl, vars = c("layer", "class"))
NT_temp_am_s <- sum(np_temp_am_s$freq)
NT_temp_am_m <- sum(np_temp_am_m$freq)
NT_temp_am_l <- sum(np_temp_am_l$freq)
NT_temp_am_xl <- sum(np_temp_am_xl$freq)

np_cold_am_s <- count(area_cold_am_s, vars = c("layer", "class"))
np_cold_am_m <- count(area_cold_am_m, vars = c("layer", "class"))
np_cold_am_l <- count(area_cold_am_l, vars = c("layer", "class"))
np_cold_am_xl <- count(area_cold_am_xl, vars = c("layer", "class"))
NT_cold_am_s <- sum(np_cold_am_s$freq)
NT_cold_am_m <- sum(np_cold_am_m$freq)
NT_cold_am_l <- sum(np_cold_am_l$freq)
NT_cold_am_xl <- sum(np_cold_am_xl$freq)

np_pola_am_s <- count(area_pola_am_s, vars = c("layer", "class"))
np_pola_am_m <- count(area_pola_am_m, vars = c("layer", "class"))
np_pola_am_l <- count(area_pola_am_l, vars = c("layer", "class"))
np_pola_am_xl <- count(area_pola_am_xl, vars = c("layer", "class"))
NT_pola_am_s <- sum(np_pola_am_s$freq)
NT_pola_am_m <- sum(np_pola_am_m$freq)
NT_pola_am_l <- sum(np_pola_am_l$freq)
NT_pola_am_xl <- sum(np_pola_am_xl$freq)

###Africa
np_trop_af_s <- count(area_trop_af_s, vars = c("layer", "class"))
np_trop_af_m <- count(area_trop_af_m, vars = c("layer", "class"))
np_trop_af_l <- count(area_trop_af_l, vars = c("layer", "class"))
np_trop_af_xl <- count(area_trop_af_xl, vars = c("layer", "class"))
NT_trop_af_s <- sum(np_trop_af_s$freq)
NT_trop_af_m <- sum(np_trop_af_m$freq)
NT_trop_af_l <- sum(np_trop_af_l$freq)
NT_trop_af_xl <- sum(np_trop_af_xl$freq)

np_arid_af_s <- count(area_arid_af_s, vars = c("layer", "class"))
np_arid_af_m <- count(area_arid_af_m, vars = c("layer", "class"))
np_arid_af_l <- count(area_arid_af_l, vars = c("layer", "class"))
np_arid_af_xl <- count(area_arid_af_xl, vars = c("layer", "class"))
NT_arid_af_s <- sum(np_arid_af_s$freq)
NT_arid_af_m <- sum(np_arid_af_m$freq)
NT_arid_af_l <- sum(np_arid_af_l$freq)
NT_arid_af_xl <- sum(np_arid_af_xl$freq)

np_temp_af_s <- count(area_temp_af_s, vars = c("layer", "class"))
np_temp_af_m <- count(area_temp_af_m, vars = c("layer", "class"))
np_temp_af_l <- count(area_temp_af_l, vars = c("layer", "class"))
np_temp_af_xl <- count(area_temp_af_xl, vars = c("layer", "class"))
NT_temp_af_s <- sum(np_temp_af_s$freq)
NT_temp_af_m <- sum(np_temp_af_m$freq)
NT_temp_af_l <- sum(np_temp_af_l$freq)
NT_temp_af_xl <- sum(np_temp_af_xl$freq)

np_cold_af_s <- count(area_cold_af_s, vars = c("layer", "class"))
np_cold_af_m <- count(area_cold_af_m, vars = c("layer", "class"))
np_cold_af_l <- count(area_cold_af_l, vars = c("layer", "class"))
np_cold_af_xl <- count(area_cold_af_xl, vars = c("layer", "class"))
NT_cold_af_s <- sum(np_cold_af_s$freq)
NT_cold_af_m <- sum(np_cold_af_m$freq)
NT_cold_af_l <- sum(np_cold_af_l$freq)
NT_cold_af_xl <- sum(np_cold_af_xl$freq)

np_pola_af_s <- count(area_pola_af_s, vars = c("layer", "class"))
np_pola_af_m <- count(area_pola_af_m, vars = c("layer", "class"))
np_pola_af_l <- count(area_pola_af_l, vars = c("layer", "class"))
np_pola_af_xl <- count(area_pola_af_xl, vars = c("layer", "class"))
NT_pola_af_s <- sum(np_pola_af_s$freq)
NT_pola_af_m <- sum(np_pola_af_m$freq)
NT_pola_af_l <- sum(np_pola_af_l$freq)
NT_pola_af_xl <- sum(np_pola_af_xl$freq)

###eurasia
np_trop_eu_s <- count(area_trop_eu_s, vars = c("layer", "class"))
np_trop_eu_m <- count(area_trop_eu_m, vars = c("layer", "class"))
np_trop_eu_l <- count(area_trop_eu_l, vars = c("layer", "class"))
np_trop_eu_xl <- count(area_trop_eu_xl, vars = c("layer", "class"))
NT_trop_eu_s <- sum(np_trop_eu_s$freq)
NT_trop_eu_m <- sum(np_trop_eu_m$freq)
NT_trop_eu_l <- sum(np_trop_eu_l$freq)
NT_trop_eu_xl <- sum(np_trop_eu_xl$freq)

mean(np_arid_am_m$freq)
np_arid_eu_s <- count(area_arid_eu_s, vars = c("layer", "class"))
np_arid_eu_m <- count(area_arid_eu_m, vars = c("layer", "class"))
np_arid_eu_l <- count(area_arid_eu_l, vars = c("layer", "class"))
np_arid_eu_xl <- count(area_arid_eu_xl, vars = c("layer", "class"))
NT_arid_eu_s <- sum(np_arid_eu_s$freq)
NT_arid_eu_m <- sum(np_arid_eu_m$freq)
NT_arid_eu_l <- sum(np_arid_eu_l$freq)
NT_arid_eu_xl <- sum(np_arid_eu_xl$freq)

np_temp_eu_s <- count(area_temp_eu_s, vars = c("layer", "class"))
np_temp_eu_m <- count(area_temp_eu_m, vars = c("layer", "class"))
np_temp_eu_l <- count(area_temp_eu_l, vars = c("layer", "class"))
np_temp_eu_xl <- count(area_temp_eu_xl, vars = c("layer", "class"))
NT_temp_eu_s <- sum(np_temp_eu_s$freq)
NT_temp_eu_m <- sum(np_temp_eu_m$freq)
NT_temp_eu_l <- sum(np_temp_eu_l$freq)
NT_temp_eu_xl <- sum(np_temp_eu_xl$freq)

np_cold_eu_s <- count(area_cold_eu_s, vars = c("layer", "class"))
np_cold_eu_m <- count(area_cold_eu_m, vars = c("layer", "class"))
np_cold_eu_l <- count(area_cold_eu_l, vars = c("layer", "class"))
np_cold_eu_xl <- count(area_cold_eu_xl, vars = c("layer", "class"))
NT_cold_eu_s <- sum(np_cold_eu_s$freq)
NT_cold_eu_m <- sum(np_cold_eu_m$freq)
NT_cold_eu_l <- sum(np_cold_eu_l$freq)
NT_cold_eu_xl <- sum(np_cold_eu_xl$freq)

np_pola_eu_s <- count(area_pola_eu_s, vars = c("layer", "class"))
np_pola_eu_m <- count(area_pola_eu_m, vars = c("layer", "class"))
np_pola_eu_l <- count(area_pola_eu_l, vars = c("layer", "class"))
np_pola_eu_xl <- count(area_pola_eu_xl, vars = c("layer", "class"))
NT_pola_eu_s <- sum(np_pola_eu_s$freq)
NT_pola_eu_m <- sum(np_pola_eu_m$freq)
NT_pola_eu_l <- sum(np_pola_eu_l$freq)
NT_pola_eu_xl <- sum(np_pola_eu_xl$freq)


S_patches <- rbind(np_trop_ww_s,np_arid_ww_s,np_temp_ww_s,np_cold_ww_s,np_pola_ww_s)
S_patches$class <- as.factor(S_patches$class)
M_patches <- rbind(np_trop_ww_m,np_arid_ww_m,np_temp_ww_m,np_cold_ww_m,np_pola_ww_m)
M_patches$class <- as.factor(M_patches$class)
L_patches <- rbind(np_trop_ww_l,np_arid_ww_l,np_temp_ww_l,np_cold_ww_l,np_pola_ww_l)
L_patches$class <- as.factor(L_patches$class)
XL_patches <- rbind(np_trop_ww_xl,np_arid_ww_xl,np_temp_ww_xl,np_cold_ww_xl,np_pola_ww_xl)
XL_patches$class <- as.factor(XL_patches$class)

S_patches_am <- rbind(np_trop_am_s,np_arid_am_s,np_temp_am_s,np_cold_am_s,np_pola_am_s)
S_patches_am$class <- as.factor(S_patches_am$class)
M_patches_am <- rbind(np_trop_am_m,np_arid_am_m,np_temp_am_m,np_cold_am_m,np_pola_am_m)
M_patches_am$class <- as.factor(M_patches_am$class)
L_patches_am <- rbind(np_trop_am_l,np_arid_am_l,np_temp_am_l,np_cold_am_l,np_pola_am_l)
L_patches_am$class <- as.factor(L_patches_am$class)
XL_patches_am <- rbind(np_trop_am_xl,np_arid_am_xl,np_temp_am_xl,np_cold_am_xl,np_pola_am_xl)
XL_patches_am$class <- as.factor(XL_patches_am$class)

S_patches_af <- rbind(np_trop_af_s,np_arid_af_s,np_temp_af_s,np_cold_af_s,np_pola_af_s)
S_patches_af$class <- as.factor(S_patches_af$class)
M_patches_af <- rbind(np_trop_af_m,np_arid_af_m,np_temp_af_m,np_cold_af_m,np_pola_af_m)
M_patches_af$class <- as.factor(M_patches_af$class)
L_patches_af <- rbind(np_trop_af_l,np_arid_af_l,np_temp_af_l,np_cold_af_l,np_pola_af_l)
L_patches_af$class <- as.factor(L_patches_af$class)
XL_patches_af <- rbind(np_trop_af_xl,np_arid_af_xl,np_temp_af_xl,np_cold_af_xl,np_pola_af_xl)
XL_patches_af$class <- as.factor(XL_patches_af$class)

S_patches_eu <- rbind(np_trop_eu_s,np_arid_eu_s,np_temp_eu_s,np_cold_eu_s,np_pola_eu_s)
S_patches_eu$class <- as.factor(S_patches_eu$class)
M_patches_eu <- rbind(np_trop_eu_m,np_arid_eu_m,np_temp_eu_m,np_cold_eu_m,np_pola_eu_m)
M_patches_eu$class <- as.factor(M_patches_eu$class)
L_patches_eu <- rbind(np_trop_eu_l,np_arid_eu_l,np_temp_eu_l,np_cold_eu_l,np_pola_eu_l)
L_patches_eu$class <- as.factor(L_patches_eu$class)
XL_patches_eu <- rbind(np_trop_eu_xl,np_arid_eu_xl,np_temp_eu_xl,np_cold_eu_xl,np_pola_eu_xl)
XL_patches_eu$class <- as.factor(XL_patches_eu$class)

ggplot(XL_patches_eu, aes(y=freq, group=class, color=class, fill=class)) + 
  geom_boxplot( outlier.shape=8,
                outlier.size=1) + ylab("Number of patches") +
  scale_color_manual(values = c("#f99900","#4b0000", "#D6666F","#6B0082","#000070" ))+
  scale_fill_manual(values = c("#f9d14a","#ab3239", "#ed968c","#7c4b73","#88a0dc" )) + ylim(c(0,50)) +
  geom_hline(yintercept=50, linetype="dashed", color = "lightgrey") +
  geom_hline(yintercept=150, linetype="dashed", color = "lightgrey") +
  theme_classic() + theme(legend.position="top",panel.border = element_blank(), axis.text.x=element_blank(), axis.ticks = element_blank())

quartz(width=9/2.54)
L_patches_af$class

#####ANOVA                                                                                                    data=S_patches)))[[1]][["Pr(>F)"]][[1]])
anovaOneW(formula = freq ~ class,
          data = XL_patches_eu,
          welchs = TRUE,
          eqv = TRUE,
          phSig = TRUE,
          phMethod = 'gamesHowell')

###################
###Calculating net fragmentation
dif_areas_ww<- as.data.frame(cbind(number_patches$layer[1:length(number_patches$layer)-1],
                                   number_patches[1:length(number_patches$freq)-1,2] - number_patches[2:length(number_patches$freq),2]))
dif_areas_am<- as.data.frame(cbind(number_America$layer[1:length(number_America$layer)-1],
                                   number_America[1:length(number_America$freq)-1,2] - number_America[2:length(number_America$freq),2])) 
dif_areas_af<- as.data.frame(cbind(number_Africa$layer[1:length(number_Africa$layer)-1],
                                   number_Africa[1:length(number_Africa$freq)-1,2] - number_Africa[2:length(number_Africa$freq),2])) 
dif_areas_eu<- as.data.frame(cbind(number_Eurasia$layer[1:length(number_Eurasia$layer)-1],
                                   number_Eurasia[1:length(number_Eurasia$freq)-1,2] - number_Eurasia[2:length(number_Eurasia$freq),2])) 
colnames(dif_areas_ww) <- c("Time", "d_patches")
colnames(dif_areas_am) <- c("Time", "d_patches")
colnames(dif_areas_af) <- c("Time", "d_patches")
colnames(dif_areas_eu) <- c("Time", "d_patches")

#####Worldwide
dif_areas_trop_ww_s<- as.data.frame(cbind(np_trop_ww_s$layer[1:length(np_trop_ww_s$layer)-1],
                                          np_trop_ww_s[1:length(np_trop_ww_s$freq)-1,3] - np_trop_ww_s[2:length(np_trop_ww_s$freq),3])) 
dif_areas_trop_ww_m<- as.data.frame(cbind(np_trop_ww_m$layer[1:length(np_trop_ww_m$layer)-1],
                                          np_trop_ww_m[1:length(np_trop_ww_m$freq)-1,3] - np_trop_ww_m[2:length(np_trop_ww_m$freq),3])) 
dif_areas_trop_ww_l<- as.data.frame(cbind(np_trop_ww_l$layer[1:length(np_trop_ww_l$layer)-1],
                                          np_trop_ww_l[1:length(np_trop_ww_l$freq)-1,3] - np_trop_ww_l[2:length(np_trop_ww_l$freq),3])) 
dif_areas_trop_ww_xl<- as.data.frame(cbind(np_trop_ww_xl$layer[1:length(np_trop_ww_xl$layer)-1],
                                           np_trop_ww_xl[1:length(np_trop_ww_xl$freq)-1,3] - np_trop_ww_xl[2:length(np_trop_ww_xl$freq),3])) 

dif_areas_arid_ww_s<- as.data.frame(cbind(np_arid_ww_s$layer[1:length(np_arid_ww_s$layer)-1],
                                          np_arid_ww_s[1:length(np_arid_ww_s$freq)-1,3] - np_arid_ww_s[2:length(np_arid_ww_s$freq),3])) 
dif_areas_arid_ww_m<- as.data.frame(cbind(np_arid_ww_m$layer[1:length(np_arid_ww_m$layer)-1],
                                          np_arid_ww_m[1:length(np_arid_ww_m$freq)-1,3] - np_arid_ww_m[2:length(np_arid_ww_m$freq),3])) 
dif_areas_arid_ww_l<- as.data.frame(cbind(np_arid_ww_l$layer[1:length(np_arid_ww_l$layer)-1],
                                          np_arid_ww_l[1:length(np_arid_ww_l$freq)-1,3] - np_arid_ww_l[2:length(np_arid_ww_l$freq),3])) 
dif_areas_arid_ww_xl<- as.data.frame(cbind(np_arid_ww_xl$layer[1:length(np_arid_ww_xl$layer)-1],
                                           np_arid_ww_xl[1:length(np_arid_ww_xl$freq)-1,3] - np_arid_ww_xl[2:length(np_arid_ww_xl$freq),3])) 

dif_areas_temp_ww_s<- as.data.frame(cbind(np_temp_ww_s$layer[1:length(np_temp_ww_s$layer)-1],
                                          np_temp_ww_s[1:length(np_temp_ww_s$freq)-1,3] - np_temp_ww_s[2:length(np_temp_ww_s$freq),3])) 
dif_areas_temp_ww_m<- as.data.frame(cbind(np_temp_ww_m$layer[1:length(np_temp_ww_m$layer)-1],
                                          np_temp_ww_m[1:length(np_temp_ww_m$freq)-1,3] - np_temp_ww_m[2:length(np_temp_ww_m$freq),3])) 
dif_areas_temp_ww_l<- as.data.frame(cbind(np_temp_ww_l$layer[1:length(np_temp_ww_l$layer)-1],
                                          np_temp_ww_l[1:length(np_temp_ww_l$freq)-1,3] - np_temp_ww_l[2:length(np_temp_ww_l$freq),3])) 
dif_areas_temp_ww_xl<- as.data.frame(cbind(np_temp_ww_xl$layer[1:length(np_temp_ww_xl$layer)-1],
                                           np_temp_ww_xl[1:length(np_temp_ww_xl$freq)-1,3] - np_temp_ww_xl[2:length(np_temp_ww_xl$freq),3])) 

dif_areas_cold_ww_s<- as.data.frame(cbind(np_cold_ww_s$layer[1:length(np_cold_ww_s$layer)-1],
                                          np_cold_ww_s[1:length(np_cold_ww_s$freq)-1,3] - np_cold_ww_s[2:length(np_cold_ww_s$freq),3])) 
dif_areas_cold_ww_m<- as.data.frame(cbind(np_cold_ww_m$layer[1:length(np_cold_ww_m$layer)-1],
                                          np_cold_ww_m[1:length(np_cold_ww_m$freq)-1,3] - np_cold_ww_m[2:length(np_cold_ww_m$freq),3])) 
dif_areas_cold_ww_l<- as.data.frame(cbind(np_cold_ww_l$layer[1:length(np_cold_ww_l$layer)-1],
                                          np_cold_ww_l[1:length(np_cold_ww_l$freq)-1,3] - np_cold_ww_l[2:length(np_cold_ww_l$freq),3])) 
dif_areas_cold_ww_xl<- as.data.frame(cbind(np_cold_ww_xl$layer[1:length(np_cold_ww_xl$layer)-1],
                                           np_cold_ww_xl[1:length(np_cold_ww_xl$freq)-1,3] - np_cold_ww_xl[2:length(np_cold_ww_xl$freq),3])) 

dif_areas_pola_ww_s<- as.data.frame(cbind(np_pola_ww_s$layer[1:length(np_pola_ww_s$layer)-1],
                                          np_pola_ww_s[1:length(np_pola_ww_s$freq)-1,3] - np_pola_ww_s[2:length(np_pola_ww_s$freq),3])) 
dif_areas_pola_ww_m<- as.data.frame(cbind(np_pola_ww_m$layer[1:length(np_pola_ww_m$layer)-1],
                                          np_pola_ww_m[1:length(np_pola_ww_m$freq)-1,3] - np_pola_ww_m[2:length(np_pola_ww_m$freq),3])) 
dif_areas_pola_ww_l<- as.data.frame(cbind(np_pola_ww_l$layer[1:length(np_pola_ww_l$layer)-1],
                                          np_pola_ww_l[1:length(np_pola_ww_l$freq)-1,3] - np_pola_ww_l[2:length(np_pola_ww_l$freq),3])) 
dif_areas_pola_ww_xl<- as.data.frame(cbind(np_pola_ww_xl$layer[1:length(np_pola_ww_xl$layer)-1],
                                           np_pola_ww_xl[1:length(np_pola_ww_xl$freq)-1,3] - np_pola_ww_xl[2:length(np_pola_ww_xl$freq),3])) 
colnames(dif_areas_trop_ww_s) <- c("Time", "d_patches")
colnames(dif_areas_trop_ww_m) <- c("Time", "d_patches")
colnames(dif_areas_trop_ww_l) <- c("Time", "d_patches")
colnames(dif_areas_trop_ww_xl) <- c("Time", "d_patches")
colnames(dif_areas_arid_ww_s) <- c("Time", "d_patches")
colnames(dif_areas_arid_ww_m) <- c("Time", "d_patches")
colnames(dif_areas_arid_ww_l) <- c("Time", "d_patches")
colnames(dif_areas_arid_ww_xl) <- c("Time", "d_patches")
colnames(dif_areas_temp_ww_s) <- c("Time", "d_patches")
colnames(dif_areas_temp_ww_m) <- c("Time", "d_patches")
colnames(dif_areas_temp_ww_l) <- c("Time", "d_patches")
colnames(dif_areas_temp_ww_xl) <- c("Time", "d_patches")
colnames(dif_areas_cold_ww_s) <- c("Time", "d_patches")
colnames(dif_areas_cold_ww_m) <- c("Time", "d_patches")
colnames(dif_areas_cold_ww_l) <- c("Time", "d_patches")
colnames(dif_areas_cold_ww_xl) <- c("Time", "d_patches")
colnames(dif_areas_pola_ww_s) <- c("Time", "d_patches")
colnames(dif_areas_pola_ww_m) <- c("Time", "d_patches")
colnames(dif_areas_pola_ww_l) <- c("Time", "d_patches")
colnames(dif_areas_pola_ww_xl) <- c("Time", "d_patches")

######America
dif_areas_trop_am_s<- as.data.frame(cbind(np_trop_am_s$layer[1:length(np_trop_am_s$layer)-1],
                                          np_trop_am_s[1:length(np_trop_am_s$freq)-1,3] - np_trop_am_s[2:length(np_trop_am_s$freq),3])) 
dif_areas_trop_am_m<- as.data.frame(cbind(np_trop_am_m$layer[1:length(np_trop_am_m$layer)-1],
                                          np_trop_am_m[1:length(np_trop_am_m$freq)-1,3] - np_trop_am_m[2:length(np_trop_am_m$freq),3])) 
dif_areas_trop_am_l<- as.data.frame(cbind(np_trop_am_l$layer[1:length(np_trop_am_l$layer)-1],
                                          np_trop_am_l[1:length(np_trop_am_l$freq)-1,3] - np_trop_am_l[2:length(np_trop_am_l$freq),3])) 
dif_areas_trop_am_xl<- as.data.frame(cbind(np_trop_am_xl$layer[1:length(np_trop_am_xl$layer)-1],
                                           np_trop_am_xl[1:length(np_trop_am_xl$freq)-1,3] - np_trop_am_xl[2:length(np_trop_am_xl$freq),3])) 

dif_areas_arid_am_s<- as.data.frame(cbind(np_arid_am_s$layer[1:length(np_arid_am_s$layer)-1],
                                          np_arid_am_s[1:length(np_arid_am_s$freq)-1,3] - np_arid_am_s[2:length(np_arid_am_s$freq),3])) 
dif_areas_arid_am_m<- as.data.frame(cbind(np_arid_am_m$layer[1:length(np_arid_am_m$layer)-1],
                                          np_arid_am_m[1:length(np_arid_am_m$freq)-1,3] - np_arid_am_m[2:length(np_arid_am_m$freq),3])) 
dif_areas_arid_am_l<- as.data.frame(cbind(np_arid_am_l$layer[1:length(np_arid_am_l$layer)-1],
                                          np_arid_am_l[1:length(np_arid_am_l$freq)-1,3] - np_arid_am_l[2:length(np_arid_am_l$freq),3])) 
dif_areas_arid_am_xl<- as.data.frame(cbind(np_arid_am_xl$layer[1:length(np_arid_am_xl$layer)-1],
                                           np_arid_am_xl[1:length(np_arid_am_xl$freq)-1,3] - np_arid_am_xl[2:length(np_arid_am_xl$freq),3])) 

dif_areas_temp_am_s<- as.data.frame(cbind(np_temp_am_s$layer[1:length(np_temp_am_s$layer)-1],
                                          np_temp_am_s[1:length(np_temp_am_s$freq)-1,3] - np_temp_am_s[2:length(np_temp_am_s$freq),3])) 
dif_areas_temp_am_m<- as.data.frame(cbind(np_temp_am_m$layer[1:length(np_temp_am_m$layer)-1],
                                          np_temp_am_m[1:length(np_temp_am_m$freq)-1,3] - np_temp_am_m[2:length(np_temp_am_m$freq),3])) 
dif_areas_temp_am_l<- as.data.frame(cbind(np_temp_am_l$layer[1:length(np_temp_am_l$layer)-1],
                                          np_temp_am_l[1:length(np_temp_am_l$freq)-1,3] - np_temp_am_l[2:length(np_temp_am_l$freq),3])) 
dif_areas_temp_am_xl<- as.data.frame(cbind(np_temp_am_xl$layer[1:length(np_temp_am_xl$layer)-1],
                                           np_temp_am_xl[1:length(np_temp_am_xl$freq)-1,3] - np_temp_am_xl[2:length(np_temp_am_xl$freq),3])) 

dif_areas_cold_am_s<- as.data.frame(cbind(np_cold_am_s$layer[1:length(np_cold_am_s$layer)-1],
                                          np_cold_am_s[1:length(np_cold_am_s$freq)-1,3] - np_cold_am_s[2:length(np_cold_am_s$freq),3])) 
dif_areas_cold_am_m<- as.data.frame(cbind(np_cold_am_m$layer[1:length(np_cold_am_m$layer)-1],
                                          np_cold_am_m[1:length(np_cold_am_m$freq)-1,3] - np_cold_am_m[2:length(np_cold_am_m$freq),3])) 
dif_areas_cold_am_l<- as.data.frame(cbind(np_cold_am_l$layer[1:length(np_cold_am_l$layer)-1],
                                          np_cold_am_l[1:length(np_cold_am_l$freq)-1,3] - np_cold_am_l[2:length(np_cold_am_l$freq),3])) 
dif_areas_cold_am_xl<- as.data.frame(cbind(np_cold_am_xl$layer[1:length(np_cold_am_xl$layer)-1],
                                           np_cold_am_xl[1:length(np_cold_am_xl$freq)-1,3] - np_cold_am_xl[2:length(np_cold_am_xl$freq),3])) 

dif_areas_pola_am_s<- as.data.frame(cbind(np_pola_am_s$layer[1:length(np_pola_am_s$layer)-1],
                                          np_pola_am_s[1:length(np_pola_am_s$freq)-1,3] - np_pola_am_s[2:length(np_pola_am_s$freq),3])) 
dif_areas_pola_am_m<- as.data.frame(cbind(np_pola_am_m$layer[1:length(np_pola_am_m$layer)-1],
                                          np_pola_am_m[1:length(np_pola_am_m$freq)-1,3] - np_pola_am_m[2:length(np_pola_am_m$freq),3])) 
dif_areas_pola_am_l<- as.data.frame(cbind(np_pola_am_l$layer[1:length(np_pola_am_l$layer)-1],
                                          np_pola_am_l[1:length(np_pola_am_l$freq)-1,3] - np_pola_am_l[2:length(np_pola_am_l$freq),3])) 
dif_areas_pola_am_xl<- as.data.frame(cbind(np_pola_am_xl$layer[1:length(np_pola_am_xl$layer)-1],
                                           np_pola_am_xl[1:length(np_pola_am_xl$freq)-1,3] - np_pola_am_xl[2:length(np_pola_am_xl$freq),3])) 
colnames(dif_areas_trop_am_s) <- c("Time", "d_patches")
colnames(dif_areas_trop_am_m) <- c("Time", "d_patches")
colnames(dif_areas_trop_am_l) <- c("Time", "d_patches")
colnames(dif_areas_trop_am_xl) <- c("Time", "d_patches")
colnames(dif_areas_arid_am_s) <- c("Time", "d_patches")
colnames(dif_areas_arid_am_m) <- c("Time", "d_patches")
colnames(dif_areas_arid_am_l) <- c("Time", "d_patches")
colnames(dif_areas_arid_am_xl) <- c("Time", "d_patches")
colnames(dif_areas_temp_am_s) <- c("Time", "d_patches")
colnames(dif_areas_temp_am_m) <- c("Time", "d_patches")
colnames(dif_areas_temp_am_l) <- c("Time", "d_patches")
colnames(dif_areas_temp_am_xl) <- c("Time", "d_patches")
colnames(dif_areas_cold_am_s) <- c("Time", "d_patches")
colnames(dif_areas_cold_am_m) <- c("Time", "d_patches")
colnames(dif_areas_cold_am_l) <- c("Time", "d_patches")
colnames(dif_areas_cold_am_xl) <- c("Time", "d_patches")
colnames(dif_areas_pola_am_s) <- c("Time", "d_patches")
colnames(dif_areas_pola_am_m) <- c("Time", "d_patches")
colnames(dif_areas_pola_am_l) <- c("Time", "d_patches")
colnames(dif_areas_pola_am_xl) <- c("Time", "d_patches")

###Africa
dif_areas_trop_af_s<- as.data.frame(cbind(np_trop_af_s$layer[1:length(np_trop_af_s$layer)-1],
                                          np_trop_af_s[1:length(np_trop_af_s$freq)-1,3] - np_trop_af_s[2:length(np_trop_af_s$freq),3])) 
dif_areas_trop_af_m<- as.data.frame(cbind(np_trop_af_m$layer[1:length(np_trop_af_m$layer)-1],
                                          np_trop_af_m[1:length(np_trop_af_m$freq)-1,3] - np_trop_af_m[2:length(np_trop_af_m$freq),3])) 
dif_areas_trop_af_l<- as.data.frame(cbind(np_trop_af_l$layer[1:length(np_trop_af_l$layer)-1],
                                          np_trop_af_l[1:length(np_trop_af_l$freq)-1,3] - np_trop_af_l[2:length(np_trop_af_l$freq),3])) 
dif_areas_trop_af_xl<- as.data.frame(cbind(np_trop_af_xl$layer[1:length(np_trop_af_xl$layer)-1],
                                           np_trop_af_xl[1:length(np_trop_af_xl$freq)-1,3] - np_trop_af_xl[2:length(np_trop_af_xl$freq),3])) 

dif_areas_arid_af_s<- as.data.frame(cbind(np_arid_af_s$layer[1:length(np_arid_af_s$layer)-1],
                                          np_arid_af_s[1:length(np_arid_af_s$freq)-1,3] - np_arid_af_s[2:length(np_arid_af_s$freq),3])) 
dif_areas_arid_af_m<- as.data.frame(cbind(np_arid_af_m$layer[1:length(np_arid_af_m$layer)-1],
                                          np_arid_af_m[1:length(np_arid_af_m$freq)-1,3] - np_arid_af_m[2:length(np_arid_af_m$freq),3])) 
dif_areas_arid_af_l<- as.data.frame(cbind(np_arid_af_l$layer[1:length(np_arid_af_l$layer)-1],
                                          np_arid_af_l[1:length(np_arid_af_l$freq)-1,3] - np_arid_af_l[2:length(np_arid_af_l$freq),3])) 
dif_areas_arid_af_xl<- as.data.frame(cbind(np_arid_af_xl$layer[1:length(np_arid_af_xl$layer)-1],
                                           np_arid_af_xl[1:length(np_arid_af_xl$freq)-1,3] - np_arid_af_xl[2:length(np_arid_af_xl$freq),3])) 

dif_areas_temp_af_s<- as.data.frame(cbind(np_temp_af_s$layer[1:length(np_temp_af_s$layer)-1],
                                          np_temp_af_s[1:length(np_temp_af_s$freq)-1,3] - np_temp_af_s[2:length(np_temp_af_s$freq),3])) 
dif_areas_temp_af_m<- as.data.frame(cbind(np_temp_af_m$layer[1:length(np_temp_af_m$layer)-1],
                                          np_temp_af_m[1:length(np_temp_af_m$freq)-1,3] - np_temp_af_m[2:length(np_temp_af_m$freq),3])) 
dif_areas_temp_af_l<- as.data.frame(cbind(np_temp_af_l$layer[1:length(np_temp_af_l$layer)-1],
                                          np_temp_af_l[1:length(np_temp_af_l$freq)-1,3] - np_temp_af_l[2:length(np_temp_af_l$freq),3])) 
dif_areas_temp_af_xl<- as.data.frame(cbind(np_temp_af_xl$layer[1:length(np_temp_af_xl$layer)-1],
                                           np_temp_af_xl[1:length(np_temp_af_xl$freq)-1,3] - np_temp_af_xl[2:length(np_temp_af_xl$freq),3])) 

dif_areas_cold_af_s<- as.data.frame(cbind(np_cold_af_s$layer[1:length(np_cold_af_s$layer)-1],
                                          np_cold_af_s[1:length(np_cold_af_s$freq)-1,3] - np_cold_af_s[2:length(np_cold_af_s$freq),3])) 
dif_areas_cold_af_m<- as.data.frame(cbind(np_cold_af_m$layer[1:length(np_cold_af_m$layer)-1],
                                          np_cold_af_m[1:length(np_cold_af_m$freq)-1,3] - np_cold_af_m[2:length(np_cold_af_m$freq),3])) 
dif_areas_cold_af_l<- as.data.frame(cbind(np_cold_af_l$layer[1:length(np_cold_af_l$layer)-1],
                                          np_cold_af_l[1:length(np_cold_af_l$freq)-1,3] - np_cold_af_l[2:length(np_cold_af_l$freq),3])) 
dif_areas_cold_af_xl<- as.data.frame(cbind(np_cold_af_xl$layer[1:length(np_cold_af_xl$layer)-1],
                                           np_cold_af_xl[1:length(np_cold_af_xl$freq)-1,3] - np_cold_af_xl[2:length(np_cold_af_xl$freq),3])) 

dif_areas_pola_af_s<- as.data.frame(cbind(np_pola_af_s$layer[1:length(np_pola_af_s$layer)-1],
                                          np_pola_af_s[1:length(np_pola_af_s$freq)-1,3] - np_pola_af_s[2:length(np_pola_af_s$freq),3])) 
dif_areas_pola_af_m<- as.data.frame(cbind(np_pola_af_m$layer[1:length(np_pola_af_m$layer)-1],
                                          np_pola_af_m[1:length(np_pola_af_m$freq)-1,3] - np_pola_af_m[2:length(np_pola_af_m$freq),3])) 
dif_areas_pola_af_l<- as.data.frame(cbind(np_pola_af_l$layer[1:length(np_pola_af_l$layer)-1],
                                          np_pola_af_l[1:length(np_pola_af_l$freq)-1,3] - np_pola_af_l[2:length(np_pola_af_l$freq),3])) 
dif_areas_pola_af_xl<- as.data.frame(cbind(np_pola_af_xl$layer[1:length(np_pola_af_xl$layer)-1],
                                           np_pola_af_xl[1:length(np_pola_af_xl$freq)-1,3] - np_pola_af_xl[2:length(np_pola_af_xl$freq),3])) 
colnames(dif_areas_trop_af_s) <- c("Time", "d_patches")
colnames(dif_areas_trop_af_m) <- c("Time", "d_patches")
colnames(dif_areas_trop_af_l) <- c("Time", "d_patches")
colnames(dif_areas_trop_af_xl) <- c("Time", "d_patches")
colnames(dif_areas_arid_af_s) <- c("Time", "d_patches")
colnames(dif_areas_arid_af_m) <- c("Time", "d_patches")
colnames(dif_areas_arid_af_l) <- c("Time", "d_patches")
colnames(dif_areas_arid_af_xl) <- c("Time", "d_patches")
colnames(dif_areas_temp_af_s) <- c("Time", "d_patches")
colnames(dif_areas_temp_af_m) <- c("Time", "d_patches")
colnames(dif_areas_temp_af_l) <- c("Time", "d_patches")
colnames(dif_areas_temp_af_xl) <- c("Time", "d_patches")
colnames(dif_areas_cold_af_s) <- c("Time", "d_patches")
colnames(dif_areas_cold_af_m) <- c("Time", "d_patches")
colnames(dif_areas_cold_af_l) <- c("Time", "d_patches")
colnames(dif_areas_cold_af_xl) <- c("Time", "d_patches")
colnames(dif_areas_pola_af_s) <- c("Time", "d_patches")
colnames(dif_areas_pola_af_m) <- c("Time", "d_patches")
colnames(dif_areas_pola_af_l) <- c("Time", "d_patches")
colnames(dif_areas_pola_af_xl) <- c("Time", "d_patches")

#####Eurasia
dif_areas_trop_eu_s<- as.data.frame(cbind(np_trop_eu_s$layer[1:length(np_trop_eu_s$layer)-1],
                                          np_trop_eu_s[1:length(np_trop_eu_s$freq)-1,3] - np_trop_eu_s[2:length(np_trop_eu_s$freq),3])) 
dif_areas_trop_eu_m<- as.data.frame(cbind(np_trop_eu_m$layer[1:length(np_trop_eu_m$layer)-1],
                                          np_trop_eu_m[1:length(np_trop_eu_m$freq)-1,3] - np_trop_eu_m[2:length(np_trop_eu_m$freq),3])) 
dif_areas_trop_eu_l<- as.data.frame(cbind(np_trop_eu_l$layer[1:length(np_trop_eu_l$layer)-1],
                                          np_trop_eu_l[1:length(np_trop_eu_l$freq)-1,3] - np_trop_eu_l[2:length(np_trop_eu_l$freq),3])) 
dif_areas_trop_eu_xl<- as.data.frame(cbind(np_trop_eu_xl$layer[1:length(np_trop_eu_xl$layer)-1],
                                           np_trop_eu_xl[1:length(np_trop_eu_xl$freq)-1,3] - np_trop_eu_xl[2:length(np_trop_eu_xl$freq),3])) 

dif_areas_arid_eu_s<- as.data.frame(cbind(np_arid_eu_s$layer[1:length(np_arid_eu_s$layer)-1],
                                          np_arid_eu_s[1:length(np_arid_eu_s$freq)-1,3] - np_arid_eu_s[2:length(np_arid_eu_s$freq),3])) 
dif_areas_arid_eu_m<- as.data.frame(cbind(np_arid_eu_m$layer[1:length(np_arid_eu_m$layer)-1],
                                          np_arid_eu_m[1:length(np_arid_eu_m$freq)-1,3] - np_arid_eu_m[2:length(np_arid_eu_m$freq),3])) 
dif_areas_arid_eu_l<- as.data.frame(cbind(np_arid_eu_l$layer[1:length(np_arid_eu_l$layer)-1],
                                          np_arid_eu_l[1:length(np_arid_eu_l$freq)-1,3] - np_arid_eu_l[2:length(np_arid_eu_l$freq),3])) 
dif_areas_arid_eu_xl<- as.data.frame(cbind(np_arid_eu_xl$layer[1:length(np_arid_eu_xl$layer)-1],
                                           np_arid_eu_xl[1:length(np_arid_eu_xl$freq)-1,3] - np_arid_eu_xl[2:length(np_arid_eu_xl$freq),3])) 

dif_areas_temp_eu_s<- as.data.frame(cbind(np_temp_eu_s$layer[1:length(np_temp_eu_s$layer)-1],
                                          np_temp_eu_s[1:length(np_temp_eu_s$freq)-1,3] - np_temp_eu_s[2:length(np_temp_eu_s$freq),3])) 
dif_areas_temp_eu_m<- as.data.frame(cbind(np_temp_eu_m$layer[1:length(np_temp_eu_m$layer)-1],
                                          np_temp_eu_m[1:length(np_temp_eu_m$freq)-1,3] - np_temp_eu_m[2:length(np_temp_eu_m$freq),3])) 
dif_areas_temp_eu_l<- as.data.frame(cbind(np_temp_eu_l$layer[1:length(np_temp_eu_l$layer)-1],
                                          np_temp_eu_l[1:length(np_temp_eu_l$freq)-1,3] - np_temp_eu_l[2:length(np_temp_eu_l$freq),3])) 
dif_areas_temp_eu_xl<- as.data.frame(cbind(np_temp_eu_xl$layer[1:length(np_temp_eu_xl$layer)-1],
                                           np_temp_eu_xl[1:length(np_temp_eu_xl$freq)-1,3] - np_temp_eu_xl[2:length(np_temp_eu_xl$freq),3])) 


dif_areas_cold_eu_s<- as.data.frame(cbind(np_cold_eu_s$layer[1:length(np_cold_eu_s$layer)-1],
                                          np_cold_eu_s[1:length(np_cold_eu_s$freq)-1,3] - np_cold_eu_s[2:length(np_cold_eu_s$freq),3])) 
dif_areas_cold_eu_m<- as.data.frame(cbind(np_cold_eu_m$layer[1:length(np_cold_eu_m$layer)-1],
                                          np_cold_eu_m[1:length(np_cold_eu_m$freq)-1,3] - np_cold_eu_m[2:length(np_cold_eu_m$freq),3])) 
dif_areas_cold_eu_l<- as.data.frame(cbind(np_cold_eu_l$layer[1:length(np_cold_eu_l$layer)-1],
                                          np_cold_eu_l[1:length(np_cold_eu_l$freq)-1,3] - np_cold_eu_l[2:length(np_cold_eu_l$freq),3])) 
dif_areas_cold_eu_xl<- as.data.frame(cbind(np_cold_eu_xl$layer[1:length(np_cold_eu_xl$layer)-1],
                                           np_cold_eu_xl[1:length(np_cold_eu_xl$freq)-1,3] - np_cold_eu_xl[2:length(np_cold_eu_xl$freq),3])) 

dif_areas_pola_eu_s<- as.data.frame(cbind(np_pola_eu_s$layer[1:length(np_pola_eu_s$layer)-1],
                                          np_pola_eu_s[1:length(np_pola_eu_s$freq)-1,3] - np_pola_eu_s[2:length(np_pola_eu_s$freq),3])) 
dif_areas_pola_eu_m<- as.data.frame(cbind(np_pola_eu_m$layer[1:length(np_pola_eu_m$layer)-1],
                                          np_pola_eu_m[1:length(np_pola_eu_m$freq)-1,3] - np_pola_eu_m[2:length(np_pola_eu_m$freq),3])) 
dif_areas_pola_eu_l<- as.data.frame(cbind(np_pola_eu_l$layer[1:length(np_pola_eu_l$layer)-1],
                                          np_pola_eu_l[1:length(np_pola_eu_l$freq)-1,3] - np_pola_eu_l[2:length(np_pola_eu_l$freq),3])) 
dif_areas_pola_eu_xl<- as.data.frame(cbind(np_pola_eu_xl$layer[1:length(np_pola_eu_xl$layer)-1],
                                           np_pola_eu_xl[1:length(np_pola_eu_xl$freq)-1,3] - np_pola_eu_xl[2:length(np_pola_eu_xl$freq),3])) 
colnames(dif_areas_trop_eu_s) <- c("Time", "d_patches")
colnames(dif_areas_trop_eu_m) <- c("Time", "d_patches")
colnames(dif_areas_trop_eu_l) <- c("Time", "d_patches")
colnames(dif_areas_trop_eu_xl) <- c("Time", "d_patches")
colnames(dif_areas_arid_eu_s) <- c("Time", "d_patches")
colnames(dif_areas_arid_eu_m) <- c("Time", "d_patches")
colnames(dif_areas_arid_eu_l) <- c("Time", "d_patches")
colnames(dif_areas_arid_eu_xl) <- c("Time", "d_patches")
colnames(dif_areas_temp_eu_s) <- c("Time", "d_patches")
colnames(dif_areas_temp_eu_m) <- c("Time", "d_patches")
colnames(dif_areas_temp_eu_l) <- c("Time", "d_patches")
colnames(dif_areas_temp_eu_xl) <- c("Time", "d_patches")
colnames(dif_areas_cold_eu_s) <- c("Time", "d_patches")
colnames(dif_areas_cold_eu_m) <- c("Time", "d_patches")
colnames(dif_areas_cold_eu_l) <- c("Time", "d_patches")
colnames(dif_areas_cold_eu_xl) <- c("Time", "d_patches")
colnames(dif_areas_pola_eu_s) <- c("Time", "d_patches")
colnames(dif_areas_pola_eu_m) <- c("Time", "d_patches")
colnames(dif_areas_pola_eu_l) <- c("Time", "d_patches")
colnames(dif_areas_pola_eu_xl) <- c("Time", "d_patches")



AppendMe <- function(dfNames) {
  do.call(rbind, lapply(dfNames, function(x) {
    cbind(get(x), source = x)
  }))
}
dif_ww_s <- AppendMe(c("dif_areas_trop_ww_s","dif_areas_arid_ww_s","dif_areas_temp_ww_s",
                  "dif_areas_cold_ww_s","dif_areas_pola_ww_s"))
dif_ww_s[dif_ww_s$source=="dif_areas_trop_ww_s",3] <- "Tropical"
dif_ww_s[dif_ww_s$source=="dif_areas_arid_ww_s",3] <- "Arid"
dif_ww_s[dif_ww_s$source=="dif_areas_temp_ww_s",3] <- "Temperate"
dif_ww_s[dif_ww_s$source=="dif_areas_cold_ww_s",3] <- "Cold"
dif_ww_s[dif_ww_s$source=="dif_areas_pola_ww_s",3] <- "Polar"
dif_ww_s$source <- as.factor(dif_ww_s$source)
dif_ww_s$d_patches[dif_ww_s$d_patches<0] <- 0
head(dif_ww_s)

ggplot(dif_ww_s, aes(x=Time, y= d_patches, group=source, color=source)) + 
  geom_line(alpha=0.3) +
  geom_smooth(method="auto", se=TRUE, fullrange=FALSE, level=0.95) +
  scale_color_manual(values = c("Tropical" = "#f9d14a",
                                "Arid" = "#ab3329",
                                "Temperate" = "#ed968c",
                                "Cold" = "#7c4b73",
                                "Polar" = "#88a0dc" ))

dif_ww_m <- AppendMe(c("dif_areas_trop_ww_m","dif_areas_arid_ww_m","dif_areas_temp_ww_m",
                       "dif_areas_cold_ww_m","dif_areas_pola_ww_m"))
dif_ww_m[dif_ww_m$source=="dif_areas_trop_ww_m",3] <- "Tropical"
dif_ww_m[dif_ww_m$source=="dif_areas_arid_ww_m",3] <- "Arid"
dif_ww_m[dif_ww_m$source=="dif_areas_temp_ww_m",3] <- "Temperate"
dif_ww_m[dif_ww_m$source=="dif_areas_cold_ww_m",3] <- "Cold"
dif_ww_m[dif_ww_m$source=="dif_areas_pola_ww_m",3] <- "Polar"
dif_ww_m$source <- as.factor(dif_ww_m$source)
dif_ww_m$d_patches[dif_ww_m$d_patches<0] <- 0
head(dif_ww_m)

ggplot(dif_ww_m, aes(x=Time, y= d_patches, group=source, color=source)) + 
  geom_line(alpha=0.3) +
  geom_smooth(method="auto", se=TRUE, fullrange=FALSE, level=0.95) +
  scale_color_manual(values = c("Tropical" = "#f9d14a",
                                "Arid" = "#ab3329",
                                "Temperate" = "#ed968c",
                                "Cold" = "#7c4b73",
                                "Polar" = "#88a0dc" ))

dif_ww_l <- AppendMe(c("dif_areas_trop_ww_l","dif_areas_arid_ww_l","dif_areas_temp_ww_l",
                       "dif_areas_cold_ww_l","dif_areas_pola_ww_l"))
dif_ww_l[dif_ww_l$source=="dif_areas_trop_ww_l",3] <- "Tropical"
dif_ww_l[dif_ww_l$source=="dif_areas_arid_ww_l",3] <- "Arid"
dif_ww_l[dif_ww_l$source=="dif_areas_temp_ww_l",3] <- "Temperate"
dif_ww_l[dif_ww_l$source=="dif_areas_cold_ww_l",3] <- "Cold"
dif_ww_l[dif_ww_l$source=="dif_areas_pola_ww_l",3] <- "Polar"
dif_ww_l$source <- as.factor(dif_ww_l$source)
dif_ww_l$d_patches[dif_ww_l$d_patches<0] <- 0
head(dif_ww_l)

ggplot(dif_ww_l, aes(x=Time, y= d_patches, group=source, color=source)) + 
  geom_line(alpha=0.3) +
  geom_smooth(method="auto", se=TRUE, fullrange=FALSE, level=0.95) +
  scale_color_manual(values = c("Tropical" = "#f9d14a",
                                "Arid" = "#ab3329",
                                "Temperate" = "#ed968c",
                                "Cold" = "#7c4b73",
                                "Polar" = "#88a0dc" ))

dif_ww_xl <- AppendMe(c("dif_areas_trop_ww_xl","dif_areas_arid_ww_xl","dif_areas_temp_ww_xl",
                       "dif_areas_cold_ww_xl","dif_areas_pola_ww_xl"))
dif_ww_xl[dif_ww_xl$source=="dif_areas_trop_ww_xl",3] <- "Tropical"
dif_ww_xl[dif_ww_xl$source=="dif_areas_arid_ww_xl",3] <- "Arid"
dif_ww_xl[dif_ww_xl$source=="dif_areas_temp_ww_xl",3] <- "Temperate"
dif_ww_xl[dif_ww_xl$source=="dif_areas_cold_ww_xl",3] <- "Cold"
dif_ww_xl[dif_ww_xl$source=="dif_areas_pola_ww_xl",3] <- "Polar"
dif_ww_xl$source <- as.factor(dif_ww_xl$source)
dif_ww_xl$d_patches[dif_ww_xl$d_patches<0] <- 0
head(dif_ww_xl)

ggplot(dif_ww_xl, aes(x=Time, y= d_patches, group=source, color=source)) + 
  geom_line(alpha=0.01) +
  geom_smooth(method="auto", se=TRUE, fullrange=FALSE, level=0.95) +
  scale_color_manual(values = c("Tropical" = "#f9d14a",
                                "Arid" = "#ab3329",
                                "Temperate" = "#ed968c",
                                "Cold" = "#7c4b73",
                                "Polar" = "#88a0dc" ))

quartz()

p <- ggplot(dif_areas_cold_ww_m, aes(x=Time, y=d_patches)) +
  geom_line() + xlim(c(1,500))
plot(p)

#####Calculating fragmentation events
eventos_fragm_ww<- length (dif_areas_ww$d_patches[dif_areas_ww$d_patches>0])
eventos_fragm_am<- length (dif_areas_am$d_patches[dif_areas_am$d_patches>0])
eventos_fragm_af<- length (dif_areas_af$d_patches[dif_areas_af$d_patches>0])
eventos_fragm_eu<- length (dif_areas_eu$d_patches[dif_areas_eu$d_patches>0])

eventos_fragm_trop_ww_s<- length (dif_areas_trop_ww_s$d_patches[dif_areas_trop_ww_s$d_patches>0])
eventos_fragm_trop_ww_m<- length (dif_areas_trop_ww_m$d_patches[dif_areas_trop_ww_m$d_patches>0])
eventos_fragm_trop_ww_l<- length (dif_areas_trop_ww_l$d_patches[dif_areas_trop_ww_l$d_patches>0])
eventos_fragm_trop_ww_xl<- length (dif_areas_trop_ww_xl$d_patches[dif_areas_trop_ww_xl$d_patches>0])

eventos_fragm_arid_ww_s<- length (dif_areas_arid_ww_s$d_patches[dif_areas_arid_ww_s$d_patches>0])
eventos_fragm_arid_ww_m<- length (dif_areas_arid_ww_m$d_patches[dif_areas_arid_ww_m$d_patches>0])
eventos_fragm_arid_ww_l<- length (dif_areas_arid_ww_l$d_patches[dif_areas_arid_ww_l$d_patches>0])
eventos_fragm_arid_ww_xl<- length (dif_areas_arid_ww_xl$d_patches[dif_areas_arid_ww_xl$d_patches>0])

eventos_fragm_temp_ww_s<- length (dif_areas_temp_ww_s$d_patches[dif_areas_temp_ww_s$d_patches>0])
eventos_fragm_temp_ww_m<- length (dif_areas_temp_ww_m$d_patches[dif_areas_temp_ww_m$d_patches>0])
eventos_fragm_temp_ww_l<- length (dif_areas_temp_ww_l$d_patches[dif_areas_temp_ww_l$d_patches>0])
eventos_fragm_temp_ww_xl<- length (dif_areas_temp_ww_xl$d_patches[dif_areas_temp_ww_xl$d_patches>0])

eventos_fragm_cold_ww_s<- length (dif_areas_cold_ww_s$d_patches[dif_areas_cold_ww_s$d_patches>0])
eventos_fragm_cold_ww_m<- length (dif_areas_cold_ww_m$d_patches[dif_areas_cold_ww_m$d_patches>0])
eventos_fragm_cold_ww_l<- length (dif_areas_cold_ww_l$d_patches[dif_areas_cold_ww_l$d_patches>0])
eventos_fragm_cold_ww_xl<- length (dif_areas_cold_ww_xl$d_patches[dif_areas_cold_ww_xl$d_patches>0])

eventos_fragm_pola_ww_s<- length (dif_areas_pola_ww_s$d_patches[dif_areas_pola_ww_s$d_patches>0])
eventos_fragm_pola_ww_m<- length (dif_areas_pola_ww_m$d_patches[dif_areas_pola_ww_m$d_patches>0])
eventos_fragm_pola_ww_l<- length (dif_areas_pola_ww_l$d_patches[dif_areas_pola_ww_l$d_patches>0])
eventos_fragm_pola_ww_xl<- length (dif_areas_pola_ww_xl$d_patches[dif_areas_pola_ww_xl$d_patches>0])

eventos_ww_df <- as.data.frame(matrix(0,5,4))
colnames(eventos_ww_df) <- c("S","M","L","XL")
rownames(eventos_ww_df) <- c("Tropical", "Arid", "Temperate", "Cold", "Polar")
eventos_ww_df$S <- c(eventos_fragm_trop_ww_s, eventos_fragm_arid_ww_s, eventos_fragm_temp_ww_s,
                     eventos_fragm_cold_ww_s, eventos_fragm_pola_ww_s)
eventos_ww_df$M <- c(eventos_fragm_trop_ww_m, eventos_fragm_arid_ww_m, eventos_fragm_temp_ww_m,
                     eventos_fragm_cold_ww_m, eventos_fragm_pola_ww_m)
eventos_ww_df$L <- c(eventos_fragm_trop_ww_l, eventos_fragm_arid_ww_l, eventos_fragm_temp_ww_l,
                     eventos_fragm_cold_ww_l, eventos_fragm_pola_ww_l)
eventos_ww_df$XL <- c(eventos_fragm_trop_ww_xl, eventos_fragm_arid_ww_xl, eventos_fragm_temp_ww_xl,
                     eventos_fragm_cold_ww_xl, eventos_fragm_pola_ww_xl)

eventos_ww_df <- as.matrix(eventos_ww_df)
barplot(height = eventos_ww_df,beside = TRUE, col=c("#f9d14a", "#ab3329", "#ed968c", "#7c4b73", "#88a0dc"),
        border="white",# legend.text=c("Tropical", "Arid", "Temperate", "Cold", "Polar"),
        ylab= "Number of fragmentation events", ylim=c(0,5000))
quartz()

####America
eventos_fragm_trop_am_s<- length (dif_areas_trop_am_s$d_patches[dif_areas_trop_am_s$d_patches>0])
eventos_fragm_trop_am_m<- length (dif_areas_trop_am_m$d_patches[dif_areas_trop_am_m$d_patches>0])
eventos_fragm_trop_am_l<- length (dif_areas_trop_am_l$d_patches[dif_areas_trop_am_l$d_patches>0])
eventos_fragm_trop_am_xl<- length (dif_areas_trop_am_xl$d_patches[dif_areas_trop_am_xl$d_patches>0])

eventos_fragm_arid_am_s<- length (dif_areas_arid_am_s$d_patches[dif_areas_arid_am_s$d_patches>0])
eventos_fragm_arid_am_m<- length (dif_areas_arid_am_m$d_patches[dif_areas_arid_am_m$d_patches>0])
eventos_fragm_arid_am_l<- length (dif_areas_arid_am_l$d_patches[dif_areas_arid_am_l$d_patches>0])
eventos_fragm_arid_am_xl<- length (dif_areas_arid_am_xl$d_patches[dif_areas_arid_am_xl$d_patches>0])

eventos_fragm_temp_am_s<- length (dif_areas_temp_am_s$d_patches[dif_areas_temp_am_s$d_patches>0])
eventos_fragm_temp_am_m<- length (dif_areas_temp_am_m$d_patches[dif_areas_temp_am_m$d_patches>0])
eventos_fragm_temp_am_l<- length (dif_areas_temp_am_l$d_patches[dif_areas_temp_am_l$d_patches>0])
eventos_fragm_temp_am_xl<- length (dif_areas_temp_am_xl$d_patches[dif_areas_temp_am_xl$d_patches>0])

eventos_fragm_cold_am_s<- length (dif_areas_cold_am_s$d_patches[dif_areas_cold_am_s$d_patches>0])
eventos_fragm_cold_am_m<- length (dif_areas_cold_am_m$d_patches[dif_areas_cold_am_m$d_patches>0])
eventos_fragm_cold_am_l<- length (dif_areas_cold_am_l$d_patches[dif_areas_cold_am_l$d_patches>0])
eventos_fragm_cold_am_xl<- length (dif_areas_cold_am_xl$d_patches[dif_areas_cold_am_xl$d_patches>0])

eventos_fragm_pola_am_s<- length (dif_areas_pola_am_s$d_patches[dif_areas_pola_am_s$d_patches>0])
eventos_fragm_pola_am_m<- length (dif_areas_pola_am_m$d_patches[dif_areas_pola_am_m$d_patches>0])
eventos_fragm_pola_am_l<- length (dif_areas_pola_am_l$d_patches[dif_areas_pola_am_l$d_patches>0])
eventos_fragm_pola_am_xl<- length (dif_areas_pola_am_xl$d_patches[dif_areas_pola_am_xl$d_patches>0])

eventos_am_df <- as.data.frame(matrix(0,5,4))
colnames(eventos_am_df) <- c("S","M","L","XL")
rownames(eventos_am_df) <- c("Tropical", "Arid", "Temperate", "Cold", "Polar")
eventos_am_df$S <- c(eventos_fragm_trop_am_s, eventos_fragm_arid_am_s, eventos_fragm_temp_am_s,
                     eventos_fragm_cold_am_s, eventos_fragm_pola_am_s)
eventos_am_df$M <- c(eventos_fragm_trop_am_m, eventos_fragm_arid_am_m, eventos_fragm_temp_am_m,
                     eventos_fragm_cold_am_m, eventos_fragm_pola_am_m)
eventos_am_df$L <- c(eventos_fragm_trop_am_l, eventos_fragm_arid_am_l, eventos_fragm_temp_am_l,
                     eventos_fragm_cold_am_l, eventos_fragm_pola_am_l)
eventos_am_df$XL <- c(eventos_fragm_trop_am_xl, eventos_fragm_arid_am_xl, eventos_fragm_temp_am_xl,
                      eventos_fragm_cold_am_xl, eventos_fragm_pola_am_xl)

eventos_am_df <- as.matrix(eventos_am_df)
barplot(height = eventos_am_df,beside = TRUE, col=c("#f9d14a", "#ab3329", "#ed968c", "#7c4b73", "#88a0dc"),
        border="white",# legend.text=c("Tropical", "Arid", "Temperate", "Cold", "Polar"),
        ylab= "Number of fragmentation events", ylim=c(0,5000))
quartz()

###Africa
eventos_fragm_trop_af_s<- length (dif_areas_trop_af_s$d_patches[dif_areas_trop_af_s$d_patches>0])
eventos_fragm_trop_af_m<- length (dif_areas_trop_af_m$d_patches[dif_areas_trop_af_m$d_patches>0])
eventos_fragm_trop_af_l<- length (dif_areas_trop_af_l$d_patches[dif_areas_trop_af_l$d_patches>0])
eventos_fragm_trop_af_xl<- length (dif_areas_trop_af_xl$d_patches[dif_areas_trop_af_xl$d_patches>0])

eventos_fragm_arid_af_s<- length (dif_areas_arid_af_s$d_patches[dif_areas_arid_af_s$d_patches>0])
eventos_fragm_arid_af_m<- length (dif_areas_arid_af_m$d_patches[dif_areas_arid_af_m$d_patches>0])
eventos_fragm_arid_af_l<- length (dif_areas_arid_af_l$d_patches[dif_areas_arid_af_l$d_patches>0])
eventos_fragm_arid_af_xl<- length (dif_areas_arid_af_xl$d_patches[dif_areas_arid_af_xl$d_patches>0])

eventos_fragm_temp_af_s<- length (dif_areas_temp_af_s$d_patches[dif_areas_temp_af_s$d_patches>0])
eventos_fragm_temp_af_m<- length (dif_areas_temp_af_m$d_patches[dif_areas_temp_af_m$d_patches>0])
eventos_fragm_temp_af_l<- length (dif_areas_temp_af_l$d_patches[dif_areas_temp_af_l$d_patches>0])
eventos_fragm_temp_af_xl<- length (dif_areas_temp_af_xl$d_patches[dif_areas_temp_af_xl$d_patches>0])

eventos_fragm_cold_af_s<- length (dif_areas_cold_af_s$d_patches[dif_areas_cold_af_s$d_patches>0])
eventos_fragm_cold_af_m<- length (dif_areas_cold_af_m$d_patches[dif_areas_cold_af_m$d_patches>0])
eventos_fragm_cold_af_l<- length (dif_areas_cold_af_l$d_patches[dif_areas_cold_af_l$d_patches>0])
eventos_fragm_cold_af_xl<- length (dif_areas_cold_af_xl$d_patches[dif_areas_cold_af_xl$d_patches>0])

eventos_fragm_pola_af_s<- length (dif_areas_pola_af_s$d_patches[dif_areas_pola_af_s$d_patches>0])
eventos_fragm_pola_af_m<- length (dif_areas_pola_af_m$d_patches[dif_areas_pola_af_m$d_patches>0])
eventos_fragm_pola_af_l<- length (dif_areas_pola_af_l$d_patches[dif_areas_pola_af_l$d_patches>0])
eventos_fragm_pola_af_xl<- length (dif_areas_pola_af_xl$d_patches[dif_areas_pola_af_xl$d_patches>0])

eventos_af_df <- as.data.frame(matrix(0,5,4))
colnames(eventos_af_df) <- c("S","M","L","XL")
rownames(eventos_af_df) <- c("Tropical", "Arid", "Temperate", "Cold", "Polar")
eventos_af_df$S <- c(eventos_fragm_trop_af_s, eventos_fragm_arid_af_s, eventos_fragm_temp_af_s,
                     eventos_fragm_cold_af_s, eventos_fragm_pola_af_s)
eventos_af_df$M <- c(eventos_fragm_trop_af_m, eventos_fragm_arid_af_m, eventos_fragm_temp_af_m,
                     eventos_fragm_cold_af_m, eventos_fragm_pola_af_m)
eventos_af_df$L <- c(eventos_fragm_trop_af_l, eventos_fragm_arid_af_l, eventos_fragm_temp_af_l,
                     eventos_fragm_cold_af_l, eventos_fragm_pola_af_l)
eventos_af_df$XL <- c(eventos_fragm_trop_af_xl, eventos_fragm_arid_af_xl, eventos_fragm_temp_af_xl,
                      eventos_fragm_cold_af_xl, eventos_fragm_pola_af_xl)

eventos_af_df <- as.matrix(eventos_af_df)
barplot(height = eventos_af_df,beside = TRUE, col=c("#f9d14a", "#ab3329", "#ed968c", "#7c4b73", "#88a0dc"),
        border="white",# legend.text=c("Tropical", "Arid", "Temperate", "Cold", "Polar"),
        ylab= "Number of fragmentation events", ylim=c(0,5000))
quartz()

####Eurasia
eventos_fragm_trop_eu_s<- length (dif_areas_trop_eu_s$d_patches[dif_areas_trop_eu_s$d_patches>0])
eventos_fragm_trop_eu_m<- length (dif_areas_trop_eu_m$d_patches[dif_areas_trop_eu_m$d_patches>0])
eventos_fragm_trop_eu_l<- length (dif_areas_trop_eu_l$d_patches[dif_areas_trop_eu_l$d_patches>0])
eventos_fragm_trop_eu_xl<- length (dif_areas_trop_eu_xl$d_patches[dif_areas_trop_eu_xl$d_patches>0])

eventos_fragm_arid_eu_s<- length (dif_areas_arid_eu_s$d_patches[dif_areas_arid_eu_s$d_patches>0])
eventos_fragm_arid_eu_m<- length (dif_areas_arid_eu_m$d_patches[dif_areas_arid_eu_m$d_patches>0])
eventos_fragm_arid_eu_l<- length (dif_areas_arid_eu_l$d_patches[dif_areas_arid_eu_l$d_patches>0])
eventos_fragm_arid_eu_xl<- length (dif_areas_arid_eu_xl$d_patches[dif_areas_arid_eu_xl$d_patches>0])

eventos_fragm_temp_eu_s<- length (dif_areas_temp_eu_s$d_patches[dif_areas_temp_eu_s$d_patches>0])
eventos_fragm_temp_eu_m<- length (dif_areas_temp_eu_m$d_patches[dif_areas_temp_eu_m$d_patches>0])
eventos_fragm_temp_eu_l<- length (dif_areas_temp_eu_l$d_patches[dif_areas_temp_eu_l$d_patches>0])
eventos_fragm_temp_eu_xl<- length (dif_areas_temp_eu_xl$d_patches[dif_areas_temp_eu_xl$d_patches>0])

eventos_fragm_cold_eu_s<- length (dif_areas_cold_eu_s$d_patches[dif_areas_cold_eu_s$d_patches>0])
eventos_fragm_cold_eu_m<- length (dif_areas_cold_eu_m$d_patches[dif_areas_cold_eu_m$d_patches>0])
eventos_fragm_cold_eu_l<- length (dif_areas_cold_eu_l$d_patches[dif_areas_cold_eu_l$d_patches>0])
eventos_fragm_cold_eu_xl<- length (dif_areas_cold_eu_xl$d_patches[dif_areas_cold_eu_xl$d_patches>0])

eventos_fragm_pola_eu_s<- length (dif_areas_pola_eu_s$d_patches[dif_areas_pola_eu_s$d_patches>0])
eventos_fragm_pola_eu_m<- length (dif_areas_pola_eu_m$d_patches[dif_areas_pola_eu_m$d_patches>0])
eventos_fragm_pola_eu_l<- length (dif_areas_pola_eu_l$d_patches[dif_areas_pola_eu_l$d_patches>0])
eventos_fragm_pola_eu_xl<- length (dif_areas_pola_eu_xl$d_patches[dif_areas_pola_eu_xl$d_patches>0])

eventos_eu_df <- as.data.frame(matrix(0,5,4))
colnames(eventos_eu_df) <- c("S","M","L","XL")
rownames(eventos_eu_df) <- c("Tropical", "Arid", "Temperate", "Cold", "Polar")
eventos_eu_df$S <- c(eventos_fragm_trop_eu_s, eventos_fragm_arid_eu_s, eventos_fragm_temp_eu_s,
                     eventos_fragm_cold_eu_s, eventos_fragm_pola_eu_s)
eventos_eu_df$M <- c(eventos_fragm_trop_eu_m, eventos_fragm_arid_eu_m, eventos_fragm_temp_eu_m,
                     eventos_fragm_cold_eu_m, eventos_fragm_pola_eu_m)
eventos_eu_df$L <- c(eventos_fragm_trop_eu_l, eventos_fragm_arid_eu_l, eventos_fragm_temp_eu_l,
                     eventos_fragm_cold_eu_l, eventos_fragm_pola_eu_l)
eventos_eu_df$XL <- c(eventos_fragm_trop_eu_xl, eventos_fragm_arid_eu_xl, eventos_fragm_temp_eu_xl,
                      eventos_fragm_cold_eu_xl, eventos_fragm_pola_eu_xl)

eventos_eu_df <- as.matrix(eventos_eu_df)
barplot(height = eventos_eu_df,beside = TRUE, col=c("#f9d14a", "#ab3329", "#ed968c", "#7c4b73", "#88a0dc"),
        border="white",# legend.text=c("Tropical", "Arid", "Temperate", "Cold", "Polar"),
        ylab= "Number of fragmentation events", ylim=c(0,5000))
quartz()

######Calculating fragmentation strengh
strengh_fragm_ww<- median(dif_areas_ww$d_patches[dif_areas_ww$d_patches>0])
strengh_fragm_am<- median(dif_areas_am$d_patches[dif_areas_am$d_patches>0])
strengh_fragm_af<- median(dif_areas_af$d_patches[dif_areas_af$d_patches>0])
strengh_fragm_ww<- median(dif_areas_ww$d_patches[dif_areas_ww$d_patches>0])

###Worldwide
strengh_fragm_trop_ww_s<- median(dif_areas_trop_ww_s$d_patches[dif_areas_trop_ww_s$d_patches>0])
strengh_fragm_trop_ww_m<- median(dif_areas_trop_ww_m$d_patches[dif_areas_trop_ww_m$d_patches>0])
strengh_fragm_trop_ww_l<- median(dif_areas_trop_ww_l$d_patches[dif_areas_trop_ww_l$d_patches>0])
strengh_fragm_trop_ww_xl<- median(dif_areas_trop_ww_xl$d_patches[dif_areas_trop_ww_xl$d_patches>0])

strengh_fragm_arid_ww_s<- median(dif_areas_arid_ww_s$d_patches[dif_areas_arid_ww_s$d_patches>0])
strengh_fragm_arid_ww_m<- median(dif_areas_arid_ww_m$d_patches[dif_areas_arid_ww_m$d_patches>0])
strengh_fragm_arid_ww_l<- median(dif_areas_arid_ww_l$d_patches[dif_areas_arid_ww_l$d_patches>0])
strengh_fragm_arid_ww_xl<- median(dif_areas_arid_ww_xl$d_patches[dif_areas_arid_ww_xl$d_patches>0])

strengh_fragm_temp_ww_s<- median(dif_areas_temp_ww_s$d_patches[dif_areas_temp_ww_s$d_patches>0])
strengh_fragm_temp_ww_m<- median(dif_areas_temp_ww_m$d_patches[dif_areas_temp_ww_m$d_patches>0])
strengh_fragm_temp_ww_l<- median(dif_areas_temp_ww_l$d_patches[dif_areas_temp_ww_l$d_patches>0])
strengh_fragm_temp_ww_xl<- median(dif_areas_temp_ww_xl$d_patches[dif_areas_temp_ww_xl$d_patches>0])

strengh_fragm_cold_ww_s<- median(dif_areas_cold_ww_s$d_patches[dif_areas_cold_ww_s$d_patches>0])
strengh_fragm_cold_ww_m<- median(dif_areas_cold_ww_m$d_patches[dif_areas_cold_ww_m$d_patches>0])
strengh_fragm_cold_ww_l<- median(dif_areas_cold_ww_l$d_patches[dif_areas_cold_ww_l$d_patches>0])
strengh_fragm_cold_ww_xl<- median(dif_areas_cold_ww_xl$d_patches[dif_areas_cold_ww_xl$d_patches>0])

strengh_fragm_pola_ww_s<- median(dif_areas_pola_ww_s$d_patches[dif_areas_pola_ww_s$d_patches>0])
strengh_fragm_pola_ww_m<- median(dif_areas_pola_ww_m$d_patches[dif_areas_pola_ww_m$d_patches>0])
strengh_fragm_pola_ww_l<- median(dif_areas_pola_ww_l$d_patches[dif_areas_pola_ww_l$d_patches>0])
strengh_fragm_pola_ww_xl<- median(dif_areas_pola_ww_xl$d_patches[dif_areas_pola_ww_xl$d_patches>0])

strengh_ww_df <- as.data.frame(matrix(0,5,4))
colnames(strengh_ww_df) <- c("S","M","L","XL")
rownames(strengh_ww_df) <- c("Tropical", "Arid", "Temperate", "Cold", "Polar")
strengh_ww_df$S <- c(strengh_fragm_trop_ww_s, strengh_fragm_arid_ww_s, strengh_fragm_temp_ww_s,
                     strengh_fragm_cold_ww_s, strengh_fragm_pola_ww_s)
strengh_ww_df$M <- c(strengh_fragm_trop_ww_m, strengh_fragm_arid_ww_m, strengh_fragm_temp_ww_m,
                     strengh_fragm_cold_ww_m, strengh_fragm_pola_ww_m)
strengh_ww_df$L <- c(strengh_fragm_trop_ww_l, strengh_fragm_arid_ww_l, strengh_fragm_temp_ww_l,
                     strengh_fragm_cold_ww_l, strengh_fragm_pola_ww_l)
strengh_ww_df$XL <- c(strengh_fragm_trop_ww_xl, strengh_fragm_arid_ww_xl, strengh_fragm_temp_ww_xl,
                      strengh_fragm_cold_ww_xl, strengh_fragm_pola_ww_xl)

strengh_ww_df <- as.matrix(strengh_ww_df)
barplot(height = strengh_ww_df,beside = TRUE, col=c("#f9d14a", "#ab3329", "#ed968c", "#7c4b73", "#88a0dc"),
        border="white", legend.text=c("Tropical", "Arid", "Temperate", "Cold", "Polar"),
        ylab= "Median number of patches increase", ylim=c(0,8))
quartz()

###America
strengh_fragm_trop_am_s<- median(dif_areas_trop_am_s$d_patches[dif_areas_trop_am_s$d_patches>0])
strengh_fragm_trop_am_m<- median(dif_areas_trop_am_m$d_patches[dif_areas_trop_am_m$d_patches>0])
strengh_fragm_trop_am_l<- median(dif_areas_trop_am_l$d_patches[dif_areas_trop_am_l$d_patches>0])
strengh_fragm_trop_am_xl<- median(dif_areas_trop_am_xl$d_patches[dif_areas_trop_am_xl$d_patches>0])

strengh_fragm_arid_am_s<- median(dif_areas_arid_am_s$d_patches[dif_areas_arid_am_s$d_patches>0])
strengh_fragm_arid_am_m<- median(dif_areas_arid_am_m$d_patches[dif_areas_arid_am_m$d_patches>0])
strengh_fragm_arid_am_l<- median(dif_areas_arid_am_l$d_patches[dif_areas_arid_am_l$d_patches>0])
strengh_fragm_arid_am_xl<- median(dif_areas_arid_am_xl$d_patches[dif_areas_arid_am_xl$d_patches>0])

strengh_fragm_temp_am_s<- median(dif_areas_temp_am_s$d_patches[dif_areas_temp_am_s$d_patches>0])
strengh_fragm_temp_am_m<- median(dif_areas_temp_am_m$d_patches[dif_areas_temp_am_m$d_patches>0])
strengh_fragm_temp_am_l<- median(dif_areas_temp_am_l$d_patches[dif_areas_temp_am_l$d_patches>0])
strengh_fragm_temp_am_xl<- median(dif_areas_temp_am_xl$d_patches[dif_areas_temp_am_xl$d_patches>0])

strengh_fragm_cold_am_s<- median(dif_areas_cold_am_s$d_patches[dif_areas_cold_am_s$d_patches>0])
strengh_fragm_cold_am_m<- median(dif_areas_cold_am_m$d_patches[dif_areas_cold_am_m$d_patches>0])
strengh_fragm_cold_am_l<- median(dif_areas_cold_am_l$d_patches[dif_areas_cold_am_l$d_patches>0])
strengh_fragm_cold_am_xl<- median(dif_areas_cold_am_xl$d_patches[dif_areas_cold_am_xl$d_patches>0])

strengh_fragm_pola_am_s<- median(dif_areas_pola_am_s$d_patches[dif_areas_pola_am_s$d_patches>0])
strengh_fragm_pola_am_m<- median(dif_areas_pola_am_m$d_patches[dif_areas_pola_am_m$d_patches>0])
strengh_fragm_pola_am_l<- median(dif_areas_pola_am_l$d_patches[dif_areas_pola_am_l$d_patches>0])
strengh_fragm_pola_am_xl<- median(dif_areas_pola_am_xl$d_patches[dif_areas_pola_am_xl$d_patches>0])

strengh_am_df <- as.data.frame(matrix(0,5,4))
colnames(strengh_am_df) <- c("S","M","L","XL")
rownames(strengh_am_df) <- c("Tropical", "Arid", "Temperate", "Cold", "Polar")
strengh_am_df$S <- c(strengh_fragm_trop_am_s, strengh_fragm_arid_am_s, strengh_fragm_temp_am_s,
                     strengh_fragm_cold_am_s, strengh_fragm_pola_am_s)
strengh_am_df$M <- c(strengh_fragm_trop_am_m, strengh_fragm_arid_am_m, strengh_fragm_temp_am_m,
                     strengh_fragm_cold_am_m, strengh_fragm_pola_am_m)
strengh_am_df$L <- c(strengh_fragm_trop_am_l, strengh_fragm_arid_am_l, strengh_fragm_temp_am_l,
                     strengh_fragm_cold_am_l, strengh_fragm_pola_am_l)
strengh_am_df$XL <- c(strengh_fragm_trop_am_xl, strengh_fragm_arid_am_xl, strengh_fragm_temp_am_xl,
                      strengh_fragm_cold_am_xl, strengh_fragm_pola_am_xl)

strengh_am_df <- as.matrix(strengh_am_df)
barplot(height = strengh_am_df,beside = TRUE, col=c("#f9d14a", "#ab3329", "#ed968c", "#7c4b73", "#88a0dc"),
        border="white", legend.text=c("Tropical", "Arid", "Temperate", "Cold", "Polar"),
        ylab= "Median number of patches increase", ylim=c(0,8))
quartz()

###Africa
strengh_fragm_trop_af_s<- median(dif_areas_trop_af_s$d_patches[dif_areas_trop_af_s$d_patches>0])
strengh_fragm_trop_af_m<- median(dif_areas_trop_af_m$d_patches[dif_areas_trop_af_m$d_patches>0])
strengh_fragm_trop_af_l<- median(dif_areas_trop_af_l$d_patches[dif_areas_trop_af_l$d_patches>0])
strengh_fragm_trop_af_xl<- median(dif_areas_trop_af_xl$d_patches[dif_areas_trop_af_xl$d_patches>0])

strengh_fragm_arid_af_s<- median(dif_areas_arid_af_s$d_patches[dif_areas_arid_af_s$d_patches>0])
strengh_fragm_arid_af_m<- median(dif_areas_arid_af_m$d_patches[dif_areas_arid_af_m$d_patches>0])
strengh_fragm_arid_af_l<- median(dif_areas_arid_af_l$d_patches[dif_areas_arid_af_l$d_patches>0])
strengh_fragm_arid_af_xl<- median(dif_areas_arid_af_xl$d_patches[dif_areas_arid_af_xl$d_patches>0])

strengh_fragm_temp_af_s<- median(dif_areas_temp_af_s$d_patches[dif_areas_temp_af_s$d_patches>0])
strengh_fragm_temp_af_m<- median(dif_areas_temp_af_m$d_patches[dif_areas_temp_af_m$d_patches>0])
strengh_fragm_temp_af_l<- median(dif_areas_temp_af_l$d_patches[dif_areas_temp_af_l$d_patches>0])
strengh_fragm_temp_af_xl<- median(dif_areas_temp_af_xl$d_patches[dif_areas_temp_af_xl$d_patches>0])

strengh_fragm_cold_af_s<- median(dif_areas_cold_af_s$d_patches[dif_areas_cold_af_s$d_patches>0])
strengh_fragm_cold_af_m<- median(dif_areas_cold_af_m$d_patches[dif_areas_cold_af_m$d_patches>0])
strengh_fragm_cold_af_l<- median(dif_areas_cold_af_l$d_patches[dif_areas_cold_af_l$d_patches>0])
strengh_fragm_cold_af_xl<- median(dif_areas_cold_af_xl$d_patches[dif_areas_cold_af_xl$d_patches>0])

strengh_fragm_pola_af_s<- median(dif_areas_pola_af_s$d_patches[dif_areas_pola_af_s$d_patches>0])
strengh_fragm_pola_af_m<- median(dif_areas_pola_af_m$d_patches[dif_areas_pola_af_m$d_patches>0])
strengh_fragm_pola_af_l<- median(dif_areas_pola_af_l$d_patches[dif_areas_pola_af_l$d_patches>0])
strengh_fragm_pola_af_xl<- median(dif_areas_pola_af_xl$d_patches[dif_areas_pola_af_xl$d_patches>0])

strengh_af_df <- as.data.frame(matrix(0,5,4))
colnames(strengh_af_df) <- c("S","M","L","XL")
rownames(strengh_af_df) <- c("Tropical", "Arid", "Temperate", "Cold", "Polar")
strengh_af_df$S <- c(strengh_fragm_trop_af_s, strengh_fragm_arid_af_s, strengh_fragm_temp_af_s,
                     strengh_fragm_cold_af_s, strengh_fragm_pola_af_s)
strengh_af_df$M <- c(strengh_fragm_trop_af_m, strengh_fragm_arid_af_m, strengh_fragm_temp_af_m,
                     strengh_fragm_cold_af_m, strengh_fragm_pola_af_m)
strengh_af_df$L <- c(strengh_fragm_trop_af_l, strengh_fragm_arid_af_l, strengh_fragm_temp_af_l,
                     strengh_fragm_cold_af_l, strengh_fragm_pola_af_l)
strengh_af_df$XL <- c(strengh_fragm_trop_af_xl, strengh_fragm_arid_af_xl, strengh_fragm_temp_af_xl,
                      strengh_fragm_cold_af_xl, strengh_fragm_pola_af_xl)

strengh_af_df <- as.matrix(strengh_af_df)
barplot(height = strengh_af_df,beside = TRUE, col=c("#f9d14a", "#ab3329", "#ed968c", "#7c4b73", "#88a0dc"),
        border="white", legend.text=c("Tropical", "Arid", "Temperate", "Cold", "Polar"),
        ylab= "Median number of patches increase", ylim=c(0,8))
quartz()

###Europe
strengh_fragm_trop_eu_s<- median(dif_areas_trop_eu_s$d_patches[dif_areas_trop_eu_s$d_patches>0])
strengh_fragm_trop_eu_m<- median(dif_areas_trop_eu_m$d_patches[dif_areas_trop_eu_m$d_patches>0])
strengh_fragm_trop_eu_l<- median(dif_areas_trop_eu_l$d_patches[dif_areas_trop_eu_l$d_patches>0])
strengh_fragm_trop_eu_xl<- median(dif_areas_trop_eu_xl$d_patches[dif_areas_trop_eu_xl$d_patches>0])

strengh_fragm_arid_eu_s<- median(dif_areas_arid_eu_s$d_patches[dif_areas_arid_eu_s$d_patches>0])
strengh_fragm_arid_eu_m<- median(dif_areas_arid_eu_m$d_patches[dif_areas_arid_eu_m$d_patches>0])
strengh_fragm_arid_eu_l<- median(dif_areas_arid_eu_l$d_patches[dif_areas_arid_eu_l$d_patches>0])
strengh_fragm_arid_eu_xl<- median(dif_areas_arid_eu_xl$d_patches[dif_areas_arid_eu_xl$d_patches>0])

strengh_fragm_temp_eu_s<- median(dif_areas_temp_eu_s$d_patches[dif_areas_temp_eu_s$d_patches>0])
strengh_fragm_temp_eu_m<- median(dif_areas_temp_eu_m$d_patches[dif_areas_temp_eu_m$d_patches>0])
strengh_fragm_temp_eu_l<- median(dif_areas_temp_eu_l$d_patches[dif_areas_temp_eu_l$d_patches>0])
strengh_fragm_temp_eu_xl<- median(dif_areas_temp_eu_xl$d_patches[dif_areas_temp_eu_xl$d_patches>0])

strengh_fragm_cold_eu_s<- median(dif_areas_cold_eu_s$d_patches[dif_areas_cold_eu_s$d_patches>0])
strengh_fragm_cold_eu_m<- median(dif_areas_cold_eu_m$d_patches[dif_areas_cold_eu_m$d_patches>0])
strengh_fragm_cold_eu_l<- median(dif_areas_cold_eu_l$d_patches[dif_areas_cold_eu_l$d_patches>0])
strengh_fragm_cold_eu_xl<- median(dif_areas_cold_eu_xl$d_patches[dif_areas_cold_eu_xl$d_patches>0])

strengh_fragm_pola_eu_s<- median(dif_areas_pola_eu_s$d_patches[dif_areas_pola_eu_s$d_patches>0])
strengh_fragm_pola_eu_m<- median(dif_areas_pola_eu_m$d_patches[dif_areas_pola_eu_m$d_patches>0])
strengh_fragm_pola_eu_l<- median(dif_areas_pola_eu_l$d_patches[dif_areas_pola_eu_l$d_patches>0])
strengh_fragm_pola_eu_xl<- median(dif_areas_pola_eu_xl$d_patches[dif_areas_pola_eu_xl$d_patches>0])

strengh_eu_df <- as.data.frame(matrix(0,5,4))
colnames(strengh_eu_df) <- c("S","M","L","XL")
rownames(strengh_eu_df) <- c("Tropical", "Arid", "Temperate", "Cold", "Polar")
strengh_eu_df$S <- c(strengh_fragm_trop_eu_s, strengh_fragm_arid_eu_s, strengh_fragm_temp_eu_s,
                     strengh_fragm_cold_eu_s, strengh_fragm_pola_eu_s)
strengh_eu_df$M <- c(strengh_fragm_trop_eu_m, strengh_fragm_arid_eu_m, strengh_fragm_temp_eu_m,
                     strengh_fragm_cold_eu_m, strengh_fragm_pola_eu_m)
strengh_eu_df$L <- c(strengh_fragm_trop_eu_l, strengh_fragm_arid_eu_l, strengh_fragm_temp_eu_l,
                     strengh_fragm_cold_eu_l, strengh_fragm_pola_eu_l)
strengh_eu_df$XL <- c(strengh_fragm_trop_eu_xl, strengh_fragm_arid_eu_xl, strengh_fragm_temp_eu_xl,
                      strengh_fragm_cold_eu_xl, strengh_fragm_pola_eu_xl)

strengh_eu_df <- as.matrix(strengh_eu_df)
barplot(height = strengh_eu_df,beside = TRUE, col=c("#f9d14a", "#ab3329", "#ed968c", "#7c4b73", "#88a0dc"),
        border="white", legend.text=c("Tropical", "Arid", "Temperate", "Cold", "Polar"),
        ylab= "Median number of patches increase", ylim=c(0,8))
quartz()
######Calculating max fragmentation event
max_fragm_ww<- max(dif_areas_ww$d_patches[dif_areas_ww$d_patches>0])
max_fragm_am<- max(dif_areas_am$d_patches[dif_areas_am$d_patches>0])
max_fragm_af<- max(dif_areas_af$d_patches[dif_areas_af$d_patches>0])
max_fragm_ww<- max(dif_areas_ww$d_patches[dif_areas_ww$d_patches>0])

###Worldwide
max_fragm_trop_ww_s<- max(dif_areas_trop_ww_s$d_patches[dif_areas_trop_ww_s$d_patches>0])
max_fragm_trop_ww_m<- max(dif_areas_trop_ww_m$d_patches[dif_areas_trop_ww_m$d_patches>0])
max_fragm_trop_ww_l<- max(dif_areas_trop_ww_l$d_patches[dif_areas_trop_ww_l$d_patches>0])
max_fragm_trop_ww_xl<- max(dif_areas_trop_ww_xl$d_patches[dif_areas_trop_ww_xl$d_patches>0])

max_fragm_arid_ww_s<- max(dif_areas_arid_ww_s$d_patches[dif_areas_arid_ww_s$d_patches>0])
max_fragm_arid_ww_m<- max(dif_areas_arid_ww_m$d_patches[dif_areas_arid_ww_m$d_patches>0])
max_fragm_arid_ww_l<- max(dif_areas_arid_ww_l$d_patches[dif_areas_arid_ww_l$d_patches>0])
max_fragm_arid_ww_xl<- max(dif_areas_arid_ww_xl$d_patches[dif_areas_arid_ww_xl$d_patches>0])

max_fragm_temp_ww_s<- max(dif_areas_temp_ww_s$d_patches[dif_areas_temp_ww_s$d_patches>0])
max_fragm_temp_ww_m<- max(dif_areas_temp_ww_m$d_patches[dif_areas_temp_ww_m$d_patches>0])
max_fragm_temp_ww_l<- max(dif_areas_temp_ww_l$d_patches[dif_areas_temp_ww_l$d_patches>0])
max_fragm_temp_ww_xl<- max(dif_areas_temp_ww_xl$d_patches[dif_areas_temp_ww_xl$d_patches>0])

max_fragm_cold_ww_s<- max(dif_areas_cold_ww_s$d_patches[dif_areas_cold_ww_s$d_patches>0])
max_fragm_cold_ww_m<- max(dif_areas_cold_ww_m$d_patches[dif_areas_cold_ww_m$d_patches>0])
max_fragm_cold_ww_l<- max(dif_areas_cold_ww_l$d_patches[dif_areas_cold_ww_l$d_patches>0])
max_fragm_cold_ww_xl<- max(dif_areas_cold_ww_xl$d_patches[dif_areas_cold_ww_xl$d_patches>0])

max_fragm_pola_ww_s<- max(dif_areas_pola_ww_s$d_patches[dif_areas_pola_ww_s$d_patches>0])
max_fragm_pola_ww_m<- max(dif_areas_pola_ww_m$d_patches[dif_areas_pola_ww_m$d_patches>0])
max_fragm_pola_ww_l<- max(dif_areas_pola_ww_l$d_patches[dif_areas_pola_ww_l$d_patches>0])
max_fragm_pola_ww_xl<- max(dif_areas_pola_ww_xl$d_patches[dif_areas_pola_ww_xl$d_patches>0])

max_ww_df <- as.data.frame(matrix(0,5,4))
colnames(max_ww_df) <- c("S","M","L","XL")
rownames(max_ww_df) <- c("Tropical", "Arid", "Temperate", "Cold", "Polar")
max_ww_df$S <- c(max_fragm_trop_ww_s, max_fragm_arid_ww_s, max_fragm_temp_ww_s,
                     max_fragm_cold_ww_s, max_fragm_pola_ww_s)
max_ww_df$M <- c(max_fragm_trop_ww_m, max_fragm_arid_ww_m, max_fragm_temp_ww_m,
                     max_fragm_cold_ww_m, max_fragm_pola_ww_m)
max_ww_df$L <- c(max_fragm_trop_ww_l, max_fragm_arid_ww_l, max_fragm_temp_ww_l,
                     max_fragm_cold_ww_l, max_fragm_pola_ww_l)
max_ww_df$XL <- c(max_fragm_trop_ww_xl, max_fragm_arid_ww_xl, max_fragm_temp_ww_xl,
                      max_fragm_cold_ww_xl, max_fragm_pola_ww_xl)

max_ww_df <- as.matrix(max_ww_df)
barplot(height = max_ww_df,beside = TRUE, col=c("#f9d14a", "#ab3329", "#ed968c", "#7c4b73", "#88a0dc"),
        border="white", legend.text=c("Tropical", "Arid", "Temperate", "Cold", "Polar"),
        ylab= "Max number of patches increase", ylim=c(0,50))

quartz()
###America
max_fragm_trop_am_s<- max(dif_areas_trop_am_s$d_patches[dif_areas_trop_am_s$d_patches>0])
max_fragm_trop_am_m<- max(dif_areas_trop_am_m$d_patches[dif_areas_trop_am_m$d_patches>0])
max_fragm_trop_am_l<- max(dif_areas_trop_am_l$d_patches[dif_areas_trop_am_l$d_patches>0])
max_fragm_trop_am_xl<- max(dif_areas_trop_am_xl$d_patches[dif_areas_trop_am_xl$d_patches>0])

max_fragm_arid_am_s<- max(dif_areas_arid_am_s$d_patches[dif_areas_arid_am_s$d_patches>0])
max_fragm_arid_am_m<- max(dif_areas_arid_am_m$d_patches[dif_areas_arid_am_m$d_patches>0])
max_fragm_arid_am_l<- max(dif_areas_arid_am_l$d_patches[dif_areas_arid_am_l$d_patches>0])
max_fragm_arid_am_xl<- max(dif_areas_arid_am_xl$d_patches[dif_areas_arid_am_xl$d_patches>0])

max_fragm_temp_am_s<- max(dif_areas_temp_am_s$d_patches[dif_areas_temp_am_s$d_patches>0]) 
max_fragm_temp_am_m<- max(dif_areas_temp_am_m$d_patches[dif_areas_temp_am_m$d_patches>0])
max_fragm_temp_am_l<- max(dif_areas_temp_am_l$d_patches[dif_areas_temp_am_l$d_patches>0])
max_fragm_temp_am_xl<- max(dif_areas_temp_am_xl$d_patches[dif_areas_temp_am_xl$d_patches>0])

max_fragm_cold_am_s<- max(dif_areas_cold_am_s$d_patches[dif_areas_cold_am_s$d_patches>0])
max_fragm_cold_am_m<- max(dif_areas_cold_am_m$d_patches[dif_areas_cold_am_m$d_patches>0])
max_fragm_cold_am_l<- max(dif_areas_cold_am_l$d_patches[dif_areas_cold_am_l$d_patches>0])
max_fragm_cold_am_xl<- max(dif_areas_cold_am_xl$d_patches[dif_areas_cold_am_xl$d_patches>0])

max_fragm_pola_am_s<- max(dif_areas_pola_am_s$d_patches[dif_areas_pola_am_s$d_patches>0])
max_fragm_pola_am_m<- max(dif_areas_pola_am_m$d_patches[dif_areas_pola_am_m$d_patches>0])
max_fragm_pola_am_l<- max(dif_areas_pola_am_l$d_patches[dif_areas_pola_am_l$d_patches>0])
max_fragm_pola_am_xl<- max(dif_areas_pola_am_xl$d_patches[dif_areas_pola_am_xl$d_patches>0])

max_am_df <- as.data.frame(matrix(0,5,4))
colnames(max_am_df) <- c("S","M","L","XL")
rownames(max_am_df) <- c("Tropical", "Arid", "Temperate", "Cold", "Polar")
max_am_df$S <- c(max_fragm_trop_am_s, max_fragm_arid_am_s, max_fragm_temp_am_s,
                 max_fragm_cold_am_s, max_fragm_pola_am_s)
max_am_df$M <- c(max_fragm_trop_am_m, max_fragm_arid_am_m, max_fragm_temp_am_m,
                 max_fragm_cold_am_m, max_fragm_pola_am_m)
max_am_df$L <- c(max_fragm_trop_am_l, max_fragm_arid_am_l, max_fragm_temp_am_l,
                 max_fragm_cold_am_l, max_fragm_pola_am_l)
max_am_df$XL <- c(max_fragm_trop_am_xl, max_fragm_arid_am_xl, max_fragm_temp_am_xl,
                  max_fragm_cold_am_xl, max_fragm_pola_am_xl)

max_am_df <- as.matrix(max_am_df)
barplot(height = max_am_df,beside = TRUE, col=c("#f9d14a", "#ab3329", "#ed968c", "#7c4b73", "#88a0dc"),
        border="white", legend.text=c("Tropical", "Arid", "Temperate", "Cold", "Polar"),
        ylab= "Max number of patches increase", ylim=c(0,50))
quartz()

###Africa
max_fragm_trop_af_s<- max(dif_areas_trop_af_s$d_patches[dif_areas_trop_af_s$d_patches>0])
max_fragm_trop_af_m<- max(dif_areas_trop_af_m$d_patches[dif_areas_trop_af_m$d_patches>0])
max_fragm_trop_af_l<- max(dif_areas_trop_af_l$d_patches[dif_areas_trop_af_l$d_patches>0])
max_fragm_trop_af_xl<- max(dif_areas_trop_af_xl$d_patches[dif_areas_trop_af_xl$d_patches>0])

max_fragm_arid_af_s<- max(dif_areas_arid_af_s$d_patches[dif_areas_arid_af_s$d_patches>0])
max_fragm_arid_af_m<- max(dif_areas_arid_af_m$d_patches[dif_areas_arid_af_m$d_patches>0])
max_fragm_arid_af_l<- max(dif_areas_arid_af_l$d_patches[dif_areas_arid_af_l$d_patches>0])
max_fragm_arid_af_xl<- max(dif_areas_arid_af_xl$d_patches[dif_areas_arid_af_xl$d_patches>0])

max_fragm_temp_af_s<- max(dif_areas_temp_af_s$d_patches[dif_areas_temp_af_s$d_patches>0])
max_fragm_temp_af_m<- max(dif_areas_temp_af_m$d_patches[dif_areas_temp_af_m$d_patches>0])
max_fragm_temp_af_l<- max(dif_areas_temp_af_l$d_patches[dif_areas_temp_af_l$d_patches>0])
max_fragm_temp_af_xl<- max(dif_areas_temp_af_xl$d_patches[dif_areas_temp_af_xl$d_patches>0])

max_fragm_cold_af_s<- max(dif_areas_cold_af_s$d_patches[dif_areas_cold_af_s$d_patches>0])
max_fragm_cold_af_m<- max(dif_areas_cold_af_m$d_patches[dif_areas_cold_af_m$d_patches>0])
max_fragm_cold_af_l<- max(dif_areas_cold_af_l$d_patches[dif_areas_cold_af_l$d_patches>0])
max_fragm_cold_af_xl<- max(dif_areas_cold_af_xl$d_patches[dif_areas_cold_af_xl$d_patches>0])

max_fragm_pola_af_s<- max(dif_areas_pola_af_s$d_patches[dif_areas_pola_af_s$d_patches>0])
max_fragm_pola_af_m<- max(dif_areas_pola_af_m$d_patches[dif_areas_pola_af_m$d_patches>0])
max_fragm_pola_af_l<- max(dif_areas_pola_af_l$d_patches[dif_areas_pola_af_l$d_patches>0])
max_fragm_pola_af_xl<- max(dif_areas_pola_af_xl$d_patches[dif_areas_pola_af_xl$d_patches>0])

max_af_df <- as.data.frame(matrix(0,5,4))
colnames(max_af_df) <- c("S","M","L","XL")
rownames(max_af_df) <- c("Tropical", "Arid", "Temperate", "Cold", "Polar")
max_af_df$S <- c(max_fragm_trop_af_s, max_fragm_arid_af_s, max_fragm_temp_af_s,
                 max_fragm_cold_af_s, max_fragm_pola_af_s)
max_af_df$M <- c(max_fragm_trop_af_m, max_fragm_arid_af_m, max_fragm_temp_af_m,
                 max_fragm_cold_af_m, max_fragm_pola_af_m)
max_af_df$L <- c(max_fragm_trop_af_l, max_fragm_arid_af_l, max_fragm_temp_af_l,
                 max_fragm_cold_af_l, max_fragm_pola_af_l)
max_af_df$XL <- c(max_fragm_trop_af_xl, max_fragm_arid_af_xl, max_fragm_temp_af_xl,
                  max_fragm_cold_af_xl, max_fragm_pola_af_xl)
max_af_df[max_af_df==-Inf] <- 0
max_af_df <- as.matrix(max_af_df)
barplot(height = max_af_df,beside = TRUE, col=c("#f9d14a", "#ab3329", "#ed968c", "#7c4b73", "#88a0dc"),
        border="white", legend.text=c("Tropical", "Arid", "Temperate", "Cold", "Polar"),
        ylab= "Max number of patches increase", ylim=c(0,50))
quartz()

###europe
max_fragm_trop_eu_s<- max(dif_areas_trop_eu_s$d_patches[dif_areas_trop_eu_s$d_patches>0])
max_fragm_trop_eu_m<- max(dif_areas_trop_eu_m$d_patches[dif_areas_trop_eu_m$d_patches>0])
max_fragm_trop_eu_l<- max(dif_areas_trop_eu_l$d_patches[dif_areas_trop_eu_l$d_patches>0])
max_fragm_trop_eu_xl<- max(dif_areas_trop_eu_xl$d_patches[dif_areas_trop_eu_xl$d_patches>0])

max_fragm_arid_eu_s<- max(dif_areas_arid_eu_s$d_patches[dif_areas_arid_eu_s$d_patches>0])
max_fragm_arid_eu_m<- max(dif_areas_arid_eu_m$d_patches[dif_areas_arid_eu_m$d_patches>0])
max_fragm_arid_eu_l<- max(dif_areas_arid_eu_l$d_patches[dif_areas_arid_eu_l$d_patches>0])
max_fragm_arid_eu_xl<- max(dif_areas_arid_eu_xl$d_patches[dif_areas_arid_eu_xl$d_patches>0])

max_fragm_temp_eu_s<- max(dif_areas_temp_eu_s$d_patches[dif_areas_temp_eu_s$d_patches>0])
max_fragm_temp_eu_m<- max(dif_areas_temp_eu_m$d_patches[dif_areas_temp_eu_m$d_patches>0])
max_fragm_temp_eu_l<- max(dif_areas_temp_eu_l$d_patches[dif_areas_temp_eu_l$d_patches>0])
max_fragm_temp_eu_xl<- max(dif_areas_temp_eu_xl$d_patches[dif_areas_temp_eu_xl$d_patches>0])

max_fragm_cold_eu_s<- max(dif_areas_cold_eu_s$d_patches[dif_areas_cold_eu_s$d_patches>0])
max_fragm_cold_eu_m<- max(dif_areas_cold_eu_m$d_patches[dif_areas_cold_eu_m$d_patches>0])
max_fragm_cold_eu_l<- max(dif_areas_cold_eu_l$d_patches[dif_areas_cold_eu_l$d_patches>0])
max_fragm_cold_eu_xl<- max(dif_areas_cold_eu_xl$d_patches[dif_areas_cold_eu_xl$d_patches>0])

max_fragm_pola_eu_s<- max(dif_areas_pola_eu_s$d_patches[dif_areas_pola_eu_s$d_patches>0])
max_fragm_pola_eu_m<- max(dif_areas_pola_eu_m$d_patches[dif_areas_pola_eu_m$d_patches>0])
max_fragm_pola_eu_l<- max(dif_areas_pola_eu_l$d_patches[dif_areas_pola_eu_l$d_patches>0])
max_fragm_pola_eu_xl<- max(dif_areas_pola_eu_xl$d_patches[dif_areas_pola_eu_xl$d_patches>0])

max_eu_df <- as.data.frame(matrix(0,5,4))
colnames(max_eu_df) <- c("S","M","L","XL")
rownames(max_eu_df) <- c("Tropical", "Arid", "Temperate", "Cold", "Polar")
max_eu_df$S <- c(max_fragm_trop_eu_s, max_fragm_arid_eu_s, max_fragm_temp_eu_s,
                 max_fragm_cold_eu_s, max_fragm_pola_eu_s)
max_eu_df$M <- c(max_fragm_trop_eu_m, max_fragm_arid_eu_m, max_fragm_temp_eu_m,
                 max_fragm_cold_eu_m, max_fragm_pola_eu_m)
max_eu_df$L <- c(max_fragm_trop_eu_l, max_fragm_arid_eu_l, max_fragm_temp_eu_l,
                 max_fragm_cold_eu_l, max_fragm_pola_eu_l)
max_eu_df$XL <- c(max_fragm_trop_eu_xl, max_fragm_arid_eu_xl, max_fragm_temp_eu_xl,
                  max_fragm_cold_eu_xl, max_fragm_pola_eu_xl)
max_eu_df <- as.matrix(max_eu_df)

max_eu_df <- as.matrix(max_eu_df)
barplot(height = max_eu_df,beside = TRUE, col=c("#f9d14a", "#ab3329", "#ed968c", "#7c4b73", "#88a0dc"),
        border="white", legend.text=c("Tropical", "Arid", "Temperate", "Cold", "Polar"),
        ylab= "Max number of patches increase", ylim=c(0,50))
quartz()

p <- ggplot(max_eu_df,                                      # Grouped barplot using ggplot2
       aes(x = colnames(max_eu_df),
           y = values,
           fill = subgroup)) +
  geom_bar(stat = "identity",
           position = "dodge")

###############################
######Maps by climate zone########
###############################
mapa <- stack_all[[1]]
mapa_tropical <- mapa
mapa_tropical[mapa_tropical!=1] <- NA
mapa_arido <- mapa
mapa_arido[mapa_arido!=2] <- NA
mapa_arido[mapa_arido==2] <- 1
mapa_templado <- mapa
mapa_templado[mapa_templado!=3] <- NA
mapa_templado[mapa_templado==3] <- 1
mapa_frio <- mapa
mapa_frio[mapa_frio!=4] <- NA
mapa_frio[mapa_frio==4] <- 1
mapa_polar <- mapa
mapa_polar[mapa_polar!=5] <- NA
mapa_polar[mapa_polar==5] <- 1
plot(mapa_polar)


mapa_E <- Eurasia_stack[[1]]
mapa_tropical_E <- mapa_E
mapa_tropical_E[mapa_tropical_E!=1] <- NA
mapa_arido_E <- mapa_E
mapa_arido_E[mapa_arido_E!=2] <- NA
mapa_arido_E[mapa_arido_E==2] <- 1
mapa_templado_E <- mapa_E
mapa_templado_E[mapa_templado_E!=3] <- NA
mapa_templado_E[mapa_templado_E==3] <- 1
mapa_frio_E <- mapa_E
mapa_frio_E[mapa_frio_E!=4] <- NA
mapa_frio_E[mapa_frio_E==4] <- 1
mapa_polar_E<- mapa_E
mapa_polar_E[mapa_polar_E!=5] <- NA
mapa_polar_E[mapa_polar_E==5] <- 1
plot(mapa_polar_E)


mapa_Am <- America_stack[[1]]
mapa_tropical_Am <- mapa_Am
mapa_tropical_Am[mapa_tropical_Am!=1] <- NA
mapa_arido_Am <- mapa_Am
mapa_arido_Am[mapa_arido_Am!=2] <- NA
mapa_arido_Am[mapa_arido_Am==2] <- 1
mapa_templado_Am <- mapa_Am
mapa_templado_Am[mapa_templado_Am!=3] <- NA
mapa_templado_Am[mapa_templado_Am==3] <- 1
mapa_frio_Am <- mapa_Am
mapa_frio_Am[mapa_frio_Am!=4] <- NA
mapa_frio_Am[mapa_frio_Am==4] <- 1
mapa_polar_Am<- mapa_Am
mapa_polar_Am[mapa_polar_Am!=5] <- NA
mapa_polar_Am[mapa_polar_Am==5] <- 1
plot(mapa_polar_Am)


mapa_Af <- Africa_stack[[1]]
mapa_tropical_Af <- mapa_Af
mapa_tropical_Af[mapa_tropical_Af!=1] <- NA
mapa_arido_Af <- mapa_Af
mapa_arido_Af[mapa_arido_Af!=2] <- NA
mapa_arido_Af[mapa_arido_Af==2] <- 1
mapa_templado_Af <- mapa_Af
mapa_templado_Af[mapa_templado_Af!=3] <- NA
mapa_templado_Af[mapa_templado_Af==3] <- 1
mapa_frio_Af <- mapa_Af
mapa_frio_Af[mapa_frio_Af!=4] <- NA
mapa_frio_Af[mapa_frio_Af==4] <- 1
mapa_polar_Af<- mapa_Af
mapa_polar_Af[mapa_polar_Af!=5] <- NA
mapa_polar_Af[mapa_polar_Af==5] <- 1
plot(mapa_polar_Af)
plot(mapa_tropical_E)


###########################################
#########Extract richnes by climate zone########
############################3#############
# rich_tropical <- Richness_moll*mapa_tropical
# rich_arido <- Richness_moll*mapa_arido
# rich_templado <- Richness_moll*mapa_templado
# rich_frio <- Richness_moll*mapa_frio
# rich_polar <- Richness_moll*mapa_polar
# 
# rich_tropical_Am <- Richness_moll*mapa_tropical_Am
# rich_arido_Am <- Richness_moll*mapa_arido_Am
# rich_templado_Am <- Richness_moll*mapa_templado_Am
# rich_frio_Am <- Richness_moll*mapa_frio_Am
# rich_polar_Am <- Richness_moll*mapa_polar_Am
# 
# rich_tropical_Af <- Richness_moll*mapa_tropical_Af
# rich_arido_Af <- Richness_moll*mapa_arido_Af
# rich_templado_Af <- Richness_moll*mapa_templado_Af
# rich_frio_Af <- Richness_moll*mapa_frio_Af
# rich_polar_Af <- Richness_moll*mapa_polar_Af
# 
# rich_tropical_E <- Richness_moll*mapa_tropical_E
# rich_arido_E <- Richness_moll*mapa_arido_E
# rich_templado_E <- Richness_moll*mapa_templado_E
# rich_frio_E <- Richness_moll*mapa_frio_E
# rich_polar_E <- Richness_moll*mapa_polar_E
# 
# stack_rich <- stack(rich_tropical, rich_arido,rich_templado,rich_frio,rich_polar,rich_tropical_Am,
#                     rich_arido_Am,rich_templado_Am,rich_frio_Am,rich_polar_Am,rich_tropical_Af,
#                     rich_arido_Af,rich_templado_Af,rich_frio_Af,rich_polar_Af,rich_tropical_E,
#                     rich_arido_E,rich_templado_E,rich_frio_E,rich_polar_E)

#####Total richness by map
cur_clim <- stack_all[[1]]
cur_clim <- projectRaster(cur_clim, crs=" +init=epsg:4326 +proj=longlat
+ellps=WGS84 +datum=WGS84 +no_defs +towgs84=0,0,0")
clim_class <- as.data.frame(cur_clim,xy=T)
 clim_class$x <- round(clim_class$x,2)
 clim_class$x <- clim_class$x-0.24
 clim_class$y <- round(clim_class$y,2)
 clim_class$y <- clim_class$y +0.08
 clim_class$Stack_reclassified_moll_1to500_1 <- round(clim_class$Stack_reclassified_moll_1to500_1,0)
which(mammals_db$Y==84.75,)
which(mammals_db$Y==-83.75,)
mammals_db$Clim_zone[7201:250560] <- clim_class$Stack_reclassified_moll_1to500_1

#### 
unique(mask_Af$Africa_stack_1)
mask_Af <-mapa_Af 
mask_Af <- projectRaster(mask_Af, crs=" +init=epsg:4326 +proj=longlat
+ellps=WGS84 +datum=WGS84 +no_defs +towgs84=0,0,0")
mask_Af <- as.data.frame(mask_Af,xy=T)
mask_Af$x <- round(mask_Af$x,2)
mask_Af$x <- mask_Af$x-0.24
mask_Af$y <- round(mask_Af$y,2)
mask_Af$y <- mask_Af$y +0.08
mask_Af[is.na(mask_Af)==F] <- 2
mask_Af[is.na(mask_Af)==T] <- 0
mammals_db$Continent <- 0
mammals_db$Continent[7201:250560] <- mask_Af$Africa_stack_1

mask_Am <-mapa_Am
mask_Am <- projectRaster(mask_Am, crs=" +init=epsg:4326 +proj=longlat
+ellps=WGS84 +datum=WGS84 +no_defs +towgs84=0,0,0")
mask_Am <- as.data.frame(mask_Am,xy=T)
mask_Am$x <- round(mask_Am$x,2)
mask_Am$x <- mask_Am$x-0.24
mask_Am$y <- round(mask_Am$y,2)
mask_Am$y <- mask_Am$y +0.08
mask_Am[is.na(mask_Am)==F] <- 1
mask_Am[is.na(mask_Am)==T] <- 0
mammals_db$Continent[7201:250560] <- mammals_db$Continent[7201:250560] + mask_Am$America_stack_1

mask_Eu <-mapa_E
mask_Eu <- projectRaster(mask_Eu, crs=" +init=epsg:4326 +proj=longlat
+ellps=WGS84 +datum=WGS84 +no_defs +towgs84=0,0,0")
mask_Eu <- as.data.frame(mask_Eu,xy=T)
mask_Eu$x <- round(mask_Eu$x,2)
mask_Eu$x <- mask_Eu$x-0.24
mask_Eu$y <- round(mask_Eu$y,2)
mask_Eu$y <- mask_Eu$y +0.08
mask_Eu[is.na(mask_Eu)==F] <- 3
mask_Eu[is.na(mask_Eu)==T] <- 0
mammals_db$Continent[7201:250560] <- mammals_db$Continent[7201:250560] + mask_Eu$Eurasia_stack_1
mammals_db$Continent[mammals_db$Continent==4] <- 1
mammals_db$Continent[mammals_db$Continent==5 ] <- 2
unique(mammals_db$Continent)
unique(mammals_db$Clim_zone)

MM_rich <- as.data.frame(matrix(0,nrow=5,ncol=4))
colnames(MM_rich) <- c("Total", "America", "Africa", "Eurasia")
rownames(MM_rich) <- c("Tropical", "Arid", "Temperate", "Cold", "Polar")

tail(names(specialist_db))
specialist_db <- mammals_db
specialist_db[,3:500] <- specialist_db[,3:500]*specialist_db$Clim_zone
specialist_db[,501:1001] <- specialist_db[,501:1001]*specialist_db$Clim_zone
specialist_db[,1002:1500] <- specialist_db[,1002:1500]*specialist_db$Clim_zone
specialist_db[,1501:2000] <- specialist_db[,1501:2000]*specialist_db$Clim_zone
specialist_db[,2001:2500] <- specialist_db[,2001:2500]*specialist_db$Clim_zone
specialist_db[,2501:3000] <- specialist_db[,2501:3000]*specialist_db$Clim_zone
specialist_db[,3001:3500] <- specialist_db[,3001:3500]*specialist_db$Clim_zone
specialist_db[,3501:4000] <- specialist_db[,3501:4000]*specialist_db$Clim_zone
specialist_db[,4001:4500] <- specialist_db[,4001:4500]*specialist_db$Clim_zone
specialist_db[,4501:5000] <- specialist_db[,4501:5000]*specialist_db$Clim_zone
specialist_db[,5001:5500] <- specialist_db[,5001:5500]*specialist_db$Clim_zone
specialist_db[,5501:5743] <- specialist_db[,5501:5741]*specialist_db$Clim_zone

dim(mammals_db)
dim(specialist_db)
specialist_db$Land_mass <- mammals_db$Land_mass
unique(specialist_db$Land_mass)

library(dplyr)
sub_espe1 <- specialist_db[3:10] %>%
  select_if(~ length(unique(na.omit(.))) == 1)
sub_espe2 <- specialist_db[1001:2000] %>%
  select_if(~ length(unique(na.omit(.))) == 1)
sub_espe3 <- specialist_db[2001:3000] %>%
  select_if(~ length(unique(na.omit(.))) == 1)
sub_espe4 <- specialist_db[3001:4000] %>%
  select_if(~ length(unique(na.omit(.))) == 1)
sub_espe5 <- specialist_db[4001:5000] %>%
  select_if(~ length(unique(na.omit(.))) == 1)
sub_espe6 <- specialist_db[5001:5741] %>%
  select_if(~ length(unique(na.omit(.))) == 1)

sub_espe <- cbind(sub_espe1, sub_espe2, sub_espe3, sub_espe4, sub_espe5, sub_espe6)
rm(sub_espe1, sub_espe2, sub_espe3, sub_espe4, sub_espe5, sub_espe6)
dim(sub_espe)

unique_values <- map(sub_espe, unique)
  for (i in 1:1417) {
    ele <- unique_values[[i]]
    ele <- ele[!sapply(ele, is.nan)]
    ele <- ele[!sapply(ele, is.na)]
    unique_values[[i]] <- ele
  }
unique_vector <- unlist(unique_values)
table(unique_vector)

specialist_Am <- mammals_db[mammals_db$Land_mass==1,]
dim(specialist_Am)
specialist_Am[,3:500] <- specialist_Am[,3:500]*specialist_Am$Clim_zone
specialist_Am[,501:1001] <- specialist_Am[,501:1001]*specialist_Am$Clim_zone
specialist_Am[,1002:1500] <- specialist_Am[,1002:1500]*specialist_Am$Clim_zone
specialist_Am[,1501:2000] <- specialist_Am[,1501:2000]*specialist_Am$Clim_zone
specialist_Am[,2001:2500] <- specialist_Am[,2001:2500]*specialist_Am$Clim_zone
specialist_Am[,2501:3000] <- specialist_Am[,2501:3000]*specialist_Am$Clim_zone
specialist_Am[,3001:3500] <- specialist_Am[,3001:3500]*specialist_Am$Clim_zone
specialist_Am[,3501:4000] <- specialist_Am[,3501:4000]*specialist_Am$Clim_zone
specialist_Am[,4001:4500] <- specialist_Am[,4001:4500]*specialist_Am$Clim_zone
specialist_Am[,4501:5000] <- specialist_Am[,4501:5000]*specialist_Am$Clim_zone
specialist_Am[,5001:5500] <- specialist_Am[,5001:5500]*specialist_Am$Clim_zone
specialist_Am[,5501:5743] <- specialist_Am[,5501:5741]*specialist_Am$Clim_zone

sub_Am1 <- specialist_Am[,3:1000] %>%
  select_if(~ length(unique(na.omit(.))) == 1)
sub_Am2 <- specialist_Am[,1001:2000] %>%
  select_if(~ length(unique(na.omit(.))) == 1)
sub_Am3 <- specialist_Am[,2001:3000] %>%
  select_if(~ length(unique(na.omit(.))) == 1)
sub_Am4 <- specialist_Am[,3001:4000] %>%
  select_if(~ length(unique(na.omit(.))) == 1)
sub_Am5 <- specialist_Am[,4001:5000] %>%
  select_if(~ length(unique(na.omit(.))) == 1)
sub_Am6 <- specialist_Am[,5001:5741] %>%
  select_if(~ length(unique(na.omit(.))) == 1)

sub_Am <- cbind(sub_Am1, sub_Am2, sub_Am3, sub_Am4, sub_Am5, sub_Am6)
rm(sub_Am1, sub_Am2, sub_Am3, sub_Am4, sub_Am5, sub_Am6)
dim(sub_Am)
rm(specialist_db)

unique_Am <- map(sub_Am, unique)
for (i in 1:455) {
  ele <- unique_Am[[i]]
  ele <- ele[!sapply(ele, is.nan)]
  ele <- ele[!sapply(ele, is.na)]
  unique_Am[[i]] <- ele
}
unique_vector_Am <- unlist(unique_Am)
table(unique_vector_Am)
rm(patch_metrics_full_names)

specialist_Af <- mammals_db[mammals_db$Land_mass==2,]
dim(specialist_Af)
specialist_Af[,3:500] <- specialist_Af[,3:500]*specialist_Af$Clim_zone
specialist_Af[,501:1001] <- specialist_Af[,501:1001]*specialist_Af$Clim_zone
specialist_Af[,1002:1500] <- specialist_Af[,1002:1500]*specialist_Af$Clim_zone
specialist_Af[,1501:2000] <- specialist_Af[,1501:2000]*specialist_Af$Clim_zone
specialist_Af[,2001:2500] <- specialist_Af[,2001:2500]*specialist_Af$Clim_zone
specialist_Af[,2501:3000] <- specialist_Af[,2501:3000]*specialist_Af$Clim_zone
specialist_Af[,3001:3500] <- specialist_Af[,3001:3500]*specialist_Af$Clim_zone
specialist_Af[,3501:4000] <- specialist_Af[,3501:4000]*specialist_Af$Clim_zone
specialist_Af[,4001:4500] <- specialist_Af[,4001:4500]*specialist_Af$Clim_zone
specialist_Af[,4501:5000] <- specialist_Af[,4501:5000]*specialist_Af$Clim_zone
specialist_Af[,5001:5500] <- specialist_Af[,5001:5500]*specialist_Af$Clim_zone
specialist_Af[,5501:5743] <- specialist_Af[,5501:5741]*specialist_Af$Clim_zone

sub_Af1 <- specialist_Af[,3:1000] %>%
  select_if(~ length(unique(na.omit(.))) == 1)
sub_Af2 <- specialist_Af[,1001:2000] %>%
  select_if(~ length(unique(na.omit(.))) == 1)
sub_Af3 <- specialist_Af[,2001:3000] %>%
  select_if(~ length(unique(na.omit(.))) == 1)
sub_Af4 <- specialist_Af[,3001:4000] %>%
  select_if(~ length(unique(na.omit(.))) == 1)
sub_Af5 <- specialist_Af[,4001:5000] %>%
  select_if(~ length(unique(na.omit(.))) == 1)
sub_Af6 <- specialist_Af[,5001:5741] %>%
  select_if(~ length(unique(na.omit(.))) == 1)

sub_Af <- cbind(sub_Af1, sub_Af2, sub_Af3, sub_Af4, sub_Af5, sub_Af6)
rm(sub_Af1, sub_Af2, sub_Af3, sub_Af4, sub_Af5, sub_Af6)
dim(sub_Af)

unique_Af <- map(sub_Af, unique)
for (i in 1:423) {
  ele <- unique_Af[[i]]
  ele <- ele[!sapply(ele, is.nan)]
  ele <- ele[!sapply(ele, is.na)]
  unique_Af[[i]] <- ele
}
unique_vector_Af <- unlist(unique_Af)
table(unique_vector_Af)


specialist_Eu <- mammals_db[mammals_db$Land_mass==3,]
dim(specialist_Eu)
specialist_Eu[,3:500] <- specialist_Eu[,3:500]*specialist_Eu$Clim_zone
specialist_Eu[,501:1001] <- specialist_Eu[,501:1001]*specialist_Eu$Clim_zone
specialist_Eu[,1002:1500] <- specialist_Eu[,1002:1500]*specialist_Eu$Clim_zone
specialist_Eu[,1501:2000] <- specialist_Eu[,1501:2000]*specialist_Eu$Clim_zone
specialist_Eu[,2001:2500] <- specialist_Eu[,2001:2500]*specialist_Eu$Clim_zone
specialist_Eu[,2501:3000] <- specialist_Eu[,2501:3000]*specialist_Eu$Clim_zone
specialist_Eu[,3001:3500] <- specialist_Eu[,3001:3500]*specialist_Eu$Clim_zone
specialist_Eu[,3501:4000] <- specialist_Eu[,3501:4000]*specialist_Eu$Clim_zone
specialist_Eu[,4001:4500] <- specialist_Eu[,4001:4500]*specialist_Eu$Clim_zone
specialist_Eu[,4501:5000] <- specialist_Eu[,4501:5000]*specialist_Eu$Clim_zone
specialist_Eu[,5001:5500] <- specialist_Eu[,5001:5500]*specialist_Eu$Clim_zone
specialist_Eu[,5501:5743] <- specialist_Eu[,5501:5741]*specialist_Eu$Clim_zone

sub_Eu1 <- specialist_Eu[,3:1000] %>%
  select_if(~ length(unique(na.omit(.))) == 1)
sub_Eu2 <- specialist_Eu[,1001:2000] %>%
  select_if(~ length(unique(na.omit(.))) == 1)
sub_Eu3 <- specialist_Eu[,2001:3000] %>%
  select_if(~ length(unique(na.omit(.))) == 1)
sub_Eu4 <- specialist_Eu[,3001:4000] %>%
  select_if(~ length(unique(na.omit(.))) == 1)
sub_Eu5 <- specialist_Eu[,4001:5000] %>%
  select_if(~ length(unique(na.omit(.))) == 1)
sub_Eu6 <- specialist_Eu[,5001:5741] %>%
  select_if(~ length(unique(na.omit(.))) == 1)

sub_Eu <- cbind(sub_Eu1, sub_Eu2, sub_Eu3, sub_Eu4, sub_Eu5, sub_Eu6)
rm(sub_Eu1, sub_Eu2, sub_Eu3, sub_Eu4, sub_Eu5, sub_Eu6)
dim(sub_Eu)
rm(unique_vector_Eu)

unique_Eu <- map(sub_Eu, unique)
for (i in 1:892) {
  ele <- unique_Eu[[i]]
  ele <- ele[!sapply(ele, is.nan)]
  ele <- ele[!sapply(ele, is.na)]
  unique_Eu[[i]] <- ele
}
unique_vector_Eu <- unlist(unique_Eu)
table(unique_vector_Eu)






mammals_trop <- mammals_db[mammals_db$Clim_zone==1,]
dim(mammals_trop)
tail(names(mammals_trop))
mammals_trop_Total <- rep(0, 5741)
mammals_trop_Total[3:500] <- colSums(mammals_trop[3:500], na.rm=T)
mammals_trop_Total[501:1500] <- colSums(mammals_trop[501:1500], na.rm=T)
mammals_trop_Total[1501:2500] <- colSums(mammals_trop[1501:2500], na.rm=T)
mammals_trop_Total[2501:3500] <- colSums(mammals_trop[2501:3500], na.rm=T)
mammals_trop_Total[3501:4500] <- colSums(mammals_trop[3501:4500], na.rm=T)
mammals_trop_Total[4501:5500] <- colSums(mammals_trop[4501:5500], na.rm=T)
mammals_trop_Total[5501:5741] <- colSums(mammals_trop[5501:5741], na.rm=T)
mammals_trop_Total[mammals_trop_Total>0] <- 1
MM_rich[1,1] <- sum(mammals_trop_Total)

mammals_trop_am <- mammals_trop[mammals_trop$Continent==1,]
length(colSums(mammals_trop_am[4501:5741], na.rm=T))
mammals_trop_am_Total <- rep(0, 5741)
mammals_trop_am_Total[3:500] <- colSums(mammals_trop_am[3:500], na.rm=T)
mammals_trop_am_Total[501:1500] <- colSums(mammals_trop_am[501:1500], na.rm=T)
mammals_trop_am_Total[1501:2500] <- colSums(mammals_trop_am[1501:2500], na.rm=T)
mammals_trop_am_Total[2501:3500] <- colSums(mammals_trop_am[2501:3500], na.rm=T)
mammals_trop_am_Total[3501:4500] <- colSums(mammals_trop_am[3501:4500], na.rm=T)
mammals_trop_am_Total[4501:5500] <- colSums(mammals_trop_am[4501:5500], na.rm=T)
mammals_trop_am_Total[5501:5741] <- colSums(mammals_trop_am[5501:5741], na.rm=T)
mammals_trop_am_Total[mammals_trop_am_Total>0] <- 1
MM_rich[1,2] <- sum(mammals_trop_am_Total)

mammals_trop_af <- mammals_trop[mammals_trop$Continent==2,]
mammals_trop_af_Total <- rep(0, 5741)
mammals_trop_af_Total[3:500] <- colSums(mammals_trop_af[3:500], na.rm=T)
mammals_trop_af_Total[501:1500] <- colSums(mammals_trop_af[501:1500], na.rm=T)
mammals_trop_af_Total[1501:2500] <- colSums(mammals_trop_af[1501:2500], na.rm=T)
mammals_trop_af_Total[2501:3500] <- colSums(mammals_trop_af[2501:3500], na.rm=T)
mammals_trop_af_Total[3501:4500] <- colSums(mammals_trop_af[3501:4500], na.rm=T)
mammals_trop_af_Total[4501:5500] <- colSums(mammals_trop_af[4501:5500], na.rm=T)
mammals_trop_af_Total[5501:5741] <- colSums(mammals_trop_af[5501:5741], na.rm=T)
mammals_trop_af_Total[mammals_trop_af_Total>0] <- 1
MM_rich[1,3] <- sum(mammals_trop_af_Total)

mammals_trop_eu <- mammals_trop[mammals_trop$Continent==3,]
mammals_trop_eu_Total <- rep(0, 5741)
mammals_trop_eu_Total[3:500] <- colSums(mammals_trop_eu[3:500], na.rm=T)
mammals_trop_eu_Total[501:1500] <- colSums(mammals_trop_eu[501:1500], na.rm=T)
mammals_trop_eu_Total[1501:2500] <- colSums(mammals_trop_eu[1501:2500], na.rm=T)
mammals_trop_eu_Total[2501:3500] <- colSums(mammals_trop_eu[2501:3500], na.rm=T)
mammals_trop_eu_Total[3501:4500] <- colSums(mammals_trop_eu[3501:4500], na.rm=T)
mammals_trop_eu_Total[4501:5500] <- colSums(mammals_trop_eu[4501:5500], na.rm=T)
mammals_trop_eu_Total[5501:5741] <- colSums(mammals_trop_eu[5501:5741], na.rm=T)
mammals_trop_eu_Total[mammals_trop_eu_Total>0] <- 1
MM_rich[1,4] <- sum(mammals_trop_eu_Total)

mammals_arid <- mammals_db[mammals_db$Clim_zone==2,]
mammals_arid_Total <- rep(0, 5741)
mammals_arid_Total[3:500] <- colSums(mammals_arid[3:500], na.rm=T)
mammals_arid_Total[501:1500] <- colSums(mammals_arid[501:1500], na.rm=T)
mammals_arid_Total[1501:2500] <- colSums(mammals_arid[1501:2500], na.rm=T)
mammals_arid_Total[2501:3500] <- colSums(mammals_arid[2501:3500], na.rm=T)
mammals_arid_Total[3501:4500] <- colSums(mammals_arid[3501:4500], na.rm=T)
mammals_arid_Total[4501:5500] <- colSums(mammals_arid[4501:5500], na.rm=T)
mammals_arid_Total[5501:5741] <- colSums(mammals_arid[5501:5741], na.rm=T)
mammals_arid_Total[mammals_arid_Total>0] <- 1
MM_rich[2,1] <- sum(mammals_arid_Total)

mammals_arid_am <- mammals_arid[mammals_arid$Continent==1,]
mammals_arid_am_Total <- rep(0, 5741)
mammals_arid_am_Total[3:500] <- colSums(mammals_arid_am[3:500], na.rm=T)
mammals_arid_am_Total[501:1500] <- colSums(mammals_arid_am[501:1500], na.rm=T)
mammals_arid_am_Total[1501:2500] <- colSums(mammals_arid_am[1501:2500], na.rm=T)
mammals_arid_am_Total[2501:3500] <- colSums(mammals_arid_am[2501:3500], na.rm=T)
mammals_arid_am_Total[3501:4500] <- colSums(mammals_arid_am[3501:4500], na.rm=T)
mammals_arid_am_Total[4501:5500] <- colSums(mammals_arid_am[4501:5500], na.rm=T)
mammals_arid_am_Total[5501:5741] <- colSums(mammals_arid_am[5501:5741], na.rm=T)
mammals_arid_am_Total[mammals_arid_am_Total>0] <- 1
MM_rich[2,2] <- sum(mammals_arid_am_Total)

mammals_arid_af <- mammals_arid[mammals_arid$Continent==2,]
mammals_arid_af_Total <- rep(0, 5741)
mammals_arid_af_Total[3:500] <- colSums(mammals_arid_af[3:500], na.rm=T)
mammals_arid_af_Total[501:1500] <- colSums(mammals_arid_af[501:1500], na.rm=T)
mammals_arid_af_Total[1501:2500] <- colSums(mammals_arid_af[1501:2500], na.rm=T)
mammals_arid_af_Total[2501:3500] <- colSums(mammals_arid_af[2501:3500], na.rm=T)
mammals_arid_af_Total[3501:4500] <- colSums(mammals_arid_af[3501:4500], na.rm=T)
mammals_arid_af_Total[4501:5500] <- colSums(mammals_arid_af[4501:5500], na.rm=T)
mammals_arid_af_Total[5501:5741] <- colSums(mammals_arid_af[5501:5741], na.rm=T)
mammals_arid_af_Total[mammals_arid_af_Total>0] <- 1
MM_rich[2,3] <- sum(mammals_arid_af_Total)

mammals_arid_eu <- mammals_arid[mammals_arid$Continent==3,]
mammals_arid_eu_Total <- rep(0, 5741)
mammals_arid_eu_Total[3:500] <- colSums(mammals_arid_eu[3:500], na.rm=T)
mammals_arid_eu_Total[501:1500] <- colSums(mammals_arid_eu[501:1500], na.rm=T)
mammals_arid_eu_Total[1501:2500] <- colSums(mammals_arid_eu[1501:2500], na.rm=T)
mammals_arid_eu_Total[2501:3500] <- colSums(mammals_arid_eu[2501:3500], na.rm=T)
mammals_arid_eu_Total[3501:4500] <- colSums(mammals_arid_eu[3501:4500], na.rm=T)
mammals_arid_eu_Total[4501:5500] <- colSums(mammals_arid_eu[4501:5500], na.rm=T)
mammals_arid_eu_Total[5501:5741] <- colSums(mammals_arid_eu[5501:5741], na.rm=T)
mammals_arid_eu_Total[mammals_arid_eu_Total>0] <- 1
MM_rich[2,4] <- sum(mammals_arid_eu_Total)

mammals_temp <- mammals_db[mammals_db$Clim_zone==3,]
mammals_temp_Total <- rep(0, 5741)
mammals_temp_Total[3:500] <- colSums(mammals_temp[3:500], na.rm=T)
mammals_temp_Total[501:1500] <- colSums(mammals_temp[501:1500], na.rm=T)
mammals_temp_Total[1501:2500] <- colSums(mammals_temp[1501:2500], na.rm=T)
mammals_temp_Total[2501:3500] <- colSums(mammals_temp[2501:3500], na.rm=T)
mammals_temp_Total[3501:4500] <- colSums(mammals_temp[3501:4500], na.rm=T)
mammals_temp_Total[4501:5500] <- colSums(mammals_temp[4501:5500], na.rm=T)
mammals_temp_Total[5501:5741] <- colSums(mammals_temp[5501:5741], na.rm=T)
mammals_temp_Total[mammals_temp_Total>0] <- 1
MM_rich[3,1] <- sum(mammals_temp_Total)

mammals_temp_am <- mammals_temp[mammals_temp$Continent==1,]
mammals_temp_am_Total <- rep(0, 5741)
mammals_temp_am_Total[3:500] <- colSums(mammals_temp_am[3:500], na.rm=T)
mammals_temp_am_Total[501:1500] <- colSums(mammals_temp_am[501:1500], na.rm=T)
mammals_temp_am_Total[1501:2500] <- colSums(mammals_temp_am[1501:2500], na.rm=T)
mammals_temp_am_Total[2501:3500] <- colSums(mammals_temp_am[2501:3500], na.rm=T)
mammals_temp_am_Total[3501:4500] <- colSums(mammals_temp_am[3501:4500], na.rm=T)
mammals_temp_am_Total[4501:5500] <- colSums(mammals_temp_am[4501:5500], na.rm=T)
mammals_temp_am_Total[5501:5741] <- colSums(mammals_temp_am[5501:5741], na.rm=T)
mammals_temp_am_Total[mammals_temp_am_Total>0] <- 1
MM_rich[3,2] <- sum(mammals_temp_am_Total)

mammals_temp_af <- mammals_temp[mammals_temp$Continent==2,]
mammals_temp_af_Total <- rep(0, 5741)
mammals_temp_af_Total[3:500] <- colSums(mammals_temp_af[3:500], na.rm=T)
mammals_temp_af_Total[501:1500] <- colSums(mammals_temp_af[501:1500], na.rm=T)
mammals_temp_af_Total[1501:2500] <- colSums(mammals_temp_af[1501:2500], na.rm=T)
mammals_temp_af_Total[2501:3500] <- colSums(mammals_temp_af[2501:3500], na.rm=T)
mammals_temp_af_Total[3501:4500] <- colSums(mammals_temp_af[3501:4500], na.rm=T)
mammals_temp_af_Total[4501:5500] <- colSums(mammals_temp_af[4501:5500], na.rm=T)
mammals_temp_af_Total[5501:5741] <- colSums(mammals_temp_af[5501:5741], na.rm=T)
mammals_temp_af_Total[mammals_temp_af_Total>0] <- 1
MM_rich[3,3] <- sum(mammals_temp_af_Total)

mammals_temp_eu <- mammals_temp[mammals_temp$Continent==3,]
mammals_temp_eu_Total <- rep(0, 5741)
mammals_temp_eu_Total[3:500] <- colSums(mammals_temp_eu[3:500], na.rm=T)
mammals_temp_eu_Total[501:1500] <- colSums(mammals_temp_eu[501:1500], na.rm=T)
mammals_temp_eu_Total[1501:2500] <- colSums(mammals_temp_eu[1501:2500], na.rm=T)
mammals_temp_eu_Total[2501:3500] <- colSums(mammals_temp_eu[2501:3500], na.rm=T)
mammals_temp_eu_Total[3501:4500] <- colSums(mammals_temp_eu[3501:4500], na.rm=T)
mammals_temp_eu_Total[4501:5500] <- colSums(mammals_temp_eu[4501:5500], na.rm=T)
mammals_temp_eu_Total[5501:5741] <- colSums(mammals_temp_eu[5501:5741], na.rm=T)
mammals_temp_eu_Total[mammals_temp_eu_Total>0] <- 1
MM_rich[3,4] <- sum(mammals_temp_eu_Total)

mammals_cold <- mammals_db[mammals_db$Clim_zone==4,]
mammals_cold_Total <- rep(0, 5741)
mammals_cold_Total[3:500] <- colSums(mammals_cold[3:500], na.rm=T)
mammals_cold_Total[501:1500] <- colSums(mammals_cold[501:1500], na.rm=T)
mammals_cold_Total[1501:2500] <- colSums(mammals_cold[1501:2500], na.rm=T)
mammals_cold_Total[2501:3500] <- colSums(mammals_cold[2501:3500], na.rm=T)
mammals_cold_Total[3501:4500] <- colSums(mammals_cold[3501:4500], na.rm=T)
mammals_cold_Total[4501:5500] <- colSums(mammals_cold[4501:5500], na.rm=T)
mammals_cold_Total[5501:5741] <- colSums(mammals_cold[5501:5741], na.rm=T)
mammals_cold_Total[mammals_cold_Total>0] <- 1
MM_rich[4,1] <- sum(mammals_cold_Total)

mammals_cold_am <- mammals_cold[mammals_cold$Continent==1,]
mammals_cold_am_Total <- rep(0, 5741)
mammals_cold_am_Total[3:500] <- colSums(mammals_cold_am[3:500], na.rm=T)
mammals_cold_am_Total[501:1500] <- colSums(mammals_cold_am[501:1500], na.rm=T)
mammals_cold_am_Total[1501:2500] <- colSums(mammals_cold_am[1501:2500], na.rm=T)
mammals_cold_am_Total[2501:3500] <- colSums(mammals_cold_am[2501:3500], na.rm=T)
mammals_cold_am_Total[3501:4500] <- colSums(mammals_cold_am[3501:4500], na.rm=T)
mammals_cold_am_Total[4501:5500] <- colSums(mammals_cold_am[4501:5500], na.rm=T)
mammals_cold_am_Total[5501:5741] <- colSums(mammals_cold_am[5501:5741], na.rm=T)
mammals_cold_am_Total[mammals_cold_am_Total>0] <- 1
MM_rich[4,2] <- sum(mammals_cold_am_Total)

mammals_cold_af <- mammals_cold[mammals_cold$Continent==2,]
mammals_cold_af_Total <- rep(0, 5741)
mammals_cold_af_Total[3:500] <- colSums(mammals_cold_af[3:500], na.rm=T)
mammals_cold_af_Total[501:1500] <- colSums(mammals_cold_af[501:1500], na.rm=T)
mammals_cold_af_Total[1501:2500] <- colSums(mammals_cold_af[1501:2500], na.rm=T)
mammals_cold_af_Total[2501:3500] <- colSums(mammals_cold_af[2501:3500], na.rm=T)
mammals_cold_af_Total[3501:4500] <- colSums(mammals_cold_af[3501:4500], na.rm=T)
mammals_cold_af_Total[4501:5500] <- colSums(mammals_cold_af[4501:5500], na.rm=T)
mammals_cold_af_Total[5501:5741] <- colSums(mammals_cold_af[5501:5741], na.rm=T)
mammals_cold_af_Total[mammals_cold_af_Total>0] <- 1
MM_rich[4,3] <- sum(mammals_cold_af_Total)

mammals_cold_eu <- mammals_cold[mammals_cold$Continent==3,]
mammals_cold_eu_Total <- rep(0, 5741)
mammals_cold_eu_Total[3:500] <- colSums(mammals_cold_eu[3:500], na.rm=T)
mammals_cold_eu_Total[501:1500] <- colSums(mammals_cold_eu[501:1500], na.rm=T)
mammals_cold_eu_Total[1501:2500] <- colSums(mammals_cold_eu[1501:2500], na.rm=T)
mammals_cold_eu_Total[2501:3500] <- colSums(mammals_cold_eu[2501:3500], na.rm=T)
mammals_cold_eu_Total[3501:4500] <- colSums(mammals_cold_eu[3501:4500], na.rm=T)
mammals_cold_eu_Total[4501:5500] <- colSums(mammals_cold_eu[4501:5500], na.rm=T)
mammals_cold_eu_Total[5501:5741] <- colSums(mammals_cold_eu[5501:5741], na.rm=T)
mammals_cold_eu_Total[mammals_cold_eu_Total>0] <- 1
MM_rich[4,4] <- sum(mammals_cold_eu_Total)

mammals_pola <- mammals_db[mammals_db$Clim_zone==5,]
mammals_pola_Total <- rep(0, 5741)
mammals_pola_Total[3:500] <- colSums(mammals_pola[3:500], na.rm=T)
mammals_pola_Total[501:1500] <- colSums(mammals_pola[501:1500], na.rm=T)
mammals_pola_Total[1501:2500] <- colSums(mammals_pola[1501:2500], na.rm=T)
mammals_pola_Total[2501:3500] <- colSums(mammals_pola[2501:3500], na.rm=T)
mammals_pola_Total[3501:4500] <- colSums(mammals_pola[3501:4500], na.rm=T)
mammals_pola_Total[4501:5500] <- colSums(mammals_pola[4501:5500], na.rm=T)
mammals_pola_Total[5501:5741] <- colSums(mammals_pola[5501:5741], na.rm=T)
mammals_pola_Total[mammals_pola_Total>0] <- 1
MM_rich[5,1] <- sum(mammals_pola_Total)

mammals_pola_am <- mammals_pola[mammals_pola$Continent==1,]
mammals_pola_am_Total <- rep(0, 5741)
mammals_pola_am_Total[3:500] <- colSums(mammals_pola_am[3:500], na.rm=T)
mammals_pola_am_Total[501:1500] <- colSums(mammals_pola_am[501:1500], na.rm=T)
mammals_pola_am_Total[1501:2500] <- colSums(mammals_pola_am[1501:2500], na.rm=T)
mammals_pola_am_Total[2501:3500] <- colSums(mammals_pola_am[2501:3500], na.rm=T)
mammals_pola_am_Total[3501:4500] <- colSums(mammals_pola_am[3501:4500], na.rm=T)
mammals_pola_am_Total[4501:5500] <- colSums(mammals_pola_am[4501:5500], na.rm=T)
mammals_pola_am_Total[5501:5741] <- colSums(mammals_pola_am[5501:5741], na.rm=T)
mammals_pola_am_Total[mammals_pola_am_Total>0] <- 1
MM_rich[5,2] <- sum(mammals_pola_am_Total)

mammals_pola_af <- mammals_pola[mammals_pola$Continent==2,]
mammals_pola_af_Total <- rep(0, 5741)
mammals_pola_af_Total[3:500] <- colSums(mammals_pola_af[3:500], na.rm=T)
mammals_pola_af_Total[501:1500] <- colSums(mammals_pola_af[501:1500], na.rm=T)
mammals_pola_af_Total[1501:2500] <- colSums(mammals_pola_af[1501:2500], na.rm=T)
mammals_pola_af_Total[2501:3500] <- colSums(mammals_pola_af[2501:3500], na.rm=T)
mammals_pola_af_Total[3501:4500] <- colSums(mammals_pola_af[3501:4500], na.rm=T)
mammals_pola_af_Total[4501:5500] <- colSums(mammals_pola_af[4501:5500], na.rm=T)
mammals_pola_af_Total[5501:5741] <- colSums(mammals_pola_af[5501:5741], na.rm=T)
mammals_pola_af_Total[mammals_pola_af_Total>0] <- 1
MM_rich[5,3] <- sum(mammals_pola_af_Total)

mammals_pola_eu <- mammals_pola[mammals_pola$Continent==3,]
mammals_pola_eu_Total <- rep(0, 5741)
mammals_pola_eu_Total[3:500] <- colSums(mammals_pola_eu[3:500], na.rm=T)
mammals_pola_eu_Total[501:1500] <- colSums(mammals_pola_eu[501:1500], na.rm=T)
mammals_pola_eu_Total[1501:2500] <- colSums(mammals_pola_eu[1501:2500], na.rm=T)
mammals_pola_eu_Total[2501:3500] <- colSums(mammals_pola_eu[2501:3500], na.rm=T)
mammals_pola_eu_Total[3501:4500] <- colSums(mammals_pola_eu[3501:4500], na.rm=T)
mammals_pola_eu_Total[4501:5500] <- colSums(mammals_pola_eu[4501:5500], na.rm=T)
mammals_pola_eu_Total[5501:5741] <- colSums(mammals_pola_eu[5501:5741], na.rm=T)
mammals_pola_eu_Total[mammals_pola_eu_Total>0] <- 1
MM_rich[5,4] <- sum(mammals_pola_eu_Total)

save(MM_rich, file = "MM_rich.RData")

# unique(mammals_pola$Loxodonta_africana)
# mammals_db$Continent[c(38199, 38200, 38201, 38202, 38203, 38204, 38205, 38206, 38207)]
# mammals_db$Clim_zone[c(38199, 38200, 38201, 38202, 38203, 38204, 38205, 38206, 38207)]


quartz()

spe_mtx_am <- matrix(data=c(309, 42, 95, 3, 6), ncol = 1, nrow = 5) 
rownames(spe_mtx_am) <- c("Tropical", "Arid", "Temperate", "Cold", "Polar")
colnames(spe_mtx_am) <- "N"

spe_mtx_af <- matrix(data=c(276, 106, 41, 0, 0), ncol = 1, nrow = 5) 
rownames(spe_mtx_af) <- c("Tropical", "Arid", "Temperate", "Cold", "Polar")
colnames(spe_mtx_af) <- "N"

spe_mtx_eu <- matrix(data=c(629, 94,134, 30, 5), ncol = 1, nrow = 5) 
rownames(spe_mtx_eu) <- c("Tropical", "Arid", "Temperate", "Cold", "Polar")
colnames(spe_mtx_eu) <- "N"


tot_mtx_am <- matrix(data=c(309*100/(1117+309), 42*100/(42+1411), 95*100/(95+1380), 3*100/(3+745), 6*100/(6+421)), ncol = 1, nrow = 5) 
rownames(tot_mtx_am) <- c("Tropical", "Arid", "Temperate", "Cold", "Polar")
colnames(tot_mtx_am) <- "%"

tot_mtx_af <- matrix(data=c(276*100/(276+907), 106*100/(106+1043), 41*100/(41+858), 0, 0), ncol = 1, nrow = 5) 
rownames(tot_mtx_af) <- c("Tropical", "Arid", "Temperate", "Cold", "Polar")
colnames(tot_mtx_af) <- "%"

tot_mtx_eu <- matrix(data=c(629*100/(629+882), 94*100/(1346+94),134*100/(134+1256), 30*100/(30+759), 5*100/(5+522)), ncol = 1, nrow = 5) 
rownames(tot_mtx_eu) <- c("Tropical", "Arid", "Temperate", "Cold", "Polar")
colnames(tot_mtx_eu) <- "%"


barplot(height = spe_mtx_af, beside = TRUE, col=c("#f9d14a", "#ab3329", "#ed968c", "#7c4b73", "#88a0dc"),
        border="white", legend.text=NULL,
        ylab= "N specialist species")
quartz()

library(waffle)


waffle(spe_mtx_am, rows = 21,
       colors=c("#f9d14a", "#ab3329", "#ed968c", "#7c4b73", "#88a0dc"))

aaa <- data.frame(spe_mtx_eu)
aaa$group <- as.factor(rownames(aaa))
aaa$fraction = aaa$N / sum(aaa$N)
aaa$ymax = cumsum(aaa$fraction)
aaa$ymin = c(0, head(aaa$ymax, n=-1))
aaa$labelPosition <- (aaa$ymax + aaa$ymin) / 2
aaa$labelPosition[4] <- 0.96
aaa$labelPosition[5] <- 0.02
aaa$label <- aaa$N

ggplot(aaa,aes(values=N, fill=group)) +
  geom_waffle(n_rows = 21, size = 0.33, colour = "white") +
  scale_fill_manual(name = NULL,
                    values = c("#ab3329", "#7c4b73", "#88a0dc","#ed968c","#f9d14a")) +
  coord_equal() +
  theme_void()

quartz()

ggplot(aaa, aes(ymax=ymax, ymin=ymin, xmax=4, xmin=3, fill=group)) +
  geom_rect() +
  geom_label( x=3.5, aes(y=labelPosition, label=label), size=6) +
  coord_polar(theta="y") + # Try to remove that to understand how the chart is built initially
  xlim(c(1, 4)) + # Try to remove that to see how to make a pie chart
  scale_fill_manual(name = NULL,
                    values = c("#ab3329", "#7c4b73", "#88a0dc","#ed968c","#f9d14a")) +
  theme_void() +
  theme(legend.position = "none")

################################
#### Correlation analyses ###
################################
library("Hmisc")
setwd("/Users/saragamboa/Dropbox/Mac/Documents/Sara Gamboa/Papers/Fragmentation")
corr_db <- read.table("rich_con2.csv",sep=";",header=T)
rownames(corr_db) <- corr_db[,1]
corr_db <- corr_db[-1]

ggplot(corr_db)+
  geom_point(aes(Spp,NT_S))+geom_smooth(aes(Spp,NT_S),method="lm",se=FALSE) +
  geom_point(aes(Spp,NT_M))+geom_smooth(aes(Spp,NT_M),method="lm",se=FALSE)+
  geom_point(aes(Spp,NT_L))+geom_smooth(aes(Spp,NT_L),method="lm",se=FALSE) +
  geom_point(aes(Spp,NT_XL))+geom_smooth(aes(Spp,NT_XL),method="lm",se=FALSE)


require(ggplot2)
require(reshape2)
corr_db2 = melt(corr_db, id.vars='Spp')

kkk <- corr_db2[1:52,]
kkk <- corr_db2[53:104,]
kkk <- corr_db2[105:156,]
kkk <- corr_db2[157:208,]
ttt <- corr_db2[209:234,]

ggplot(kkk) +
  geom_jitter(aes(value,Spp, colour=variable),) +
  geom_smooth(aes(value,Spp, colour=variable), method=lm, se=FALSE) +
  facet_wrap(~variable, scales="free_x") +
  labs(x = "Fragmentation", y = "Number of species (N)")


ggplot(as.data.frame(residuals(modelo)),aes(sample = residuals(modelo))) +  geom_qq() 

#####WLS#####
fit=lm(Spp~P_Total,data=corr_db)
plot(fitted(fit), resid(fit), xlab='Fitted Values', ylab='Residuals')
library(lmtest)
bptest(fit)
wt <- 1 / lm(abs(fit$residuals) ~ fit$fitted.values)$fitted.values^2
wls_model <- lm(Spp ~ P_Total, data = corr_db, weights=wt)
summary(wls_model)

ggplot(kkk,aes(y=Spp,x=value,color=variable))+geom_point()+stat_smooth(method="lm",se=FALSE)

fit=glm(Spp~N_eventos_S+N_eventos_M+N_eventos_L+N_eventos_XL+NT_XL,data=corr_db)
fit=glm(Spp~N_eventos_S+N_eventos_M+N_eventos_L+N_eventos_XL+NT_XL+NT_L+NT_M+NT_S+Strengh_XL+Strengh_L+Strengh_M+Strengh_S+Max_S+Max_M+Max_L+Max_XL,data=corr_db)

fit=glm(Spp~N_eventos_S+N_eventos_M+N_eventos_L+N_eventos_XL+NT_XL+NT_L+NT_M+NT_S+Strengh_M+Strengh_S+Max_S+Max_M,data=corr_db)
fit=glm(Spp~N_eventos_S+N_eventos_M+N_eventos_L+N_eventos_XL+NT_XL+NT_L+NT_M+NT_S+Strengh_L+Strengh_XL+Max_L+Max_XL,data=corr_db)


######Let's include climate
setwd("~/Documents/Sara Gamboa/Papers/Fragmentation/clima")
Temp <- Time_0_temp_array
Prec <- Time_0_prec_array

dim(Prec)

library(abind)

Temp2<- apply(Temp, 1:2, mean)
Prec2 <- apply(Prec, 1:2, sum)


toMol <- function(matrix){
  r_raster <- raster(matrix)
  crs(r_raster) <- CRS("+init=epsg:4326")
  extent(r_raster) <- extent(-180,180,-90,90)
  crs.mol <- CRS("+proj=moll +lon_0=0 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs")
  raster_mol <- projectRaster(r_raster, crs=crs.mol)
  return(raster_mol)
}

Temp_mol <- toMol(Temp2)
Prec_mol <- toMol(Prec2)
plot(Temp_mol)
plot(mapa_Af)

CR_data <- function(Clim,Area){
salida <- vector(mode="double", length = 5)
for (i in 1:5) {
  kkk <- Clim[Area==i]
  salida[i] <- mean(kkk)
  print(i)
}
return(salida)
}

Africa_T <- CR_data(Temp_mol,mapa_Af)
Africa_P <- CR_data(Prec_mol,mapa_Af)
America_T <- CR_data(Temp_mol,mapa_Am)
America_P <- CR_data(Prec_mol,mapa_Am)
Euroc_T <- CR_data(Temp_mol,mapa_E)
Euroc_P <- CR_data(Prec_mol,mapa_E)

setwd("/Users/saragamboa/Dropbox/Mac/Documents/Sara Gamboa/Papers/Fragmentation")
corr_db <- read.table("rich_con2.csv",sep=";",header=T)
fit=glm(Spp ~ T_Mean + P_Total ,data=corr_db) ##72%
fit=glm(Spp ~ T_Mean ,data=corr_db) ##72%
fit=glm(Spp ~ P_Total ,data=corr_db) ##23%
fit=glm(Spp ~ N_eventos_S + N_eventos_M + N_eventos_L + N_eventos_XL + NT_XL, data = corr_db)##85%
fit=glm(Spp ~ N_eventos_S + N_eventos_M + N_eventos_L + N_eventos_XL + NT_XL + T_Mean + P_Total, data = corr_db) ##92%

####Especialistas
setwd("/Users/saragamboa/Dropbox/Mac/Documents/Sara Gamboa/Papers/Fragmentation")
corr_db <- read.table("espe_con.csv",sep=";",header=T)
rownames(corr_db) <- corr_db[,1]
corr_db <- corr_db[-1]


fit=glm(Spp ~ T_Mean + P_Total ,data=corr_db) ##72%
fit=glm(Spp ~ T_Mean ,data=corr_db) ##72%
fit=glm(Spp ~ P_Total ,data=corr_db) ##23%
fit=glm(Spp ~ N_eventos_S + N_eventos_M + N_eventos_L + N_eventos_XL + NT_XL, data = corr_db)##85%
fit=glm(Spp ~ N_eventos_S + N_eventos_M + N_eventos_L + N_eventos_XL + NT_XL + T_Mean + P_Total, data = corr_db) #


fit=lm(Spp~NT_S,data=corr_db)
wt <- 1 / lm(abs(fit$residuals) ~ fit$fitted.values)$fitted.values^2
wls_model <- lm(Spp ~ NT_S, data = corr_db, weights=wt)
summary(wls_model)

fit <- glm(Spp ~ T_Mean + P_Total ,data=corr_db) ##76%
fit <- glm(Spp ~ T_Mean ,data=corr_db) ##40%
fit <- glm(Spp ~ P_Total ,data=corr_db) ##72%
fit <- glm(Spp ~ N_eventos_S + N_eventos_M + NT_S + NT_M, data = corr_db)##92%
fit <- glm(Spp ~ N_eventos_S + N_eventos_M + NT_S + NT_M + Strengh_S + Strengh_M + Max_S + Max_M, data = corr_db)##92%

modelo <- step(glm(Spp ~ N_eventos_S + N_eventos_M + NT_S + NT_M + Strengh_S + Strengh_M + Max_S + Max_M, data = corr_db), direction="both")

fit <- glm(Spp ~ N_eventos_S +  NT_M, data = corr_db)
fit <- glm(Spp ~ N_eventos_S +  NT_M +T_Mean, data = corr_db)
summary(fit)

row.names(corr_db)

plot(corr_db$NT_M, corr_db$Spp, col=c(1,2,3,4,5,1,2,3,1,2,3,4,5), pch=16)

quartz()

res <- cor(corr_db)
res_k <- cor(corr_db_k)
cor(corr_db, use = "complete.obs")
res <- round(res, 2)
res_k <- round(res_k, 2)

dput(corr_db)
names(corr_db_k) <- c("Spp","Patches_S","Patches_M","Patches_L","Patches_XL","Events_S","Events_M","Events_L","Events_XL","Strengh_S","Strength_M","Strength_L","Strength_XL")

corr_db_k <- corr_db[,1:13]
corr_m <- as.matrix(corr_db)
corr_m_k <- as.matrix(corr_db_k)

res_p <- rcorr(corr_m, type = c("pearson","spearman"))
res_p_k <- rcorr(corr_m_k, type = c("pearson","spearman"))
save(res_p,file = "correlation_p_matrix.Rdata")

flattenCorrMatrix <- function(cormat, pmat) {
  ut <- upper.tri(cormat)
  data.frame(
    row = rownames(cormat)[row(cormat)[ut]],
    column = rownames(cormat)[col(cormat)[ut]],
    cor  =(cormat)[ut],
    p = pmat[ut]
  )
}
FCM <- flattenCorrMatrix(res_p_k$r, res_p_k$P)
Spp_FCM <- FCM[FCM$row=="Spp",]

symnum(res, abbr.colnames = FALSE)

library(corrplot)
res_p$P[is.na(res_p$P)==T] <- 0
res_p_k$P[is.na(res_p_k$P)==T] <- 0
corrplot(res_p$r, type="upper", order="hclust", 
         p.mat = res_p$P, sig.level = 0.05, insig = "blank")
corrplot(res_p$r, type="upper", order="hclust", 
         p.mat = res_p$P, sig.level = 0.05, insig = "blank", method="pie")
corrplot(res_p$r, type="upper", order="hclust", 
         p.mat = res_p$P, sig.level = 0.05, insig = "blank", method="color")
corrplot(res_p_k$r, type="upper", order="hclust", 
         p.mat = res_p_k$P, sig.level = 0.05, insig = "blank", method="color")

coloursv <- Up_Down$q.val

Spp_FCM
library(corrr)
corr_db_k %>% correlate()
corr_db_k %>% correlate() %>% focus(Spp) %>%
  ggplot(aes(x = term, y = Spp)) +
  geom_bar(stat = "identity") +
  ylab("Correlation with number of species") +
  xlab("Variable")

quartz()
  ggplot(Spp_FCM, aes(x = column, y = cor, fill=p)) +
  geom_bar(stat = "identity") +
  ylab("Correlation with number of species") +
  xlab("Variable")
  citation("landscapemetrics")
  

