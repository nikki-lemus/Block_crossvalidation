################################################################################
# Project: Data partitioning methods and model transferability
# Script: Initial exploration in geographic regions
# Author(s): Marlon E. Cobos 
# Date: 10/07/2025  # dd/mm/yyyy
################################################################################

# Note: If you are working on this project, I recommend you make a copy of this 
# script to start exploring with your examples.


# Setting R up -----------------------------------------------------------------
# install packages (if needed)
install.packages("geodata")
install.packages("mop")
install.packages("scales")

install.packages("remotes")
remotes::install_github("marlonecobos/kuenm2")

# load packages needed
library(geodata)
library(mop)
library(kuenm2)
library(scales)
# ------------------------------------------------------------------------------


# Data -------------------------------------------------------------------------
## environmental data
dir.create("Data")

wrclim <- worldclim_global(var = "bio", res = 5, path = "Data")
# ------------------------------------------------------------------------------


# Data preparation and initial analyses ----------------------------------------
# defining clipping area South America
extsa <- ext(c(-82, -34, -56, 13))

sawc <- crop(wrclim, extsa)

plot(sawc$wc2.1_5m_bio_1)

# a smaller area Cerrado-Amazon transition
extca <- ext(c(-57, -43, -20, -5.5))

lines(extca)

cawc <- crop(sawc, extca)

plot(cawc$wc2.1_5m_bio_1)

# partitions
b1 <- ext(c(-57, mean(c(-57, -43)), mean(c(-20, -5.5)), -5.5))
lines(b1)

b2 <- ext(c(mean(c(-57, -43)), -43, mean(c(-20, -5.5)), -5.5))
lines(b2)

b3 <- ext(c(-57, mean(c(-57, -43)), -20, mean(c(-20, -5.5))))
lines(b3)

b4 <- ext(c(mean(c(-57, -43)), -43, -20, mean(c(-20, -5.5))))
lines(b4)          


# PCA
capca <- perform_pca(raster_variables = cawc, center = TRUE, scale = TRUE)
capca2 <- capca$env[[1:2]]
# ------------------------------------------------------------------------------


# Exploratory analysis ---------------------------------------------------------
# MOP from three blocks to all
## prepare area for three blocks
b1_3v <- union(union(vect(b1), vect(b2)), vect(b3))
crs(b1_3v) <- crs(capca2)

plot(capca2$PC1)
plot(b1_3v, col = 2, add = TRUE)

## run MOP
mop1_3_all <- mop(m = crop(capca2, b1_3v, mask = TRUE), g = capca2, 
                  type = "detailed", calculate_distance = TRUE, 
                  where_distance = "all")


# plots to visualize results
x11()
par(mfrow = c(2, 3))

plot(capca2$PC1, main = "PC1")
plot(capca2$PC2, main = "PC2")

plot(mop1_3_all$mop_distances, main = "Distance")
plot(mop1_3_all$mop_basic, main = "Outside ranges")
plot(mop1_3_all$mop_simple, main = "N variables outside")

#plot(mop1_3_all$mop_detailed$towards_low_combined,  # nothing is outside low end
#     main = "Variables outside low end")
plot(mop1_3_all$mop_detailed$towards_high_combined, 
     main = "Variables outside high end")
# ------------------------------------------------------------------------------
