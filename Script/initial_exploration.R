################################################################################
# Project: Data partitioning methods and model transferability
# Script: Initial exploration in geographic regions
# Author(s): Marlon E. Cobos 
# Modified: 17/07/2025  # dd/mm/yyyy
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

# loading functions
source("Script/functions.R")
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
blist <- block_extents4(ext = extca)

# PCA
capca <- perform_pca(raster_variables = cawc, center = TRUE, scale = TRUE)
capca2 <- capca$env[[1:2]]
# ------------------------------------------------------------------------------


# Exploratory analysis ---------------------------------------------------------
# MOP from three blocks to all
area1 <- mop_comb(block_list = blist, vars = capca2, calculate_distance = TRUE)

# comparing mop distances for every block left out
ref <- area1[[1]]$b3_comb
int <- blist[[4]]

disarea_ref <- crop(area1[[1]]$mop_distances, ref, mask = TRUE)
disarea_int <- crop(area1[[1]]$mop_distances, int)

disval_ref <- as.data.frame(disarea_ref, cell = FALSE)
disval_int <- as.data.frame(disarea_int, cell = FALSE)

bxvals <- rbind(disval_ref, disval_int)
bxvals$Areas <- c(rep("Reference", nrow(disval_ref)), rep("Left_out", nrow(disval_int)))
bxvals$Areas <- as.factor(bxvals$Areas)
bxvals$Distance <- as.numeric(bxvals$Distance)

colnames(bxvals)[1] <- "Distance"

head(bxvals)

boxplot(data.frame(Reference = disval_ref[, 1], Left_out = disval_int[, 1]), 
        xlab = "Areas", ylab = "Distances", las = 1)

colMeans(data.frame(Reference = disval_ref[, 1], Left_out = disval_int[, 1]))

# plot in environmental space blocks 1-3 vs block4
x11()
plot(as.data.frame(crop(capca2, area1$Comb_123$b3_comb, mask = TRUE), cells = FALSE),
     col = "gray2", pch = 16, xlim = c(-7.3, 8.6), ylim = c(-6.1, 7.5))
points(as.data.frame(crop(capca2, blist$block_4), cells = FALSE),
       col = alpha("red", 0.25), pch = 3)


# plots to visualize results
x11()
par(mfrow = c(2, 3))

plot(capca2$PC1, main = "PC1")
plot(capca2$PC2, main = "PC2")

plot(area1$Comb_123$mop_distances, main = "Distance")
plot(area1$Comb_123$mop_basic, main = "Outside ranges")
plot(area1$Comb_123$mop_simple, main = "N variables outside")

plot(area1$Comb_123$mop_detailed$towards_low_combined,  
     main = "Variables outside low end")
plot(area1$Comb_123$mop_detailed$towards_high_combined, 
     main = "Variables outside high end")
# ------------------------------------------------------------------------------
