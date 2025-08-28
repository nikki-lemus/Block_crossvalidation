remotes::install_github("marlonecobos/evniche")

library(evniche)
library(geodata)

# downloading environmental variables
wrclim <- worldclim_global(var = "bio", res = 5, path = "Data")

# keeping only temperature and precipitation and masking to SA (BIO 1 and 12)
## SA extent
SAext <- ext(c(-85, -30, -60, 15))

## cropping and renaming
bio_SA <- crop(wrclim[[c(1, 12)]], SAext)
names(bio_SA) <- c("Temperature", "Precipitation")

plot(bio_SA)

# preparing data for analyses (from raster to matrix)
data_T_P <- as.data.frame(bio_SA, xy = TRUE)

# plot to explore convinient ranges
evniche:::plot_2D(background = data_T_P[, 3:4], 
                  lp_list = list(lty = 1, lwd = 2, col_l = "red",
                                 alpha_l = 1, pch = 1, cex_p = 1,
                                 col_p = "white", alpha_p = 0.5))


# niche range (limits)
host_niche_range <- cbind(Temperature = c(12, 26), Precipitation = c(700, 2800))
host_niche_range

# host virtual niche
## variances
vars <- evniche:::var_from_range(range = host_niche_range)

## covariance limits
cov_lim <- evniche:::covariance_limits(range = host_niche_range)

## variance-covariance matrix
cov <- cov_lim$max_covariance * 0.1 # covariance selected
varcov <- evniche:::var_cov_matrix(variances = vars, covariances = cov) 

## centroid
cent <- evniche:::centroid(range = host_niche_range)

## ellipsoid characteristics (virtual niche)
host_niche <- ell_features(centroid = cent, covariance_matrix = varcov,
                           level = 0.99)

## niche plot in E space
evniche:::plot_2D(features = host_niche, background = data_T_P[, 3:4], 
                  lp_list = list(lty = 1, lwd = 2, col_l = "red",
                                 alpha_l = 1, pch = 1, cex_p = 1,
                                 col_p = "white", alpha_p = 0.5))

# generate host data
## predict suitability on bioclimatic data based on ellipsoids
pred_geo <- evniche:::ell_predict(data = bio_SA, features = host_niche)

plot(pred_geo$suitability_trunc)

## write prediction (I have a folder named Results in my directory)
writeRaster(pred_geo$suitability_trunc, filename = "Results/prediction.tif")

## perform explorations to define ecorregions for accessible area in a GIS
## layer with selected ecoregions (I exported selected ecoregions in QGIS)
ecor <- vect("Data/vector/selected_eco_mec.gpkg")

## 50 km buffer 
ecor50 <- buffer(ecor, width = 50000, quadsegs = 100)
ecor50 <- aggregate(ecor50)

writeVector(ecor50, filename = "Data/vector/accessible_area.gpkg")

## check with a plot
plot(pred_geo$suitability_trunc)
lines(ecor50, col = "red", lwd = 2)

## mask variables to buffer
bio_acc <- crop(bio_SA, ecor50, mask = TRUE)

## bios to data frame
data_acc <- as.data.frame(bio_acc, xy = TRUE)

## predict again on data
pred_acc <- evniche:::ell_predict(data = data_acc, longitude = "x", 
                                  latitude = "y", features = host_niche)

## generate new data
## based on the ellipsoid and available conditions
set.seed(1)
vd_pre_host <- virtual_data(features = host_niche, from = "prediction",
                            data = data_acc, prediction = pred_acc, n = 200)

points(vd_pre_host[, 1:2], col = "red")

## save data
write.csv(vd_pre_host, file = "Results/data_virtual_species.csv", 
          row.names = FALSE)


