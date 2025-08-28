#remotes::install_github("marlonecobos/evniche")

library(evniche)
library(geodata)

# downloading environmental variables
wrclim <- worldclim_global(var = "bio", res = 5, path = "Data")



# keeping only temperature and precipitation and masking to SA (BIO 1 and 12)
## SA extent
#long,long, lat, lat?
AUext <- ext(c(110, 160, -40, -9.5))

## cropping and renaming
bio_AU <- crop(wrclim[[c(1, 12)]], AUext)
names(bio_AU) <- c("Temperature", "Precipitation")

plot(bio_AU)

# preparing data for analyses (from raster to matrix)
data_T_P <- as.data.frame(bio_AU, xy = TRUE)

plot(data_T_P[, 3:4])


# niche range (limits)
host_niche_range <- cbind(Temperature = c(20, 30), Precipitation = c(500, 2000))
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



# host data
## predict suitability on bioclimatic data based on ellipsoids
pred_geo <- evniche:::ell_predict(data = bio_AU, features = host_niche)

plot(pred_geo$suitability_trunc)

pred_host <- evniche:::ell_predict(data = data_T_P, features = host_niche, 
                                   longitude = "x", latitude = "y")

## generate new data
## based on the ellipsoid and available conditions
set.seed(1)
vd_pre_host <- virtual_data(features = host_niche, from = "prediction",
                            data = data_T_P, prediction = pred_host, n = 200)

points(vd_pre_host[, 1:2], col = "red")

