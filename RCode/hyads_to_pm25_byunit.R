rm( list = ls())

source( '~/Dropbox/Harvard/RFMeval_Local/HyADS_to_pm25/RCode/hyads_to_pm25_functions.R')

#coordinate reference system projection string for spatial data
p4s <- "+proj=lcc +lat_1=33 +lat_2=45 +lat_0=40 +lon_0=-97 +a=6370000 +b=6370000"

#======================================================================#
## Load saved object of cmaq-ddm / hyads models from hyads_to_pm25_month.R
#======================================================================#
load( '~/Dropbox/Harvard/RFMeval_Local/HyADS_to_pm25/RData/hyads_to_cmaq_models.RData')

# load the models 
hyads.ann.model <- preds.ann.hyads06w05$model.gam
idwe.ann.model <- preds.ann.idwe06w05$model.gam

# set up the dataset
dat.coords <- coordinates( dats2006.a)
dats2006raw.dt <- data.table( cbind( dat.coords, values( dats2006.a)))

# do the predictions
hyads2006.pred <- predict( hyads.ann.model, newdata = dats2006raw.dt)
idwe2006.pred <- predict( idwe.ann.model, newdata = dats2006raw.dt)

# rasterize
dats2006.r <- rasterFromXYZ( data.table( dat.coords, hyads2006.pred, idwe2006.pred), crs = p4s)

#======================================================================#
## get total state populations
#======================================================================#
grid_popwgt.xyz <- fread( '~/Dropbox/Harvard/RFMeval_Local/HyADS_to_pm25/HyADS_grid/population/hyads_grid_population.csv',
                          drop = 'V1')

grid_popwgt.r <- rasterFromXYZ( grid_popwgt.xyz, crs = p4s)
grid_popwgt.sf <- st_as_sf( rasterToPolygons( grid_popwgt.r))

# get total state populations
us_states <- st_transform( USAboundaries::us_states(), p4s)
us_states.p <- st_interpolate_aw( grid_popwgt.sf, us_states, extensive = T)

# merge back to state information
us_states.pop <- cbind( us_states.p, us_states[us_states.p$Group.1,])
us_states.pop$geometry.1 <- NULL

#======================================================================#
## area weight over states
#======================================================================#
# combine exposure and population
dats2006_popwgt.r <- project_and_stack( dats2006.r, grid_popwgt.r)

# weight by 2006 population
dats2006_popwgt.names <- names( dats2006_popwgt.r)
dats2006_popwgt.r <- dats2006_popwgt.r * dats2006_popwgt.r$X2006
names( dats2006_popwgt.r) <- dats2006_popwgt.names

# take over states
dats2006_popwgt.sf <- st_as_sf( rasterToPolygons( dats2006_popwgt.r))
plot( dats2006_popwgt.sf[is.na( dats2006_popwgt.sf$hyads2006.pred),]$geometry, add = T)

#NA's coming from 
dats2006_popwgt.sf
dats2006_popwgt.states <- st_interpolate_aw( dats2006_popwgt.sf, us_states, extensive = F)

## need to update masking in all functions - cells with centroid not covered are cropped
##   https://gis.stackexchange.com/questions/255025/r-raster-masking-a-raster-by-polygon-also-remove-cells-partially-covered
## should probably update usa mask too - need to use USAboundaries for consistency
SpP_ras <- rasterize(SpP, r, getCover=TRUE)
SpP_ras[SpP_ras==0] <- NA




