library( censusapi)
library( data.table)
library( USAboundaries)
library( ggplot2)
library( viridis)
library( scales)
library( raster)
library( sf)
library( fasterize)
# https://cran.r-project.org/web/packages/censusapi/vignettes/getting-started.html
# https://www.census.gov/data/developers/data-sets/popest-popproj/popest.2000-2010_Intercensals.html

Sys.setenv( CENSUS_KEY= '6120c82cfe5ed2b6f8dbed17dc946e7b05307a39')
readRenviron("~/.Renviron")

# listCensusMetadata(name = "2000/pep/int_population", type = "variables")
# listCensusMetadata(name = "2000/pep/int_population", type = "geography")

## =========================================================== ##
# get population estimates for counties 
## =========================================================== ##
pop <- data.table( getCensus(name = "2000/pep/int_population", region = "county:*",
                             vars = c("DATE_", "GEONAME", "DATE_DESC", "POP")))
pop[, POP := as.numeric( POP)]
pop.2006 <- pop[grep( '2006 pop', DATE_DESC)][, year := 2006]
pop.2011 <- pop[grep( '2010 pop', DATE_DESC)][, year := 2011]

## =========================================================== ##
# link with spatial data
## =========================================================== ##
us_counties.in <- us_counties()
us_counties.pop06 <- data.table( merge( us_counties.in, pop.2006,
                                        by.x = c( 'statefp', 'countyfp'), by.y = c( 'state', 'county')))
us_counties.pop11 <- data.table( merge( us_counties.in, pop.2011,
                                        by.x = c( 'statefp', 'countyfp'), by.y = c( 'state', 'county')))

## =========================================================== ##
#just lower 48 states
## =========================================================== ##
states.lower <- state.name[ !( state.name %in% c( 'Hawaii', 'Alaska'))]
us_counties.pop06 <- us_counties.pop06[ state_name %in% states.lower]
us_counties.pop11 <- us_counties.pop11[ state_name %in% states.lower]

## =========================================================== ##
# read annual HyADS in 2005
## =========================================================== ##
p4s <- "+proj=lcc +lat_1=33 +lat_2=45 +lat_0=40 +lon_0=-97 +a=6370000 +b=6370000"
load( '~/Dropbox/Harvard/RFMeval_Local/HyADS_to_pm25/RData/hyads_to_cmaq_models.RData')
# 
# # load the models 
# hyads.ann.model <- preds.ann.hyads06w05$model.gam
# idwe.ann.model <- preds.ann.idwe06w05$model.gam
# 
# # set up the dataset
# dat.coords <- coordinates( dats2006.a)
# dats2006raw.dt <- data.table( cbind( dat.coords, values( dats2006.a)))
# 
# # do the predictions
# hyads2005.pred <- predict( hyads.ann.model, newdata = dats2006raw.dt)
# idwe2005.pred <- predict( idwe.ann.model, newdata = dats2006raw.dt)
# 
# # rasterize
# dats2005.r <- rasterFromXYZ( data.table( dat.coords, hyads2005.pred, idwe2005.pred), crs = p4s)

## =========================================================== ##
## # get everything on the same crs
## =========================================================== ##
us_counties06.sf <- st_transform( st_as_sf( us_counties.pop06), crs = p4s)
us_counties11.sf <- st_transform( st_as_sf( us_counties.pop11), crs = p4s)
dats2005.sf <- st_as_sf( rasterToPolygons( dats2006.a))

## =========================================================== ##
## # area weight to hyads grid
## =========================================================== ##
grid_popwgt06 <- st_interpolate_aw( us_counties06.sf['POP'], 
                                          dats2005.sf, extensive = T)
grid_popwgt11 <- st_interpolate_aw( us_counties11.sf['POP'], 
                                          dats2005.sf, extensive = T)

## =========================================================== ##
## # onvert to raster and save
## =========================================================== ##
grid_popwgt06.r <- fasterize( grid_popwgt06, raster( dats2006.a), field = 'POP')
grid_popwgt11.r <- fasterize( grid_popwgt11, raster( dats2006.a), field = 'POP')
grid_popwgt.xyz <- dcast( rbind( data.table( rasterToPoints( grid_popwgt06.r))[, year := 2006],
                         data.table( rasterToPoints( grid_popwgt11.r))[, year := 2011]),
                         x + y ~ year, value.var = 'layer')

write.csv( grid_popwgt.xyz, file = '~/Dropbox/Harvard/RFMeval_Local/HyADS_to_pm25/HyADS_grid/population/hyads_grid_population.csv')


## =========================================================== ##
## #plot
## =========================================================== ##
ggplot( data = us_counties.pop06[state_name %in% states.lower]) + 
  facet_wrap( . ~ year, ncol = 1) +
  geom_sf( aes( fill = as.numeric( POP) / aland, geometry = geometry),
           color = NA) + 
  scale_fill_viridis( limits = c( 0, .0007), oob = squish)
ggplot( data = hyads_grid_popwgt06) + 
  # facet_wrap( . ~ year, ncol = 1) +
  geom_sf( aes( fill = POP, geometry = geometry),
           color = NA) + 
  scale_fill_viridis( limits = c( 0,400000))
