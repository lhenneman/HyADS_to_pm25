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
hyads2005.dt <- fread( '~/Dropbox/Harvard/RFMeval_Local/HyADS_to_pm25/HyADS_grid/2005grid/HyADS_grid_2005.csv', drop = 'V1')
hyads2005 <- rasterFromXYZ( hyads2005.dt[, .( x, y, hyads)], crs = p4s)

## =========================================================== ##
## # get everything on the same crs
## =========================================================== ##
us_counties06.sf <- st_transform( st_as_sf( us_counties.pop06), crs = p4s)
us_counties11.sf <- st_transform( st_as_sf( us_counties.pop11), crs = p4s)
hyads2005.sf <- st_as_sf( rasterToPolygons( hyads2005))

## =========================================================== ##
## # area weight to hyads grid
## =========================================================== ##
hyads_grid_popwgt06 <- st_interpolate_aw( us_counties06.sf['POP'], 
                                          hyads2005.sf, extensive = T)
hyads_grid_popwgt11 <- st_interpolate_aw( us_counties11.sf['POP'], 
                                          hyads2005.sf, extensive = T)

## =========================================================== ##
## # onvert to raster and save
## =========================================================== ##
hyads_grid06.r <- fasterize( hyads_grid_popwgt06, hyads2005, field = 'POP')
hyads_grid11.r <- fasterize( hyads_grid_popwgt11, hyads2005, field = 'POP')
hyads_grid.xyz <- rbind( data.table( rasterToPoints( hyads_grid06.r))[, year := 2006],
                         data.table( rasterToPoints( hyads_grid11.r))[, year := 2011])

setnames( hyads_grid.xyz, 'layer', 'population')
write.csv( hyads_grid.xyz, file = '~/Dropbox/Harvard/RFMeval_Local/HyADS_to_pm25/HyADS_grid/population/hyads_grid_population.csv')


## =========================================================== ##
## #plot
## =========================================================== ##
ggplot( data = us_counties.pop[state_name %in% states.lower]) + 
  facet_wrap( . ~ year, ncol = 1) +
  geom_sf( aes( fill = as.numeric( POP) / aland, geometry = geometry),
           color = NA) + 
  scale_fill_viridis( limits = c( 0, .0007), oob = squish)
ggplot( data = hyads_grid_popwgt06) + 
  # facet_wrap( . ~ year, ncol = 1) +
  geom_sf( aes( fill = POP, geometry = geometry),
           color = NA) + 
  scale_fill_viridis( limits = c( 0,400000))
