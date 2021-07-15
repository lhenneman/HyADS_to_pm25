library( sf)
library( raster)
library( data.table)
library( fst)
library( areal)
library( tigris)
library( ggplot2)
library( viridis)

## CHANGE THIS ##
# gridded hyads file location
hyads_file_loc <- '~/Dropbox/Harvard/Meetings_and_People/JoanCasey/EJ_PowerPlants/HyADS_census_tracts'

#coordinate reference system projection string for spatial data
p4s <- "+proj=lcc +lat_1=33 +lat_2=45 +lat_0=40 +lon_0=-97 +a=6370000 +b=6370000"

## ==================================================== ##
##          load Ohio census tracts
## ==================================================== ##
tracts.oh <- tracts( 'OH') %>% 
  st_as_sf() %>%
  st_transform( p4s)

# check out land area by tract
plot( tracts.oh[,'ALAND'], border = NA)

## ==================================================== ##
##          Define function
##          Read in hysplit raster, aggregate to areas
## ==================================================== ##
grids_to_counties <- function( n, files.in, years, area.sf){
  # get correct file and year
  file.in <- files.in[n]
  year.in <- years[n]
  
  print( file.in)
  # read in the file
  in.g <- read.fst( file.in, as.data.table = T)
  in.g[ is.na( in.g)] <- 0
  
  # remove units with all zeros
  in.unames <- names( in.g)[ !( names( in.g) %in% c( 'x', 'y'))]
  hyads.sums <- suppressWarnings( melt( in.g[, lapply( .SD, sum), .SDcols = ( in.unames)]))
  hyads.sums.use <- as.character( hyads.sums$variable)[hyads.sums$value > 0]
  
  # trim the dataset to include only units that are not 0
  in.g.trim <- in.g[, c( 'x', 'y', hyads.sums.use), with = F]
  
  # convert to raster
  in.r <- rasterFromXYZ( in.g.trim)
  crs( in.r) <- p4s
  names( in.r) <- names(in.g.trim)[!( names( in.g.trim) %in% c( 'x', 'y'))]
  if( length( names( in.r)) == 1)
    names( in.r) <- 'hyads'
  
  # convert to sf object
  ncin_spatpoly <- rasterToPolygons( in.r)
  ncin_sf <- st_as_sf( ncin_spatpoly)
  ncin_sf <- st_transform( ncin_sf, crs( area.sf))
  ncin_sf$GID <- 1:nrow( ncin_sf)
  
  # convert to data table and melt
  ncin.dt <- data.table( ncin_sf)[, geometry := NULL]
  ncin.m <- melt( ncin.dt, id.vars = 'GID',
                  variable.name = 'uID', value.name = 'pm25')
  
  # define small in variable
  ncin.train <- ncin_sf[,c( 'GID')]
  
  # take areas of zips over grids
  weights <- aw_intersect( area.sf, source = ncin.train, areaVar = "area")
  weights.dt <- data.table( weights)
  
  # take total areas of zip fores
  area.sf.a <- data.table( GEOID = area.sf$GEOID,
                           area.area = as.vector( st_area( area.sf)))
  
  # merge together for weighting dataset
  weights.m <- merge( weights.dt, area.sf.a, by = 'GEOID')
  weights.m[, areaWeight := area / area.area]
  
  # intersect weights and grids
  ncin.m.intersect <- merge( weights.m, ncin.m, by = 'GID', allow.cartesian = TRUE)
  
  # multiple weights by grid pm25, sum by 
  ncin.m.intersect[, pm25a := areaWeight * pm25]
  areas.pm <- ncin.m.intersect[, .( pm25 = sum( pm25a)), by = .( GEOID, uID)]
  
  # cast for smaller file size
  areas.pm.c <- dcast( areas.pm, GEOID ~ uID, value.var = 'pm25')
  
  # assign year column
  areas.pm.c[, year := year.in]
  
  return( areas.pm.c)
}

## ==================================================== ##
##          Run the function
## ==================================================== ##
# get the names of the gridded HyADS output files
grid.files.yr <- list.files( hyads_file_loc,
                             pattern = 'grids_.*\\d{4}\\.fst',
                             full.names = TRUE)

# fun function above to aggregate gridded values to census tracts each year
tracts_oh.hyads <- lapply( seq_along( grid.files.yr), 
                           grids_to_counties,
                           files.in = grid.files.yr,
                           years = 1999:2018,
                           area.sf = tracts.oh) %>% rbindlist

## ==================================================== ##
##          Plot the results
## ==================================================== ##
# merge back with spatial info
tracts_oh.hyads.sf <- merge( tracts.oh, tracts_oh.hyads, by = 'GEOID')

# quick plot of every four years (takes a few minutes to plot)
plotyears <- seq( 1999, 2018, 4)
map_hyads <- ggplot( tracts_oh.hyads.sf[tracts_oh.hyads.sf$year %in% plotyears,],
        aes( fill = hyads)) + 
  geom_sf( color = NA) + 
  coord_sf() +
  scale_fill_viridis() + 
  facet_wrap( . ~ year) + 
  theme_bw()

# save the plot
ggsave( file.path( hyads_file_loc, "oh_tracts_pm25_total_1999-2018_20200602.png"), map_hyads,
        height = 8, width = 8, scale = 1)

## ==================================================== ##
##          save the results
## ==================================================== ##
write.fst( tracts_oh.hyads, 
           path = file.path( hyads_file_loc, "oh_tracts_pm25_total_1999-2018_20200602.fst"))


