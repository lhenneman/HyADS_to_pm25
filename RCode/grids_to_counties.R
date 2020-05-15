library( sf)
library( raster)
library( data.table)
library( fst)
library( areal)

#coordinate reference system projection string for spatial data
p4s <- "+proj=lcc +lat_1=33 +lat_2=45 +lat_0=40 +lon_0=-97 +a=6370000 +b=6370000"


## ==================================================== ##
##          load counties - just first 30 in TX
## ==================================================== ##
states_48 <- state.name[!(state.name %in% c( 'Hawaii', 'Alaska'))]
counties.us <- USAboundaries::us_counties()
counties.us <- st_transform( counties.us, p4s)[counties.us$state_name %in% states_48,]


## ==================================================== ##
##          Read in hysplit raster
## ==================================================== ##
grids_to_counties <- function( file.in){
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
  ncin_sf <- st_transform( ncin_sf, crs( counties.us))
  ncin_sf$GID <- 1:nrow( ncin_sf)
  
  # convert to data table and melt
  ncin.dt <- data.table( ncin_sf)[, geometry := NULL]
  ncin.m <- melt( ncin.dt, id.vars = 'GID',
                  variable.name = 'uID', value.name = 'pm25')
  
  # define small in variable
  ncin.train <- ncin_sf[,c( 'GID')]
  
  # take areas of zips over grids
  weights <- aw_intersect( counties.us, source = ncin.train, areaVar = "area")
  weights.dt <- data.table( weights)
  
  # take total areas of zip fores
  counties.us.a <- data.table( geoid = counties.us$geoid,
                        county.area = as.vector( st_area( counties.us)))
  
  # merge together for weighting dataset
  weights.m <- merge( weights.dt, counties.us.a, by = 'geoid')
  weights.m[, areaWeight := area / county.area]
  
  # intersect weights and grids
  ncin.m.intersect <- merge( weights.m, ncin.m, by = 'GID', allow.cartesian = TRUE)
  
  # multiple weights by grid pm25, sum by 
  ncin.m.intersect[, pm25a := areaWeight * pm25]
  counties.pm <- ncin.m.intersect[, .( pm25 = sum( pm25a)), by = .( geoid, uID)]
  
  # cast for smaller file size
  counties.pm.c <- dcast( counties.pm, geoid ~ uID, value.var = 'pm25')
  
  # define output file, write
  file.out <- gsub( 'grids_', 'counties_', file.in)
  write.fst( counties.pm.c, path = file.out)
  
  return( file.out)
}

## ==================================================== ##
##          Run the function
## ==================================================== ##
# grid.files <- list.files( '/n/zigler_lab/lhenneman/diseperseR/main/output/exp25/',
#                           pattern = 'grids_.*_\\d{2}\\.fst',
#                           full.names = TRUE)
# grids_pm.list <- lapply( grid.files[181:192], grids_to_zips)

grid.files.yr <- list.files( '/n/zigler_lab/lhenneman/diseperseR/main/output/exp25/',
                             pattern = 'grids_.*\\d{4}\\.fst',
                             full.names = TRUE)
grids_pm.list <- lapply( grid.files.yr, grids_to_counties)

# gather the total county files into single file
count_files_total <- list.files( '/n/zigler_lab/lhenneman/diseperseR/main/output/exp25/',
                                 pattern = 'counties_pm25_total_\\d{4}\\.fst',
                                 full.names = TRUE)
county_allyears <- rbindlist( lapply( 1:20,
                                      function( n, f, y) {
                                        read_fst( f[n], as.data.table = T)[, year := y[n]]
                                      }, count_files_total, 1999:2018))

write.fst( county_allyears, 
           path = "/n/zigler_lab/lhenneman/diseperseR/main/output/exp25/counties_pm25_total_1999-2018.fst")

## ==================================================== ##
##          Run the function for raw hyads
## ==================================================== ##
grid.files2005 <- list.files( '/n/zigler_lab/lhenneman/diseperseR/main/output/exp/',
                              pattern = 'grids_.*total_2005.fst',
                              full.names = TRUE)
grids_pm.list <- lapply( grid.files2005, grids_to_zips)

grid.files2005 <- list.files( '/n/zigler_lab/lhenneman/diseperseR/main/output/exp25/',
                              pattern = 'grids_.*total_2005.fst',
                              full.names = TRUE)

grids_pm.list <- lapply( grid.files2005, grids_to_zips)

