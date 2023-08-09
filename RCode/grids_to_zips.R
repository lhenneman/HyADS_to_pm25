library( sf)
library( raster)
library( data.table)
library( fst)
library( areal)

#coordinate reference system projection string for spatial data
p4s <- "+proj=lcc +lat_1=33 +lat_2=45 +lat_0=40 +lon_0=-97 +a=6370000 +b=6370000"

# read zcta shapefile and crosswalk
zip_sf_reader <- function( d = direct.dat){
  # zcta file downloaded from 'ftp://ftp2.census.gov/geo/tiger/GENZ2016/shp/cb_2016_us_zcta510_500k.zip'
  zcta_shapefile <- file.path( d, 'cb_2017_us_zcta510_500k.shp')
  # zcta_shapefile <- file.path( d, 'cb_2016_us_zcta510_500k.shp')
  
  # zcta-ZIP crosswalk file downloaded from 'http://mcdc2.missouri.edu/data/corrlst/'
  cw <- disperseR::crosswalk
  # cw <- fread( crosswalk_csv)
  # make sure ZCTA's are 5 digits to merge on zcta ID
  # cw$ZCTA <- formatC( cw$ZCTA, width = 5, format = "d", flag = "0") 
  
  zips <- st_read(zcta_shapefile)
  zips <- st_transform( zips, p4s)
  setnames( zips, 'ZCTA5CE10', 'ZCTA')
  zips <- merge( zips, cw, by = "ZCTA", all = F, allow.cartesian = TRUE)
  # make sure ZIPs are 5 digits to merge on zcta ID
  zips$ZIP <- formatC( zips$ZIP, width = 5, format = "d", flag = "0")
  
  return( zips)
}

po_box <- st_read('~//Dropbox/GeorgeMason/Research/Data/GIS/USA_Zip_Code_Points/v10/zip_usa.gdb')

## ==================================================== ##
##          load counties - just first 30 in TX
## ==================================================== ##
direct.dat <- '/n/zigler_lab/lhenneman/diseperseR/main/input/zcta_500k'
# direct.dat <- '~/Dropbox/Harvard/Manuscripts/Energy_Transitions/Data_and_code/data/gis'
zips <- zip_sf_reader()[, 'ZIP']

## ==================================================== ##
##          Read in hysplit raster
## ==================================================== ##
grids_to_zips <- function( file.in){
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
  
  # convert to sf object
  ncin_spatpoly <- rasterToPolygons( in.r)
  ncin_sf <- st_as_sf( ncin_spatpoly)
  ncin_sf <- st_transform( ncin_sf, crs( zips))
  ncin_sf$GID <- 1:nrow( ncin_sf)
  
  # convert to data table and melt
  ncin.dt <- data.table( ncin_sf)[, geometry := NULL]
  ncin.m <- melt( ncin.dt, id.vars = 'GID',
                  variable.name = 'uID', value.name = 'pm25')
  
  # define small in variable
  ncin.train <- ncin_sf[,c( 'GID')]
  
  # take areas of zips over grids
  weights <- aw_intersect( zips, source = ncin.train, areaVar = "area")
  weights.dt <- data.table( weights)
  
  # take total areas of zip fores
  zips.a <- data.table( ZIP = zips$ZIP,
                        ZIP.area = as.vector( st_area( zips)))
  
  # merge together for weighting dataset
  weights.m <- merge( weights.dt, zips.a, by = 'ZIP')
  weights.m[, areaWeight := area / ZIP.area]
  
  # intersect weights and grids
  ncin.m.intersect <- merge( weights.m, ncin.m, by = 'GID', allow.cartesian = TRUE)
  
  # multiple weights by grid pm25, sum by 
  ncin.m.intersect[, pm25a := areaWeight * pm25]
  zips.pm <- ncin.m.intersect[, .( pm25 = sum( pm25a)), by = .( ZIP, uID)]
  
  # cast for smaller file size
  zips.pm.c <- dcast( zips.pm, ZIP ~ uID, value.var = 'pm25')
  
  # define output file, write
  file.out <- gsub( 'grids_', 'zips_', file.in)
  write.fst( zips.pm.c, path = file.out)
  
  return( file.out)
}

## ==================================================== ##
##          Run the function
## ==================================================== ##
grid.files <- list.files( '/n/zigler_lab/lhenneman/diseperseR/main/output/exp25/',
                          pattern = 'grids_.*_\\d{2}\\.fst',
                          full.names = TRUE)
grids_pm.list <- lapply( grid.files[c( 1:12, 205:240)], grids_to_zips)

grid.files.yr <- list.files( '/n/zigler_lab/lhenneman/diseperseR/main/output/exp25/',
                             pattern = 'grids_.*\\d{4}\\.fst',
                             full.names = TRUE)
grids_pm.list <- lapply( grid.files.yr[c( 1, 18:21, 38:40)], grids_to_zips)

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


## ==================================================== ##
##          Read in hysplit raster
## ==================================================== ##
in.f <- '/n/zigler_lab/lhenneman/diseperseR/main/output/exp25/grids_pm25_byunit_2000_01.fst'
in.f <- '~/Dropbox/Harvard/RFMeval_Local/HyADS_to_pm25/HyADS_grid/grid_to_zips/grids_pm25_byunit_2000_01.fst'
in.g <- read.fst( in.f, as.data.table = T)
in.g[ is.na( in.g)] <- 0
in.r <- rasterFromXYZ( in.g)
crs( in.r) <- p4s
names( in.r) <- names(in.g)[!( names(in.g) %in% c( 'x', 'y'))]

ncin_spatpoly <- rasterToPolygons( in.r)
ncin_sf <- st_as_sf( ncin_spatpoly)
ncin_sf <- st_transform( ncin_sf, crs( zips))
ncin_sf$GID <- 1:nrow( ncin_sf)


## ==================================================== ##
##          Exploring
## ==================================================== ##
summary( ncin_sf[,1:5])

# scales linearly by time

system.time(
  in.zips1 <- aw_interpolate( zips, tid = ZIP, sid = GID, weight = 'sum', output = 'sf',
                              source = ncin_sf, intensive = names( ncin_sf)[1])
)

aw_preview_weights(zips, tid = ZIP, source = ncin_sf, sid = GID, 
                   type = "intensive")

system.time(
  in.zips5 <- aw_interpolate( zips, tid = ZIP, sid = GID, weight = 'sum', output = 'sf',
                              source = ncin_sf, intensive = names( ncin_sf)[1:5])
)
system.time(
  in.zips50 <- aw_interpolate( zips, tid = ZIP, sid = GID, weight = 'sum', output = 'sf',
                               source = ncin_sf, intensive = names( ncin_sf)[1:50])
)


in.zips <- ar_validate( source = ncin_sf, target = zips, 
                        varList = names( ncin_sf)[1:5], verbose = TRUE)

## ==================================================== ##
##          Function for the manual way
## ==================================================== ##
# sum of 'area' variable matches actual area of ZIP
# 

ncin.train <- ncin_sf[,c( 'GID')]
ncin.train

weights.pv <- aw_preview_weights(zips, tid = ZIP, source = ncin.train, sid = GID, 
                                 type = "intensive")

weights.ex <-  zips %>%
  aw_intersect(source = ncin.train, areaVar = "area") %>%
  aw_total(source = ncin.train, id = GID, areaVar = "area", totalVar = "totalArea",
           type = "extensive", weight = "sum") %>%
  aw_weight(areaVar = "area", totalVar = "totalArea", 
            areaWeight = "areaWeight")
weights.in <-  zips %>%
  aw_intersect(source = ncin.train, areaVar = "area") %>%
  aw_total(source = ncin.train, id = GID, areaVar = "area", totalVar = "totalArea",
           type = "intensive", weight = "sum") %>%
  aw_weight(areaVar = "area", totalVar = "totalArea", 
            areaWeight = "areaWeight")

int.86047 <- st_intersection( zips[zips$ZIP == '86047',], ncin.train)
int.86047a <- st_area( int.86047)

weights <-  zips %>%
  aw_intersect(source = ncin.train, areaVar = "area") %>% data.table
zips.a <- data.table( ZIP = zips$ZIP,
                      ZIP.area = as.vector( st_area( zips)))

weights.m <- merge( weights, zips.a, by = 'ZIP')
weights.m[, areaWeight := area / ZIP.area]
weights.m[ZIP == '86047']

weights.ex.dt <- data.table( weights.ex)#[, .(GID, ZIP, areaWeight)]
weights.in.dt <- data.table( weights.in)#[, .(GID, ZIP, areaWeight)]
weights.ex.dt[ZIP == '86047']
weights.in.dt[ZIP == '86047']


cols.n <- 1000
ncin.dt <- data.table( ncin_sf[,c( names( ncin_sf)[1:cols.n], 'GID')])[, geometry := NULL]
ncin.m <- melt( ncin.dt, id.vars = 'GID',
                variable.name = 'uID', value.name = 'pm25')

ncin.m.intersect <- merge( weights.m, ncin.m, by = 'GID', allow.cartesian = TRUE)

ncin.m.intersect[, pm25a := areaWeight * pm25]
zips.pm <- ncin.m.intersect[, .( pm25 = sum( pm25a)), by = .( ZIP, uID)]

zips.pm[pm25 > .2]
in.zips1 <- aw_interpolate( zips, tid = ZIP, sid = GID, weight = 'sum', output = 'sf',
                            source = ncin_sf, intensive = c( 'X113.3'))
data.table( in.zips1)[ZIP == '86047']
ncin.m.intersect[ZIP == '86047' & uID == 'X113.3']
weights.dt[ZIP == '86047']

zips.area <- ncin.m.intersect[, .( area = sum( areaWeight)), by = .( ZIP, uID)]
wgts.area <- weights.dt[, .( area = sum( areaWeight)), by = .( ZIP)]

ncin.m.intersect

## ==================================================== ##
##          Try the manual way
## ==================================================== ##
ncin.train <- ncin_sf[,c( 'X10.2', 'GID')]

zips %>%
  aw_intersect(source = ncin.train, areaVar = "area") -> intersect

data.table( intersect)[ZIP == '01002']

intersect %>%
  aw_total(source = ncin.train, id = GID, areaVar = "area", totalVar = "totalArea",
           type = "intensive", weight = "sum") -> intersect

intersect %>%
  aw_weight(areaVar = "area", totalVar = "totalArea", 
            areaWeight = "areaWeight") -> intersect

data.table( intersect)[ZIP == '39153']

intersect %>%
  aw_calculate(value = 'X10.2', areaWeight = "areaWeight") -> intersect



intersect %>%
  aw_aggregate(target = zips, tid = ZIP, interVar = X10.2) -> result

result.dt <- data.table( result)
result.dt[ZIP == '39153']


ncin.dt <- data.table( ncin_sf[,c( names( ncin_sf)[1:5], 'GID')])[, geometry := NULL]
ncin.m <- melt( ncin.dt, id.vars = 'GID',
                variable.name = 'uID', value.name = 'pm25')

ncin.m.intersect <- merge( intersect[, .( GID, ZIP, areaWeight)], 
                           ncin.m, by = 'GID', allow.cartesian = TRUE)

ncin.m.intersect[, pm25a := areaWeight * pm25]
zips.pm <- ncin.m.intersect[, .( pm25 = sum( pm25a)), by = .( ZIP, uID)]

ncin.m.intersect[ZIP == '39153' & uID == 'X10.2']
zips.pm[ZIP == '39153']
## ==================================================== ##
##          Define function for overing counties
## ==================================================== ##
species_over_zips <- function( file.in,
                               species,
                               .zips = zips){
  print( file.in)
  
  ## define file name and open ncdf file
  ncin_raster <- raster( file.in)
  
  ## plot, check out projection
  # plot( ncin_raster)
  
  ## read in zips
  zips.sp <- sf::as_Spatial( .zips)
  zips.sp <- spTransform(zips.sp, proj4string( ncin_raster))
  # zips.sf <- st_transform(.zips, p4s)
  
  # limit to lower 48
  states.use <- state.abb[!(state.abb %in% c( 'HI', 'AK'))]
  
  # run one state at a time
  out.dt <- data.table()
  for( s in states.use){
    print( s)
    zips.sp.use <- zips.sp[ zips.sp$STATE == s,]
    zips.sf.use <- .zips[ .zips$STATE == s,]
    
    ## trim ncin_raster
    extent.use <- extent( extent( zips.sp.use)@xmin - .1,
                          extent( zips.sp.use)@xmax + .1,
                          extent( zips.sp.use)@ymin - .1,
                          extent( zips.sp.use)@ymax + .1)
    ncin_raster.crop <- crop( ncin_raster, extent.use)
    
    # convert to spatial data frame if you want, change the projection, etc...
    ncin_spatpoly <- rasterToPolygons( ncin_raster.crop)
    ncin_sf <- st_as_sf( ncin_spatpoly)
    ncin_sf <- st_transform(ncin_sf, st_crs( .zips))
    
    # get county average concentrations
    # zips.oversp <- over( zips.sp.use[1:100,], ncin_spatpoly)
    zips.oversf <- st_interpolate_aw( ncin_sf, zips.sf.use, extensive = F)
    zips.out <- data.table( as.data.table( zips.oversf), 
                            ZIP = zips.sf.use$ZIP, STATE = zips.sf.use$STATE)
    
    # house cleaning, labeling, etc
    .year <- as( gsub( '^.*NA_|\\d{2}_.*', '', file.in), 'numeric')
    .Month <- gsub( '^.*NA_\\d{4}|_.*', '', file.in) 
    .month <-as( gsub( '^.*NA_\\d{4}|_.*', '', file.in), 'numeric')
    zips.out[, `:=` ( year = .year,
                      month = .month,
                      species = species,
                      Group.1 = NULL,
                      geometry = NULL)]
    
    out.dt <- rbind( out.dt, zips.out)
  }
  fname <- paste0( '~/Dropbox/Harvard/RFMeval_Local/PollutantFields_RM/pm25_zips2/so4_zips', .year, '_', .Month, '.csv')
  # fname <- paste0( '/n/scratchlfs/zigler_lab/lhenneman/pm25_zip_link/pm25_zips2/pm25_zips', .year, '_', .Month, '.csv')
  write.csv( file = fname, out.dt)
  
  return( out.dt)
}

## ==================================================== ##
##          Run function for overing counties
## ==================================================== ##
over_counties.pm25 <- rbindlist( lapply( filesall05.so4,
                                         species_over_zips,
                                         species = 'so4'))





