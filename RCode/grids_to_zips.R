library( sf)
library( raster)
library( data.table)
library( fst)

# read zcta shapefile and crosswalk
zip_sf_reader <- function( d = direct.dat){
  # zcta file downloaded from 'ftp://ftp2.census.gov/geo/tiger/GENZ2016/shp/cb_2016_us_zcta510_500k.zip'
  zcta_shapefile <- file.path( d, 'cb_2017_us_zcta510_500k.shp')
  
  # zcta-ZIP crosswalk file downloaded from 'http://mcdc2.missouri.edu/data/corrlst/'
  cw <- disperseR::crosswalk
  # cw <- fread( crosswalk_csv)
  # make sure ZCTA's are 5 digits to merge on zcta ID
  # cw$ZCTA <- formatC( cw$ZCTA, width = 5, format = "d", flag = "0") 
  
  zips <- st_read(zcta_shapefile)
  setnames( zips, 'ZCTA5CE10', 'ZCTA')
  zips <- merge( zips, cw, by = "ZCTA", all = F, allow.cartesian = TRUE)
  # make sure ZIPs are 5 digits to merge on zcta ID
  zips$ZIP <- formatC( zips$ZIP, width = 5, format = "d", flag = "0")
  
  return( zips)
}

## ==================================================== ##
##          load counties - just first 30 in TX
## ==================================================== ##
direct.dat <- '/n/zigler_lab/lhenneman/diseperseR/main/input/zcta_500k'
zips <- zip_sf_reader()[, 'ZIP']

## ==================================================== ##
##          Read in hysplit raster
## ==================================================== ##
in.f <- '/n/zigler_lab/lhenneman/diseperseR/main/output/exp25/grids_pm25_byunit_2000_01.fst'
in.g <- read.fst( in.f)
in.r <- rasterFromXYZ( in.g)

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





