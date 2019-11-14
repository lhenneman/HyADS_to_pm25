library( disperseR)
library( raster)
library( data.table)

load( '~/Dropbox/Harvard/RFMeval_Local/HyADS/hyads_unwgted_nopbl/hyspdisp_grid_2005.RData')

june <- data.table( MAP8.2005)

dims <- 3:50 #dim( june)[2]

june.sum <- rowSums( june[, ..dims], na.rm = T)

june2 <- cbind( june[,.(x,y)], june.sum)
june_10017.2 <- cbind( june[,.(x,y)], june$`10017.2`)


june.r <- rasterFromXYZ( june2)
june_10017.2.r <- rasterFromXYZ( june_10017.2)
plot( june.r)
plot( june_10017.2.r)


# check out the boundary layer height
Sys.setenv(TZ='UTC')
pbl <- brick('~/Desktop/run_hyspdisp/hpbl/hpbl.mon.mean.nc')
plot( pbl$X2005.06.01.00.56.02)

# usa sub
usa <- rnaturalearth::ne_countries(scale = 110, type = "countries", country = "United States of America",
                                   geounit = NULL, sovereignty = NULL,
                                   returnclass = c("sp"))
usa.sub <- disaggregate(usa)[6,]
crop.extent <- usa.sub
p4s = "+proj=aea +lat_1=20 +lat_2=60 +lat_0=40 +lon_0=-96 +x_0=0 +y_0=0 +ellps=GRS80 +datum=NAD83 +units=m"
crop.extent.proj <- projectExtent( crop.extent, p4s)

# check out the raw gridded files
zpc_dir <- '~/Dropbox/Harvard/RFMeval_Local/HyADS/hyads_unwgted/june_2005_test'

ym <- '20056'

year.h <- substr( ym, 1, 4)
month.m <- as.integer( substr( ym, 5, 6))
month.h <- formatC( month.m, width = 2, format = "d", flag = "0")
pattern <- paste0( 'gridlinks.*', year.h,'-', month.h)

files.month <- list.files( path = zpc_dir,
                           pattern = pattern,
                           full.names = T)

unitnames <- gsub( paste0( '.*gridlinks_|_', year.h, '-', month.h, '.*csv$'),
                   '', files.month)
names(files.month) <- unitnames

pb <- txtProgressBar(min = 0, max = length( files.month), style = 3)
data.h <- lapply( seq_along(files.month),
                  function( i,
                            files){
                    gc()
                    d <- fread( files[i],
                                drop = 'V1', showProgress = FALSE)
                    setnames( d, 'hyspdisp', names(files)[i])
                    # d[, `:=` (uID = names(files)[i] )]
                    
                    r <- rasterFromXYZ( d)
                    crs( r) <- p4s
                    rm( d)
                    if( !is.null( crop.extent))
                      r <- crop( r, crop.extent.proj)
                    
                    setTxtProgressBar(pb, i)
                    
                    return(r)
                  },
                  files.month)


plot( data.h[[44]])
plot( data.h[[5]]) #######
plot( usa.sub, add = T)

project_and_stack <- function( ..., mask.use = NULL){
  if( is.list( ...)){
    list.r <- unlist( ...)
  } else{
    list.r <- list( ...)
  }
  
  for( r in 2: length( list.r)){
    list.r[[r]] <- projectRaster( list.r[[r]], list.r[[1]])
  }
  
  # mask over usa
  if( !is.null( mask.use)){
    list.out <- lapply( list.r, mask, mask.use)
  } else
    list.out <- list.r
  
  return( stack( list.out))
}

data.h.stack <- project_and_stack(data.h)
plot( sum( data.h.stack, na.rm = T))
plot( sum( data.h.stack[[c( 5, 44:46)]], na.rm = T)) #44:46

usa.sub.t <- spTransform( usa.sub, p4s)
plot( usa.sub.t, add = T)

units <- c( '10017.2', '10143.AAB01',   '10151.1A',   '10151.1B' )

## =================================================== ##
## do the linking for unit 5

hyo_dir <- '~/Dropbox/Harvard/RFMeval_Local/HyADS/rawHYSPLIT'
zpc_dir <- '~/Dropbox/Harvard/RFMeval_Local/HyADS/linkHYSPLIT'
month_YYYYMM <- '20056'
start_date <- NULL
end_date <- NULL
unit <- units2005[uID == '10017.2']
duration_run_hours <- 240
hpbl_raster <- pbl
overwrite <- F

if( is.null( start_date) | is.null( end_date)){
  start_date <- as.Date( paste( substr( month_YYYYMM, 1, 4),
                                substr( month_YYYYMM, 5, 6),
                                '01', sep = '-'))
  
  end_date <- seq( start_date,
                   by = paste (1, "months"),
                   length = 2)[2] - 1
}

## name the eventual output file
grid_output_file <- file.path( zpc_dir,
                               paste0("gridlinks_",
                                      unit$ID, "_",
                                      start_date, "_",
                                      end_date,
                                      ".csv"))

vec_dates <- seq.Date( as.Date( start_date),
                       as.Date( end_date),
                       by = '1 day')
vec_filedates <- seq.Date( from = as.Date( start_date) - ceiling( duration_run_hours / 24),
                           to = as.Date( end_date),
                           by = '1 day')

## list the files
pattern.file <- paste0( '_', gsub( '[*]', '[*]', unit$ID), '_(', paste(vec_filedates, collapse = '|'), ')')
files.read <- list.files(path = hyo_dir,
                         pattern = pattern.file,
                         full.names = T)

# read the files
l <- lapply(files.read,
            fread)

## Combine all parcels into single data table
d <- rbindlist(l)

## Trim dates & first hour
d <- d[d$Pdate %in% as( c( vec_dates), "character") &
         hour > 1,]

ggplot( data = d) +
  geom_hex( aes( x = lon, y = lat), bins = 1000) + 
  # scale_fill_manual( limits = c( 1000, 5000)) +
  coord_cartesian( xlim = c( -80, -60), ylim = c(35, 45))

# check extent
d_xmin <- min( d$lon)
e_xmin <- extent( hpbl_raster)[1]
if( d_xmin < e_xmin - 5)
  hpbl_raster <- rotate( hpbl_raster)

# trim pbl
d_trim <- trim_pbl( d,
                    rasterin = hpbl_raster)
d_notrim <- d
ggplot( data = d_trim) +
  geom_hex( aes( x = lon, y = lat), bins = 1000) + 
  # scale_fill_manual( limits = c( 1000, 5000)) +
  coord_cartesian( xlim = c( -80, -60), ylim = c(35, 45))

## Link to grid
disp_df_link <- link_zip( d. = d_notrim,
                          p4string = p4s,
                          rasterin = hpbl_raster,
                          return.grid = T)

## deeper into zip link
xy <- d.[,.( lon, lat)]
spdf.in <- SpatialPointsDataFrame( coords = xy,
                                   data = d.,
                                   proj4string = CRS( "+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"))
spdf <- spTransform( spdf.in, p4string)

spdf.in.sf <- st_as_sf( spdf.in)

rasterize( spdf.in)

r <- rasterFromXYZ( disp_df_link)
plot( r, ylim = c( -2e6, 2e6), xlim = c(1e6, 4e6))

# extract data layer from raster, disaggregate to .1°x.1°
pbl_layer <- subset_nc_date(hpbl_brick = rasterin,
                            vardate = d$Pdate[1])
pbl_layer.t <- projectRaster( pbl_layer,
                              crs = CRS( proj4string( spdf)))

# aim for a resolution of 12 km
pbl_resolution <- res( pbl_layer.t)
x_fact <- floor( pbl_resolution[1] / 12000)
y_fact <- floor( pbl_resolution[2] / 12000)
pbl_layer.d <- disaggregate( pbl_layer.t,
                             fact = c( x_fact, y_fact))

# count number of particles in each cell,
# find original raster cells, divide number by pbl
r <- pbl_layer.d
values( r) <- NA
cells <- cellFromXY( r, spdf)
tab <- table( cells)
pbls <- pbl_layer.d[as.numeric( names( tab))]
r[as.numeric( names( tab))] <- tab #/ pbls

# crop around point locations for faster extracting
e <- extent(spdf)
r2 <- crop( trim(r,
                 padding = 1),
            e)
plot( r2, ylim = c( -1e6, 1e6), xlim = c(1e6, 3e6))
plot( pbl_layer.d, ylim = c( -1e6, 1e6), xlim = c(1e6, 3e6))
