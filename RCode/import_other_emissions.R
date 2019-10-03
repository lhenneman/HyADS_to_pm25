library( data.table)
library( raster)

#coordinate reference system projection string for spatial data
p4s <- "+proj=lcc +lat_1=33 +lat_2=45 +lat_0=40 +lon_0=-97 +a=6370000 +b=6370000"

d_nonegu.f <- "~/Dropbox/Harvard/RFMeval_Local/CMAQ_DDM/COAL IMPACTS/INVENTORY/NONEGU COAL/ptinv_ptnonipm_xportfrac_cap2005v2_2005cs_orl_06jan2011_v4_orl_COAL.txt"

d_nonegu <- fread( d_nonegu.f, skip = "06029", header = F)[,1:63]
d_nonegu.names <- unlist( fread( d_nonegu.f, skip = 'FIPS,PLANTID,', header = F, nrows = 1))
names( d_nonegu) <- d_nonegu.names
d_nonegu.slim <- d_nonegu[ POLCODE == 'SO2', .( XLOC, YLOC, ANN_EMIS)]

## Convert to spatial object
d_nonegu.sp <- SpatialPointsDataFrame( d_nonegu.slim[, .( XLOC, YLOC)], 
                                       data.frame( d_nonegu.slim[, ANN_EMIS]),
                                       proj4string = CRS( "+proj=longlat +datum=WGS84 +no_defs"))
d_nonegu.sp <- spTransform( d_nonegu.sp, CRS( p4s))

# as raster
e <- extent( -2736000, 2592000, -2088000, 1944000)
r <- raster( nrows = 112, ncols = 148, ext = e, crs = crs( p4s))
d_nonegu.r <- rasterize( d_nonegu.sp, r)$d_nonegu.slim...ANN_EMIS.
d_nonegu.r[is.na(d_nonegu.r[])] <- 0

# plot
plot( d_nonegu.r)

# over USA
usa <- rnaturalearth::ne_countries(scale = 110, type = "countries", country = "United States of America", 
                                   geounit = NULL, sovereignty = NULL,
                                   returnclass = c("sp"))
usa.sub <- disaggregate(usa)[6,]
usa.sub <- spTransform(usa.sub, CRSobj = crs( p4s)) #proj4string( ))
plot( usa.sub, add = T)

# crop over USA
d_nonegu.r.mask <- mask( d_nonegu.r, usa.sub)
d_nonegu.r.crop <- crop( d_nonegu.r.mask, usa.sub)
plot( d_nonegu.r.crop)

