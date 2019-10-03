library( data.table)
library( raster)
library( sf)

## ========================================================= ##
##                Read in emissions data
## ========================================================= ##
d_cems_cmaq.f <- "~/Dropbox/Harvard/RFMeval_Local/CMAQ_DDM/COAL IMPACTS/INVENTORY/CEM/2005_cemsum.txt"
d_nonegu.f <- "~/Dropbox/Harvard/RFMeval_Local/CMAQ_DDM/COAL IMPACTS/INVENTORY/NONEGU COAL/ptinv_ptnonipm_xportfrac_cap2005v2_2005cs_orl_06jan2011_v4_orl_COAL.txt"

d_cmaq <- fread( d_cems_cmaq.f)
d_nonegu <- fread( d_nonegu.f, skip = "06029", header = F)[,1:63]
d_nonegu.names <- unlist( fread( d_nonegu.f, skip = 'FIPS,PLANTID,', header = F, nrows = 1))
names( d_nonegu) <- d_nonegu.names

d_nonegu.slim <- d_nonegu[ POLCODE == 'SO2', .( XLOC, YLOC, ANN_EMIS)]


## ========================================================= ##
##    Convert to spatial object, take over CMAQ raster
## ========================================================= ##
d_nonegu.sp <- st_as_sf( SpatialPointsDataFrame( d_nonegu.slim[, .( XLOC, YLOC)], 
                                       data.frame( d_nonegu.slim[, ANN_EMIS])))
plot( d_nonegu.sp)
