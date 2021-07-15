#srun -p test --mem 21g -t 0-06:00 -c 1 -N 1 --pty /bin/bash
rm( list = ls())

do.raw <- TRUE
do.annual <- FALSE
do.xb <- FALSE

#coordinate reference system projection string for spatial data
p4s <- "+proj=lcc +lat_1=33 +lat_2=45 +lat_0=40 +lon_0=-97 +a=6370000 +b=6370000"

# load import functions/datasets
source( '~/repos/HyADS_to_pm25/RCode/hyads_to_pm25_functions.R')
# source( '~/Dropbox/Harvard/RFMeval_Local/HyADS_to_pm25/RCode/hyads_to_pm25_functions.R')
load( '/n/zigler_lab/lhenneman/HyADS_to_pm25/rdata/hyads_to_cmaq_models.RData')
# load( '/Users/lhenneman/Dropbox/Harvard/RFMeval_Local/HyADS/disperseR_hyads/hyads_to_cmaq_models.RData')

# point to hyads exposures (raw, unitless)
hyads.dir <- '/n/zigler_lab/lhenneman/diseperseR/main/output/exp'
hyadsPM25.dir <- '/n/zigler_lab/lhenneman/diseperseR/main/output/exp25'
# hyads.dir <- '/Users/lhenneman/Dropbox/Harvard/RFMeval_Local/HyADS/disperseR_hyads'
# hyadsPM25.dir <- '~/Google Drive/Harvard/HyADS/'

# get usa mask for masking
# download USA polygon from rnaturalearth
us_states.names <- state.abb[!(state.abb %in% c( 'HI', 'AK'))]
us_states <- st_transform( USAboundaries::us_states(), p4s)
mask.usa <- sf::as_Spatial(us_states)[ us_states$state_abbr %in% us_states.names,]


#======================================================================#
## Monthly conversions
#======================================================================#
#mapply( hyads_to_pm25_unit, rep( 2013:2015, each = 12), 1:12,
mapply( hyads_to_pm25_unit, rep( 2017:2018, each = 2), c( 2,6),
        MoreArgs = list(
                fstart = file.path( hyads.dir, 'grids_exposures_byunit_'),
                fstart_out = file.path( hyadsPM25.dir, 'grids_pm25_byunit_'),
                model.dataset = preds.mon.hyads06w05,
                model.name = 'model.cv', #'model.gam'
                name.x = 'hyads',
                mask.use = mask.usa
        )
)

#======================================================================#
## Monthly conversions -- total
#======================================================================#
#mapply( hyads_to_pm25_unit, rep( 2013:2015, each = 12), 1:12,
mapply( hyads_to_pm25_unit, rep( 1999:2018, each = 12), 1:12,
        MoreArgs = list(
                fstart = file.path( hyads.dir, 'grids_exposures_total_'),
                fstart_out = file.path( hyadsPM25.dir, 'grids_pm25_total_'),
                model.dataset = preds.mon.hyads06w05,
                model.name = 'model.cv', #'model.gam'
                name.x = 'hyads',
                mask.use = mask.usa,
                total = T
        )
)

#======================================================================#
## Annual conversions
#======================================================================#
#mapply( hyads_to_pm25_unit, rep( 2013:2015, each = 12), 1:12,
lapply( 2017:2018,
        hyads_to_pm25_unit,
        fstart = file.path( hyads.dir, 'grids_exposures_byunit_'),
        fstart_out = file.path( hyadsPM25.dir, 'grids_pm25_byunit_'),
        model.dataset = preds.ann.hyads06w05,
        model.name = 'model.cv', #'model.gam'
        name.x = 'hyads',
        mask.use = mask.usa
)

#======================================================================#
## Annual conversions -- total
#======================================================================#
#mapply( hyads_to_pm25_unit, rep( 2013:2015, each = 12), 1:12,
lapply( c( 2017:2018),
        hyads_to_pm25_unit,
        fstart = file.path( hyads.dir, 'grids_exposures_total_'),
        fstart_out = file.path( hyadsPM25.dir, 'grids_pm25_total_'),
        model.dataset = preds.ann.hyads06w05,
        model.name = 'model.cv', #'model.gam'
        name.x = 'hyads',
        mask.use = mask.usa,
        total = T
)

#======================================================================#
## Annual total conversions (just need 2005, so can use predictions)
#======================================================================#
preds.r <- preds.ann.hyads06w05$Y.ho.hat.raster$y.hat.lm.cv
dat.coords <- coordinates( preds.r)
dat_2005.dt <- data.table( cbind( dat.coords, values( preds.r)))
setnames( dat_2005.dt, 'V3', 'hyads_pm25')

write.fst( na.omit( dat_2005.dt), 
           path = "/n/zigler_lab/lhenneman/diseperseR/main/output/exp25/grids_pm25_total_2005.fst")


