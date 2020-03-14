#srun -p test --mem 21g -t 0-06:00 -c 1 -N 1 --pty /bin/bash
rm( list = ls())

do.raw <- TRUE
do.annual <- FALSE
do.xb <- FALSE

#coordinate reference system projection string for spatial data
p4s <- "+proj=lcc +lat_1=33 +lat_2=45 +lat_0=40 +lon_0=-97 +a=6370000 +b=6370000"

# load import functions/datasets
source( '~/repos/HyADS_to_pm25/RCode/hyads_to_pm25_functions.R')
load( '/n/zigler_lab/lhenneman/HyADS_to_pm25/rdata/hyads_to_cmaq_models.RData')

# point to hyads exposures (raw, unitless)
hyads.dir <- '/n/zigler_lab/lhenneman/diseperseR/main/output/exp'
hyadsPM25.dir <- '/n/zigler_lab/lhenneman/diseperseR/main/output/exp25'

# get usa mask for masking
# download USA polygon from rnaturalearth
us_states.names <- state.abb[!(state.abb %in% c( 'HI', 'AK'))]
us_states <- st_transform( USAboundaries::us_states(), p4s)
mask.usa <- sf::as_Spatial(us_states)[ us_states$state_abbr %in% us_states.names,]


#======================================================================#
## Prepare population data for functions below
#======================================================================#
mapply( hyads_to_pm25_unit, rep( 2010:2015, each = 12), 1:12,
        MoreArgs = list(
          fstart = file.path( hyads.dir, 'grids_exposures_byunit_'),
          fstart_out = file.path( hyadsPM25.dir, 'grids_pm25_byunit_'),
          model.dataset = preds.mon.hyads06w05,
          model.name = 'model.cv', #'model.gam'
          name.x = 'hyads',
          mask.use = mask.usa
        )
)



