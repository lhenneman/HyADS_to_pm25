#srun -p test --mem 100g -t 0-06:00 -c 1 -N 1 --pty /bin/bash
rm( list = ls())

#source( '~/Dropbox/Harvard/RFMeval_Local/HyADS_to_pm25/RCode/hyads_to_pm25_functions.R')
source( '~/repos/HyADS_to_pm25/RCode/hyads_to_pm25_functions.R')
#coordinate reference system projection string for spatial data
p4s <- "+proj=lcc +lat_1=33 +lat_2=45 +lat_0=40 +lon_0=-97 +a=6370000 +b=6370000"

#======================================================================#
## Load saved object of cmaq-ddm / hyads models from hyads_to_pm25_month.R
#======================================================================#
#load( '~/Dropbox/Harvard/RFMeval_Local/HyADS_to_pm25/RData/hyads_to_cmaq_models.RData')
load( '/n/scratchlfs/zigler_lab/lhenneman/run_hyspdisp/output_ampd_dists/hyads_to_cmaq_models.RData')
# load the models 
hyads.ann.model <- preds.ann.hyads06w05$model.gam
idwe.ann.model <- preds.ann.idwe06w05$model.gam

# set up the dataset
dat.coords <- coordinates( dats2006.a)
dats2006raw.dt <- data.table( cbind( dat.coords, values( dats2006.a)))

# do the predictions
hyads2006.pred <- predict( hyads.ann.model, newdata = dats2006raw.dt)
idwe2006.pred <- predict( idwe.ann.model, newdata = dats2006raw.dt)

# rasterize
dats2006.r <- rasterFromXYZ( data.table( dat.coords, hyads2006.pred, idwe2006.pred), crs = p4s)
dats2006.r[is.na( dats2006.r)] <- 0
#======================================================================#
## get total state populations
#======================================================================#
# grid_popwgt.xyz <- fread( '~/Dropbox/Harvard/RFMeval_Local/HyADS_to_pm25/HyADS_grid/population/hyads_grid_population.csv',
#                          drop = 'V1')

grid_popwgt.xyz <- fread( '/n/home03/lhenneman/inputdata/census_population/hyads_grid_population.csv',
                          drop = 'V1')
grid_popwgt.r <- rasterFromXYZ( grid_popwgt.xyz, crs = p4s)
grid_popwgt2006.r <- grid_popwgt.r$X2006

#ggplot.a.raster( unstack( grid_popwgt.r), bounds = c( 0, .5e6), facet.names = c( '2006', '2011'))

grid_popwgt.sf <- st_as_sf( rasterToPolygons( grid_popwgt.r))

# get total state populations
us_states <- st_transform( USAboundaries::us_states(), p4s)
us_states.p <- st_interpolate_aw( grid_popwgt.sf, us_states, extensive = T)

# merge back to state information
us_states.pop <- cbind( us_states.p, us_states[us_states.p$Group.1,])
us_states.pop$geometry.1 <- NULL
us_states.pop.dt <- data.table( us_states.pop)

#======================================================================#
## area weight over states
#======================================================================#
# combine exposure and population
dats2006_popwgt.r <- project_and_stack( dats2006.r, grid_popwgt.r)

# weight by 2006 population
dats2006_popwgt.names <- names( dats2006.r)
dats2006_popwgtexp.r <- subset( dats2006_popwgt.r, dats2006_popwgt.names)
dats2006_popwgt.r <- dats2006_popwgtexp.r * dats2006_popwgt.r$X2006
names( dats2006_popwgt.r) <- dats2006_popwgt.names

# take over states
dats2006_popwgt.sf <- st_as_sf( rasterToPolygons( dats2006_popwgt.r))

#NA's coming from 
dats2006_popwgt.states <- st_interpolate_aw( dats2006_popwgt.sf, us_states, extensive = F)

# get usa mask for masking
# download USA polygon from rnaturalearth
us_states.names <- state.abb[!(state.abb %in% c( 'HI', 'AK'))]
mask.usa <- sf::as_Spatial(us_states)[ us_states$state_abbr %in% us_states.names,]

#======================================================================#
## Do the conversions
#======================================================================#
states.use <- c( 'PA', 'KY', 'GA', 'WI', 'TX', 'CO', 'CA')
#fstart.idwe <- '~/Dropbox/Harvard/RFMeval_Local/HyADS_to_pm25/RData/ampd_dists_sox_weighted'
#fstart.hyads <- '~/Dropbox/Harvard/RFMeval_Local/HyADS_to_pm25/HyADS_grid/2006grid/GRIDexposures_byunit_'

fstart.idwe <- '/n/scratchlfs/zigler_lab/lhenneman/run_hyspdisp/output_ampd_dists/ampd_dists_sox_weighted'
fstart.hyads <- '/n/scratchlfs/zigler_lab/lhenneman/run_hyspdisp/output_HyADS/GRIDexposures_byunit_'

saveloc.idwe  <- '/n/scratchlfs/zigler_lab/lhenneman/run_hyspdisp/output_popweighted/popwgt_idwe'
saveloc.hyads <- '/n/scratchlfs/zigler_lab/lhenneman/run_hyspdisp/output_popweighted/popwgt_hyads'

array_num <- as.numeric( Sys.getenv("SLURM_ARRAY_TASK_ID"))
array_num <- ifelse( array_num == ''|is.na( array_num), 1, array_num)


# IDWE
if( array_num == 1){
idwe_exp06 <- rbindlist( lapply( 1:12, 
                                 state_exposurer,
                                 fstart = fstart.idwe,
                                 year.m = 2006,
                                 model.dataset = preds.mon.idwe06w05,
                                 model.name = 'model.cv', #'model.gam'
                                 name.x = 'idwe',
                                 mask.use = mask.usa)) #[ mask.usa$state_abbr %in% states.use,]))

  write.csv( file = paste0( saveloc.idwe, '2006.csv'), idwe_exp06)
}

if( array_num == 2){
idwe_exp11 <- rbindlist( lapply( 1:12, 
                                 state_exposurer,
                                 fstart = fstart.idwe,
                                 year.m = 2011,
                                 model.dataset = preds.mon.idwe06w05,
                                 model.name = 'model.cv', #'model.gam'
                                 name.x = 'idwe',
                                 mask.use = mask.usa))
  write.csv( file = paste0( saveloc.idwe, '2011.csv'), idwe_exp11)
}

if( array_num %in% 3:14){
mon <- array_num - 2
# HyADS
hyads_exp06 <- rbindlist( lapply( mon, 
                                  state_exposurer,
                                  fstart = fstart.hyads,
                                  year.m = 2006,
                                  model.dataset = preds.mon.hyads06w05,
                                  model.name = 'model.gam',
                                  name.x = 'hyads',
                                  mask.use = mask.usa))
  write.csv( file = paste0( saveloc.hyads, 2006, '_', mon, '.csv'), hyads_exp06)
}

if( array_num %in% 15:26){
mon <- array_num - 14
hyads_exp11 <- rbindlist( lapply( mon, 
                                  state_exposurer,
                                  fstart = fstart.hyads,
                                  year.m = 2011,
                                  model.dataset = preds.mon.hyads06w05,
                                  model.name = 'model.gam',
                                  name.x = 'hyads',
                                  mask.use = mask.usa))

  write.csv( file = paste0( saveloc.hyads, 2011, '_', mon, '.csv'), hyads_exp11)
}

