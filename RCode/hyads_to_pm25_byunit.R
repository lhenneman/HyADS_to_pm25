#srun -p test --mem 100g -t 0-06:00 -c 1 -N 1 --pty /bin/bash
rm( list = ls())

platform <- c( 'mac', 'cannon')[2]
do.annual <- FALSE

#coordinate reference system projection string for spatial data
p4s <- "+proj=lcc +lat_1=33 +lat_2=45 +lat_0=40 +lon_0=-97 +a=6370000 +b=6370000"

if( platform == 'mac'){
  source( '~/Dropbox/Harvard/RFMeval_Local/HyADS_to_pm25/RCode/hyads_to_pm25_functions.R')
  load( '~/Dropbox/Harvard/RFMeval_Local/HyADS_to_pm25/RData/hyads_to_cmaq_models.RData')
  grid_popwgt.xyz <- fread( '~/Dropbox/Harvard/RFMeval_Local/HyADS_to_pm25/HyADS_grid/population/hyads_grid_population.csv',
                            drop = 'V1')
  fstart.idwe <- '~/Dropbox/Harvard/RFMeval_Local/HyADS_to_pm25/RData/ampd_dists_sox_weighted'
  fstart.hyads <- '~/Dropbox/Harvard/RFMeval_Local/HyADS_to_pm25/HyADS_grid/2006grid/GRIDexposures_byunit_'
  
  fname2006.hyads <- '~/Dropbox/Harvard/RFMeval_Local/HyADS_to_pm25/HyADS_grid/2006grid/HyADSunit_grid_annual_nopbl_2006.csv'
  fname2011.hyads <- '~/Dropbox/Harvard/RFMeval_Local/HyADS_to_pm25/HyADS_grid/2011grid/HyADSunit_grid_annual_nopbl_2011.csv'
  fname2006.idwe <- '~/Dropbox/Harvard/RFMeval_Local/HyADS_to_pm25/RData/ampd_dists_sox_weighted_2006_unit.csv'
  fname2011.idwe <- '~/Dropbox/Harvard/RFMeval_Local/HyADS_to_pm25/RData/ampd_dists_sox_weighted_2011_unit.csv'
  
  saveloc.idwe  <- '~/Dropbox/Harvard/RFMeval_Local/HyADS_to_pm25/RData/popwgt_idwe'
  saveloc.hyads <- '~/Dropbox/Harvard/RFMeval_Local/HyADS_to_pm25/RData/popwgt_hyads'
} else{
  source( '~/repos/HyADS_to_pm25/RCode/hyads_to_pm25_functions.R')
  load( '/n/scratchlfs/zigler_lab/lhenneman/run_hyspdisp/output_ampd_dists/hyads_to_cmaq_models.RData')
  grid_popwgt.xyz <- fread( '/n/home03/lhenneman/inputdata/census_population/hyads_grid_population.csv',
                            drop = 'V1')
  fstart.idwe <- '/n/scratchlfs/zigler_lab/lhenneman/run_hyspdisp/output_ampd_dists/ampd_dists_sox_weighted'
  fstart.hyads <- '/n/scratchlfs/zigler_lab/lhenneman/run_hyspdisp/output_HyADS/GRIDexposures_byunit_'
  
  fname2006.hyads <- '/n/scratchlfs/zigler_lab/lhenneman/run_hyspdisp/output_HyADS/2006grid/HyADSunit_grid_annual_nopbl_2006.csv'
  fname2011.hyads <- '/n/scratchlfs/zigler_lab/lhenneman/run_hyspdisp/output_HyADS/2011grid/HyADSunit_grid_annual_nopbl_2011.csv'
  fname.ann.idwe <- '/n/scratchlfs/zigler_lab/lhenneman/run_hyspdisp/output_HyADS/2006grid/ampd_dists_sox_weighted_annual_unit.csv'
  
  saveloc.idwe  <- '/n/scratchlfs/zigler_lab/lhenneman/run_hyspdisp/output_popweighted/popwgt_idwe'
  saveloc.hyads <- '/n/scratchlfs/zigler_lab/lhenneman/run_hyspdisp/output_popweighted/popwgt_hyads'
  
}


#======================================================================#
## Prepare population data for functions below
#======================================================================#
# grid-weighted population as raster
grid_popwgt.r <- rasterFromXYZ( grid_popwgt.xyz, crs = p4s)
grid_popwgt.sf <- st_as_sf( rasterToPolygons( grid_popwgt.r))

# get total state populations
us_states <- st_transform( USAboundaries::us_states(), p4s)
us_states.p <- st_interpolate_aw( grid_popwgt.sf, us_states, extensive = T)

# merge back to state information
us_states.pop <- cbind( us_states.p, us_states[us_states.p$Group.1,])
us_states.pop$geometry.1 <- NULL
us_states.pop.dt <- data.table( us_states.pop)

# get usa mask for masking
# download USA polygon from rnaturalearth
us_states.names <- state.abb[!(state.abb %in% c( 'HI', 'AK'))]
mask.usa <- sf::as_Spatial(us_states)[ us_states$state_abbr %in% us_states.names,]

#======================================================================#
## Do the annual conversions - no differencing
#======================================================================#
# annual
if( do.annual){
  hyads_06a <- state_exposurer.year( fname = fname2006.hyads,
                                     year.m = 2006,
                                     model.use = preds.ann.hyads06w05$model.gam,
                                     name.x = 'hyads',
                                     mask.use = mask.usa,
                                     dat.a = dats2006.a,
                                     grid_pop.r = grid_popwgt.r,
                                     state_pops = copy( us_states.pop.dt))
  write.csv( file = paste0( saveloc.hyads, '_annual2006.csv'), hyads_06a)
  
  hyads_11a <- state_exposurer.year( fname = fname2011.hyads,
                                     year.m = 2011,
                                     model.use = preds.ann.hyads06w05$model.gam,
                                     name.x = 'hyads',
                                     mask.use = mask.usa,
                                     dat.a = dats2011.a,
                                     grid_pop.r = grid_popwgt.r,
                                     state_pops = copy( us_states.pop.dt))
  write.csv( file = paste0( saveloc.hyads, '_annual2011.csv'), hyads_11a)
  
  idwe_06a <- state_exposurer.year( fname = fname2011.idwe,
                                    year.m = 2011,
                                    model.use = preds.ann.idwe06w05$model.cv,
                                    name.x = 'idwe',
                                    mask.use = mask.usa,
                                    dat.a = dats2011.a,
                                    grid_pop.r = grid_popwgt.r,
                                    state_pops = copy( us_states.pop.dt))
  write.csv( file = paste0( saveloc.idwe, '_annual2006.csv'), idwe_06a)
  
  idwe_11a <- state_exposurer.year( fname = fname2011.idwe,
                                    year.m = 2011,
                                    model.use = preds.ann.idwe06w05$model.cv,
                                    name.x = 'idwe',
                                    mask.use = mask.usa,
                                    dat.a = dats2011.a,
                                    grid_pop.r = grid_popwgt.r,
                                    state_pops = copy( us_states.pop.dt))
  write.csv( file = paste0( saveloc.idwe, '_annual2011.csv'), idwe_11a)
  
  #======================================================================#
  ## Do the annual conversions - with differencing
  #======================================================================#
  # annual
  hyads_06a.d <- state_exposurer.year( fname = fname2006.hyads,
                                       year.m = 2006,
                                       model.use = preds.ann.hyads06w05$model.gam,
                                       name.x = 'hyads',
                                       mask.use = mask.usa,
                                       dat.a = dats2006.a,
                                       grid_pop.r = grid_popwgt.r,
                                       state_pops = copy( us_states.pop.dt),
                                       take.diff = T)
  write.csv( file = paste0( saveloc.hyads, '_annual_diff2006.csv'), hyads_06a.d)
  
  hyads_11a.d <- state_exposurer.year( fname = fname2011.hyads,
                                       year.m = 2011,
                                       model.use = preds.ann.hyads06w05$model.gam,
                                       name.x = 'hyads',
                                       mask.use = mask.usa,
                                       dat.a = dats2011.a,
                                       grid_pop.r = grid_popwgt.r,
                                       state_pops = copy( us_states.pop.dt),
                                       take.diff = T)
  write.csv( file = paste0( saveloc.hyads, '_annual_diff2011.csv'), hyads_11a.d)
  
  idwe_06a.d <- state_exposurer.year( fname = fname2006.idwe,
                                      year.m = 2006,
                                      model.use = preds.ann.idwe06w05$model.cv,
                                      name.x = 'idwe',
                                      mask.use = mask.usa,
                                      dat.a = dats2006.a,
                                      grid_pop.r = grid_popwgt.r,
                                      state_pops = copy( us_states.pop.dt),
                                      take.diff = T)
  write.csv( file = paste0( saveloc.idwe, '_annual_diff2006.csv'), idwe_06a.d)
  
  idwe_11a.d <- state_exposurer.year( fname = fname2011.idwe,
                                      year.m = 2011,
                                      model.use = preds.ann.idwe06w05$model.cv,
                                      name.x = 'idwe',
                                      mask.use = mask.usa,
                                      dat.a = dats2011.a,
                                      grid_pop.r = grid_popwgt.r,
                                      state_pops = copy( us_states.pop.dt),
                                      take.diff = T)
  write.csv( file = paste0( saveloc.idwe, '_annual_diff2011.csv'), idwe_11a.d)
  
  ## now idwe with gam!
  idwe_06a.dgam <- state_exposurer.year( fname = fname2006.idwe,
                                         year.m = 2006,
                                         model.use = preds.ann.idwe06w05$model.gam,
                                         name.x = 'idwe',
                                         mask.use = mask.usa,
                                         dat.a = dats2006.a,
                                         grid_pop.r = grid_popwgt.r,
                                         state_pops = copy( us_states.pop.dt),
                                         take.diff = T)
  write.csv( file = paste0( saveloc.idwe, '_annual_diffgam2006.csv'), idwe_06a.dgam)
  
  idwe_11a.dgam <- state_exposurer.year( fname = fname2011.idwe,
                                         year.m = 2011,
                                         model.use = preds.ann.idwe06w05$model.gam,
                                         name.x = 'idwe',
                                         mask.use = mask.usa,
                                         dat.a = dats2011.a,
                                         grid_pop.r = grid_popwgt.r,
                                         state_pops = copy( us_states.pop.dt),
                                         take.diff = T)
  write.csv( file = paste0( saveloc.idwe, '_annual_diffgam2011.csv'), idwe_11a.dgam)
  
}
#======================================================================#
## Do the monthly conversions
#======================================================================#
array_num <- as.numeric( Sys.getenv("SLURM_ARRAY_TASK_ID"))
array_num <- ifelse( array_num == ''|is.na( array_num), 1, array_num)
mon <- array_num %% 12

# IDWE
if( array_num == 1:12){
  idwe_exp06 <- rbindlist( lapply( mon, 
                                   state_exposurer,
                                   fstart = fstart.idwe,
                                   year.m = 2006,
                                   model.dataset = preds.mon.idwe06w05,
                                   model.name = 'model.gam', #'model.gam'
                                   name.x = 'idwe',
                                   mask.use = mask.usa,
                                   take.diff = T)) #[ mask.usa$state_abbr %in% states.use,]))
  
  write.csv( file = paste0( saveloc.idwe, '2006.csv'), idwe_exp06)
}

if( array_num == 13:24){
  idwe_exp11 <- rbindlist( lapply( mon, 
                                   state_exposurer,
                                   fstart = fstart.idwe,
                                   year.m = 2011,
                                   model.dataset = preds.mon.idwe06w05,
                                   model.name = 'model.gam', #'model.gam'
                                   name.x = 'idwe',
                                   mask.use = mask.usa,
                                   take.diff = T))
  write.csv( file = paste0( saveloc.idwe, '2011.csv'), idwe_exp11)
}

if( array_num %in% 25:36){
  # HyADS
  hyads_exp06 <- rbindlist( lapply( mon, 
                                    state_exposurer,
                                    fstart = fstart.hyads,
                                    year.m = 2006,
                                    model.dataset = preds.mon.hyads06w05,
                                    model.name = 'model.gam',
                                    name.x = 'hyads',
                                    mask.use = mask.usa,
                                    take.diff = T))
  write.csv( file = paste0( saveloc.hyads, 2006, '_', mon, '.csv'), hyads_exp06)
}

if( array_num %in% 37:48){
  hyads_exp11 <- rbindlist( lapply( mon, 
                                    state_exposurer,
                                    fstart = fstart.hyads,
                                    year.m = 2011,
                                    model.dataset = preds.mon.hyads06w05,
                                    model.name = 'model.gam',
                                    name.x = 'hyads',
                                    mask.use = mask.usa,
                                    take.diff = T))
  
  write.csv( file = paste0( saveloc.hyads, 2011, '_', mon, '.csv'), hyads_exp11)
}




#======================================================================#
## Load saved object of cmaq-ddm / hyads models from hyads_to_pm25_month.R
#======================================================================#
# load the models 
hyads.ann.model <- preds.ann.hyads06w05$model.gam
idwe.ann.model <- preds.ann.idwe06w05$model.gam

# set up the dataset
dat.coords <- coordinates( dats2006.a)
dats2006raw.dt <- data.table( cbind( dat.coords, values( dats2006.a)))
dats2006raw.dt <- data.table( cbind( dat.coords, values( dats2011.a)))

# do the predictions
hyads2006.pred <- predict( hyads.ann.model, newdata = dats2006raw.dt)
idwe2006.pred <- predict( idwe.ann.model, newdata = dats2006raw.dt)

# rasterize
dats2006.r <- rasterFromXYZ( data.table( dat.coords, hyads2006.pred, idwe2006.pred), crs = p4s)
dats2006.r[is.na( dats2006.r)] <- 0
#======================================================================#
## get total state populations
#======================================================================#
grid_popwgt.r <- rasterFromXYZ( grid_popwgt.xyz, crs = p4s)
grid_popwgt2006.r <- grid_popwgt.r$X2006

#ggplot.a.raster( unstack( grid_popwgt.r), bounds = c( 0, .5e6), facet.names = c( '2006', '2011'))



#======================================================================#
## Extra stuff -- 2006 only
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

