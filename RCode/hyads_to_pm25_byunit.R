#srun -p test --mem 100g -t 0-06:00 -c 1 -N 1 --pty /bin/bash
rm( list = ls())

platform <- c( 'mac', 'cannon')[2]
do.annual <- FALSE
do.xb <- FALSE

#coordinate reference system projection string for spatial data
p4s <- "+proj=lcc +lat_1=33 +lat_2=45 +lat_0=40 +lon_0=-97 +a=6370000 +b=6370000"

if( platform == 'mac'){
  source( '~/Dropbox/Harvard/RFMeval_Local/HyADS_to_pm25/RCode/hyads_to_pm25_functions.R')
  load( '~/Dropbox/Harvard/RFMeval_Local/HyADS_to_pm25/RData/hyads_to_cmaq_models3.RData')
  if( do.xb)
    load( '~/Dropbox/Harvard/RFMeval_Local/HyADS_to_pm25/RData/hyads_to_cmaq_models_xg.RData')
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
  load( '/n/scratchlfs/zigler_lab/lhenneman/run_hyspdisp/output_ampd_dists/hyads_to_cmaq_models2.RData')
  grid_popwgt.xyz <- fread( '/n/home03/lhenneman/inputdata/census_population/hyads_grid_population.csv',
                            drop = 'V1')
  fstart.idwe <- '/n/scratchlfs/zigler_lab/lhenneman/run_hyspdisp/output_ampd_dists/ampd_dists_sox_weighted'
  fstart.hyads <- '/n/scratchlfs/zigler_lab/lhenneman/run_hyspdisp/output_HyADS/gridexposures/GRIDexposures_byunit_'
  
  # fname2006.hyads <- '/n/scratchlfs/zigler_lab/lhenneman/run_hyspdisp/output_HyADS/2006grid/HyADSunit_grid_annual_nopbl_2006.csv'
  # fname2011.hyads <- '/n/scratchlfs/zigler_lab/lhenneman/run_hyspdisp/output_HyADS/2011grid/HyADSunit_grid_annual_nopbl_2011.csv'
  # fname.ann.idwe <- '/n/scratchlfs/zigler_lab/lhenneman/run_hyspdisp/output_HyADS/2006grid/ampd_dists_sox_weighted_annual_unit.csv'
  
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
  # hyads_06a <- state_exposurer.year( fname = fname2006.hyads,
  #                                    year.m = 2006,
  #                                    model.use = preds.ann.hyads06w05$model.gam,
  #                                    name.x = 'hyads',
  #                                    mask.use = mask.usa,
  #                                    dat.a = dats2006.a,
  #                                    grid_pop.r = grid_popwgt.r,
  #                                    state_pops = copy( us_states.pop.dt))
  # write.csv( file = paste0( saveloc.hyads, '_annual2006.csv'), hyads_06a)
  # 
  # hyads_11a <- state_exposurer.year( fname = fname2011.hyads,
  #                                    year.m = 2011,
  #                                    model.use = preds.ann.hyads06w05$model.gam,
  #                                    name.x = 'hyads',
  #                                    mask.use = mask.usa,
  #                                    dat.a = dats2011.a,
  #                                    grid_pop.r = grid_popwgt.r,
  #                                    state_pops = copy( us_states.pop.dt))
  # write.csv( file = paste0( saveloc.hyads, '_annual2011.csv'), hyads_11a)
  # 
  # idwe_06a <- state_exposurer.year( fname = fname2011.idwe,
  #                                   year.m = 2006,
  #                                   model.use = preds.ann.idwe06w05$model.cv,
  #                                   name.x = 'idwe',
  #                                   mask.use = mask.usa,
  #                                   dat.a = dats2006.a,
  #                                   grid_pop.r = grid_popwgt.r,
  #                                   state_pops = copy( us_states.pop.dt))
  # write.csv( file = paste0( saveloc.idwe, '_annual2006.csv'), idwe_06a)
  # 
  # idwe_11a <- state_exposurer.year( fname = fname2011.idwe,
  #                                   year.m = 2011,
  #                                   model.use = preds.ann.idwe06w05$model.cv,
  #                                   name.x = 'idwe',
  #                                   mask.use = mask.usa,
  #                                   dat.a = dats2011.a,
  #                                   grid_pop.r = grid_popwgt.r,
  #                                   state_pops = copy( us_states.pop.dt))
  # write.csv( file = paste0( saveloc.idwe, '_annual2011.csv'), idwe_11a)
  
  #======================================================================#
  ## Do the annual conversions - with differencing
  #======================================================================#
  # 1-pchisq(deviance( preds.ann.hyads06w05$model.cv)-deviance( preds.ann.hyads06w05$model.gam),1)
  # 1-pchisq(deviance( preds.ann.idwe06w05$model.cv)-deviance( preds.ann.idwe06w05$model.gam),1)
  # annual
  hyads_11a.d <- state_exposurer.year( fname = fname2011.hyads,
                                       year.m = 2011,
                                       model.use = preds.ann.hyads06w05$model.cv,
                                       name.x = 'hyads',
                                       mask.use = mask.usa,
                                       dat.a = dats2011.a,
                                       grid_pop.r = grid_popwgt.r,
                                       state_pops = copy( us_states.pop.dt),
                                       take.diff = T)
  write.csv( file = paste0( saveloc.hyads, '_annual_diff2011_3.csv'), hyads_11a.d$popwgt_states)
  
  hyads_06a.d <- state_exposurer.year( fname = fname2006.hyads,
                                       year.m = 2006,
                                       model.use = preds.ann.hyads06w05$model.cv,
                                       name.x = 'hyads',
                                       mask.use = mask.usa,
                                       dat.a = dats2006.a,
                                       grid_pop.r = grid_popwgt.r,
                                       state_pops = copy( us_states.pop.dt),
                                       take.diff = T)
  write.csv( file = paste0( saveloc.hyads, '_annual_diff2006_3.csv'), hyads_06a.d$popwgt_states)
  
  
  ## now idwe with gam!
  idwe_06a.dgam <- state_exposurer.year( fname = fname2006.idwe,
                                         year.m = 2006,
                                         model.use = preds.ann.idwe06w05$model.cv,
                                         name.x = 'idwe',
                                         mask.use = mask.usa,
                                         dat.a = dats2006.a,
                                         grid_pop.r = grid_popwgt.r,
                                         state_pops = copy( us_states.pop.dt),
                                         take.diff = T)
  write.csv( file = paste0( saveloc.idwe, '_annual_diff2006_3.csv'), idwe_06a.dgam$popwgt_states)
  # 
  idwe_11a.dgam <- state_exposurer.year( fname = fname2011.idwe,
                                         year.m = 2011,
                                         model.use = preds.ann.idwe06w05$model.cv,
                                         name.x = 'idwe',
                                         mask.use = mask.usa,
                                         dat.a = dats2011.a,
                                         grid_pop.r = grid_popwgt.r,
                                         state_pops = copy( us_states.pop.dt),
                                         take.diff = T)
  write.csv( file = paste0( saveloc.idwe, '_annual_diff2011_3.csv'), idwe_11a.dgam$popwgt_states)
  
  # gather the sums of all units
  hyads_06allu.r <- sum( hyads_06a.d$pred_pm.r[[1:( dim( hyads_06a.d$pred_pm.r)[3] - 1)]])
  hyads_11allu.r <- sum( hyads_11a.d$pred_pm.r[[1:( dim( hyads_11a.d$pred_pm.r)[3] - 1)]])
  idwe_06allu.r <- sum( idwe_06a.dgam$pred_pm.r[[1:( dim( idwe_06a.dgam$pred_pm.r)[3] - 1)]])
  idwe_11allu.r <- sum( idwe_11a.dgam$pred_pm.r[[1:( dim( idwe_11a.dgam$pred_pm.r)[3] - 1)]])
  
  # gather the total impacts
  hyads_06tot.r <- hyads_06a.d$zero_out.r + hyads_06allu.r
  hyads_11tot.r <- hyads_11a.d$zero_out.r + hyads_11allu.r
  idwe_06tot.r <- idwe_06a.dgam$zero_out.r + idwe_06allu.r
  idwe_11tot.r <- idwe_11a.dgam$zero_out.r + idwe_11allu.r
  
  # save the output
  hyads_06zero.dt <- data.table( rasterToPoints( hyads_06a.d$zero_out.r))[, `:=`( year = 2006, field = 'zero', model = 'hyads')]
  hyads_11zero.dt <- data.table( rasterToPoints( hyads_11a.d$zero_out.r))[, `:=`( year = 2011, field = 'zero', model = 'hyads')]
  idwe_06zero.dt  <- data.table( rasterToPoints( idwe_06a.dgam$zero_out.r))[,  `:=`( year = 2006, field = 'zero', model = 'idwe')]
  idwe_11zero.dt  <- data.table( rasterToPoints( idwe_11a.dgam$zero_out.r))[,  `:=`( year = 2011, field = 'zero', model = 'idwe')]
  setnames( hyads_06zero.dt, 'dat.pred0', 'layer')
  setnames( hyads_11zero.dt, 'dat.pred0', 'layer')
  setnames( idwe_06zero.dt,  'dat.pred0', 'layer')
  setnames( idwe_11zero.dt,  'dat.pred0', 'layer')
  
  hyads_06allu.dt <- data.table( rasterToPoints( hyads_06allu.r))[, `:=`( year = 2006, field = 'allunits', model = 'hyads')]
  hyads_11allu.dt <- data.table( rasterToPoints( hyads_11allu.r))[, `:=`( year = 2011, field = 'allunits', model = 'hyads')]
  idwe_06allu.dt  <- data.table( rasterToPoints( idwe_06allu.r))[,  `:=`( year = 2006, field = 'allunits', model = 'idwe')]
  idwe_11allu.dt  <- data.table( rasterToPoints( idwe_11allu.r))[,  `:=`( year = 2011, field = 'allunits', model = 'idwe')]
  
  hyads_06tot.dt <- data.table( rasterToPoints( hyads_06tot.r))[, `:=`( year = 2006, field = 'total', model = 'hyads')]
  hyads_11tot.dt <- data.table( rasterToPoints( hyads_11tot.r))[, `:=`( year = 2011, field = 'total', model = 'hyads')]
  idwe_06tot.dt  <- data.table( rasterToPoints( idwe_06tot.r))[,  `:=`( year = 2006, field = 'total', model = 'idwe')]
  idwe_11tot.dt  <- data.table( rasterToPoints( idwe_11tot.r))[,  `:=`( year = 2011, field = 'total', model = 'idwe')]
  
  annual_output <- rbind( hyads_06zero.dt, hyads_11zero.dt, idwe_06zero.dt, idwe_11zero.dt,
                          hyads_06allu.dt, hyads_11allu.dt, idwe_06allu.dt, idwe_11allu.dt,
                          hyads_06tot.dt, hyads_11tot.dt, idwe_06tot.dt, idwe_11tot.dt)
  write.csv( file = "~/Dropbox/Harvard/RFMeval_Local/HyADS_to_pm25/RData/annual_fields_hyads_idwe3.csv", annual_output)
  
}

#======================================================================#
## Do the annual conversions with xboost
#======================================================================#
if( do.xb){
  hyads_06a_xb <- state_exposurer.year( fname = fname2006.hyads,
                                        year.m = 2006,
                                        model.use = preds.ann.xgboost06w05$model.xg,
                                        name.x = 'hyads',
                                        mask.use = mask.usa,
                                        dat.a = dats2006.a,
                                        grid_pop.r = grid_popwgt.r,
                                        state_pops = copy( us_states.pop.dt),
                                        take.diff = T,
                                        xboost = T)
  write.csv( file = paste0( saveloc.hyads, '_annual2006_xb.csv'), hyads_06a_xb)
  
  hyads_11a_xb <- state_exposurer.year( fname = fname2011.hyads,
                                        year.m = 2011,
                                        model.use = preds.ann.xgboost06w05$model.xg,
                                        name.x = 'hyads',
                                        mask.use = mask.usa,
                                        dat.a = dats2011.a,
                                        grid_pop.r = grid_popwgt.r,
                                        state_pops = copy( us_states.pop.dt),
                                        take.diff = T,
                                        xboost = T)
  write.csv( file = paste0( saveloc.hyads, '_annual2011_xb.csv'), hyads_11a_xb)
  
  idwe_06a_xb <- state_exposurer.year( fname = fname2011.idwe,
                                       year.m = 2011,
                                       model.use = preds.ann.xgboost06w05$model.xg,
                                       name.x = 'idwe',
                                       mask.use = mask.usa,
                                       dat.a = dats2011.a,
                                       grid_pop.r = grid_popwgt.r,
                                       state_pops = copy( us_states.pop.dt),
                                       take.diff = T,
                                       xboost = T)
  write.csv( file = paste0( saveloc.idwe, '_annual2006_xb.csv'), idwe_06a_xb)
  
  idwe_11a_xb <- state_exposurer.year( fname = fname2011.idwe,
                                       year.m = 2011,
                                       model.use = preds.ann.xgboost06w05$model.xg,
                                       name.x = 'idwe',
                                       mask.use = mask.usa,
                                       dat.a = dats2011.a,
                                       grid_pop.r = grid_popwgt.r,
                                       state_pops = copy( us_states.pop.dt),
                                       take.diff = T,
                                       xboost = T)
  write.csv( file = paste0( saveloc.idwe, '_annual2011_xb.csv'), idwe_11a_xb)
  
  
}
#======================================================================#
## Do the monthly conversions
#======================================================================#
array_num <- as.numeric( Sys.getenv("SLURM_ARRAY_TASK_ID"))
array_num <- ifelse( array_num == ''|is.na( array_num), 1, array_num)
mon <- array_num %% 12
mon <- ifelse( mon == 0, 12, mon)

# IDWE
if( array_num %in% 1:12){
  idwe_exp06 <- state_exposurer( month.n = mon, 
                                 fstart = fstart.idwe,
                                 year.m = 2006,
                                 model.dataset = preds.mon.idwe06w05,
                                 model.name = 'model.cv', #'model.gam'
                                 name.x = 'idwe',
                                 mask.use = mask.usa,
                                 take.diff = T) #[ mask.usa$state_abbr %in% states.use,]))
  
  write.csv( file = paste0( saveloc.idwe, 2006, '_', mon, '_3.csv'), idwe_exp06$popwgt_states)
}

if( array_num %in% 13:24){
  idwe_exp11 <- state_exposurer( month.n = mon, 
                                 fstart = fstart.idwe,
                                 year.m = 2011,
                                 model.dataset = preds.mon.idwe06w05,
                                 model.name = 'model.cv', #'model.gam'
                                 name.x = 'idwe',
                                 mask.use = mask.usa,
                                 take.diff = T)
  write.csv( file = paste0( saveloc.idwe, 2011, '_', mon, '_3.csv'), idwe_exp11$popwgt_states)
}

if( array_num %in% 25:36){
  # HyADS
  hyads_exp06 <- state_exposurer( month.n = mon, 
                                  fstart = fstart.hyads,
                                  year.m = 2006,
                                  model.dataset = preds.mon.hyads06w05,
                                  model.name = 'model.cv',
                                  name.x = 'hyads',
                                  mask.use = mask.usa,
                                  take.diff = T)
  write.csv( file = paste0( saveloc.hyads, 2006, '_', mon, '_3.csv'), hyads_exp06$popwgt_states)
}

if( array_num %in% 37:48){
  hyads_exp11 <- state_exposurer( month.n = mon, 
                                  fstart = fstart.hyads,
                                  year.m = 2011,
                                  model.dataset = preds.mon.hyads06w05,
                                  model.name = 'model.cv',
                                  name.x = 'hyads',
                                  mask.use = mask.usa,
                                  take.diff = T)
  
  write.csv( file = paste0( saveloc.hyads, 2011, '_', mon, '_3.csv'), hyads_exp11$popwgt_states)
}




#======================================================================#
## Load saved object of cmaq-ddm / hyads models from hyads_to_pm25_month.R
#======================================================================#
# load the models 
#hyads.ann.model <- preds.ann.hyads06w05$model.gam
#idwe.ann.model <- preds.ann.idwe06w05$model.gam

# set up the dataset
#dat.coords <- coordinates( dats2006.a)
#dats2006raw.dt <- data.table( cbind( dat.coords, values( dats2006.a)))
#dats2006raw.dt <- data.table( cbind( dat.coords, values( dats2011.a)))

# do the predictions
#hyads2006.pred <- predict( hyads.ann.model, newdata = dats2006raw.dt)
#idwe2006.pred <- predict( idwe.ann.model, newdata = dats2006raw.dt)

# rasterize
#dats2006.r <- rasterFromXYZ( data.table( dat.coords, hyads2006.pred, idwe2006.pred), crs = p4s)
#dats2006.r[is.na( dats2006.r)] <- 0
#======================================================================#
## get total state populations
#======================================================================#
#grid_popwgt.r <- rasterFromXYZ( grid_popwgt.xyz, crs = p4s)
#grid_popwgt2006.r <- grid_popwgt.r$X2006

#ggplot.a.raster( unstack( grid_popwgt.r), bounds = c( 0, .5e6), facet.names = c( '2006', '2011'))



#======================================================================#
## Extra stuff -- 2006 only
#======================================================================#
# combine exposure and population
# dats2006_popwgt.r <- project_and_stack( dats2006.r, grid_popwgt.r)
# 
# # weight by 2006 population
# dats2006_popwgt.names <- names( dats2006.r)
# dats2006_popwgtexp.r <- subset( dats2006_popwgt.r, dats2006_popwgt.names)
# dats2006_popwgt.r <- dats2006_popwgtexp.r * dats2006_popwgt.r$X2006
# names( dats2006_popwgt.r) <- dats2006_popwgt.names
# 
# # take over states
# dats2006_popwgt.sf <- st_as_sf( rasterToPolygons( dats2006_popwgt.r))
# 
# #NA's coming from 
# dats2006_popwgt.states <- st_interpolate_aw( dats2006_popwgt.sf, us_states, extensive = F)

