# copied from evaluate_RFMs/plot_unit_ranks

library( data.table)
library( ggplot2)
library( viridis)
library( grid)
library( ggpubr)
library( cowplot)
library( raster)
library( sf)

## =========================================================== ##
## Define facilities plotting data
## =========================================================== ##

plot_ranked_facs <- function( ranks.dt,
                              state.abbrev = 'US',
                              rank.name = 'hyspdisp.rank',
                              size.var = 'hyspdisp.py.mean',
                              size.name = 'Pop-weighted HyADS',
                              toprank = 10,
                              geo.name = NULL,
                              xlims = NULL,
                              ylims = NULL,
                              latlonrange.year = 2005,
                              dist.scalebar = 400){
  
  # -- limit data table to units under the rank -- #
  ranks.dt.trim <- ranks.dt[get( rank.name) <= toprank & state %in% state.abbrev]
  
  # -- set name of variable size variable -- #
  setnames( ranks.dt.trim, size.var, 'size.var')
  
  # -- link with PP data if not already  -- #
  if( !( 'Longitude' %in% names( ranks.dt.trim) & 'Latitude' %in% names( ranks.dt.trim))){
    D05 <- fread("~/Dropbox/Harvard/ARP/inMAP/Merge_AMPD_NEI/final_merge_nei_ampd_all_units_2005.csv")[, year := 2005]
    D06 <- fread("~/Dropbox/Harvard/ARP/inMAP/Merge_AMPD_NEI/final_merge_nei_ampd_all_units_2006.csv")[, year := 2006]
    D11 <- fread("~/Dropbox/Harvard/ARP/inMAP/Merge_AMPD_NEI/final_merge_nei_ampd_all_units_2011.csv")[, year := 2011]
    D12 <- fread("~/Dropbox/Harvard/ARP/inMAP/Merge_AMPD_NEI/final_merge_nei_ampd_all_units_2012.csv")[, year := 2012]
    D <- rbind( D05, D06, D11, D12)
    D[, uID := gsub('_|-|\\*', '.', ID)]
    ranks.dt.units <- merge( ranks.dt.trim, D, by = c( 'uID', 'year'))
  } else
    ranks.dt.units <- ranks.dt.trim
  
  # -- find lat/lon range  -- #
  if( is.null( xlims) & is.null( ylims)){
    latlonrange <- data.table( xlim = c( min( ranks.dt.units$Longitude) - .1,
                                         max( ranks.dt.units$Longitude) + .1),
                               ylim = c( min( ranks.dt.units$Latitude - .5),
                                         max( ranks.dt.units$Latitude + .1)),
                               year = latlonrange.year)
  } else
    latlonrange <- data.table( xlim = xlims,
                               ylim = ylims,
                               year = latlonrange.year)
  
  
  # -- download states  -- #
  states <- data.table( map_data("state"))
  state <- tolower( state.name[ match( state.abbrev, state.abb)])
  state.poly <- states[region %in% state]
  
  # -- make the plot  -- #
  if( is.null( geo.name)){
    plot.title <- NULL
  } else
    plot.title <- paste("Top facilities by pop-weighted exposure on", geo.name)
  
  # -- make the plot  -- #
  gg_coal <- ggplot() + 
    theme_bw() + 
    labs(title = plot.title) +
    facet_wrap( . ~ year, ncol = 2) +
    geom_polygon(data = states, 
                 aes(x = long, y = lat, group = group),
                 fill = 'white', 
                 color = "black",
                 size = .25) + 
    geom_polygon(data = state.poly, 
                 aes(x = long, y = lat, group = group),
                 fill = 'thistle1', 
                 color = "black",
                 size = .25) + 
    coord_sf(
      xlim = latlonrange$xlim,
      ylim = latlonrange$ylim,
      datum = NA
    ) +
    geom_point(data = ranks.dt.units, 
               aes(x = Longitude, 
                   y = Latitude, 
                   size = size.var),
               stroke = 1,
               shape = 1,
               color = '#479ddd') +
    scale_size_area(guide = guide_legend(title.position = "top"),
                    name = paste( size.name, 'exposure')#,
                    # range = c(.5,5.0)
    ) +
    theme(
      plot.title = element_text(size = 16, hjust = 0.5), #element_blank(), #
      axis.title = element_text(size = 24),
      axis.text = element_blank(),
      axis.title.x = element_blank(),
      axis.title.y = element_blank(),
      legend.title = element_text(size = 10),
      legend.title.align = 0.5,
      legend.position = "bottom", #c(.22, .15),
      legend.text = element_text(size = 8, angle = 0),
      legend.background = element_rect(fill = 'transparent'),
      legend.key.size = unit(.05, 'npc'),
      legend.direction = 'horizontal',
      # rect = element_blank(), #( fill = 'transparent'),
      strip.text = element_text( size = 14),
      strip.background = element_rect( fill = 'white')
    )  + 
    geom_rect( data = latlonrange,
               aes(xmin = xlim[1] - 5, 
                   xmax = xlim[1] + (xlim[2] - xlim[1]) / 2, 
                   ymin = ylim[1] - 5, 
                   ymax = ylim[1] + .5), 
               fill = 'white', 
               color = NA) +
    ggsn::scalebar( location = 'bottomleft',
                    anchor = c( x = latlonrange$xlim[1] + .2, y = latlonrange$ylim[1] + .4), 
                    x.min = latlonrange$xlim[1],
                    y.min = latlonrange$ylim[1],
                    x.max = latlonrange$xlim[2],
                    y.max = latlonrange$ylim[2],
                    dist = dist.scalebar / 2, 
                    height = 0.02, 
                    st.dist = 0.04, 
                    st.size = 3, 
                    dd2km = TRUE, 
                    model = 'WGS84',
                    facet.var = 'year',
                    facet.lev = as.character( latlonrange.year))
  
  print( gg_coal)
  return( list( dt = ranks.dt.units,
                plot = gg_coal,
                latlonrange = copy( latlonrange)[,year := NULL]))
  
}

#======================================================================#
# Create multiple plots
#======================================================================#
plot_multipleconnect4 <- function(state.use = 'GA',
                                  dt = ranks_adj_all, 
                                  plot.metrics =  c( 'Adjoint' = 'SOx exposure',
                                                     'HyADS' = 'hyspdisp.py.mean', 
                                                     'IDW' = 'sox.inverse.dist.py.mean'),
                                  years.use = c( 2005, 2012)
){
  
  if( state.use == 'allstates') {
    legend.pos <- 'bottom'
    dt_use <- na.omit( copy( dt[ year %in% years.use]))
  } else {
    legend.pos <- 'none'
    dt_use <- na.omit( copy( dt[ state %in% state.use & year %in% years.use]))
  }
  
  
  plot.metrics.names <- names( plot.metrics)
  setnames(dt_use, plot.metrics, names( plot.metrics))
  
  take.corrs <- function( plot.metric.name,
                          corr.type = 'pearson'){
    corrs <- dt_use[, lapply( plot.metrics.names,
                              function(name.met) 
                                cor( get( name.met), 
                                     get( plot.metric.name),
                                     use = 'complete.obs',
                                     method = corr.type)),
                    by = .( state.factor, year)]
    setnames(corrs, paste0( 'V', 1:length( plot.metrics.names)), plot.metrics.names)
    corrs[, `:=` (base.model = plot.metric.name,
                  type = corr.type)]
  }
  
  corrs.dt.p <- rbindlist( lapply( plot.metrics.names, take.corrs))
  corrs.dt.s <- rbindlist( lapply( plot.metrics.names, take.corrs, 'spearman'))
  corrs.dt.m <- melt( rbind( corrs.dt.p, corrs.dt.s),
                      id.vars = c( 'state.factor', 'year', 'base.model', 'type'))
  
  corrs.dt.mu <- unique( corrs.dt.m, by = c( 'base.model', 'type', 'state.factor', 'year'))[base.model != variable]
  
  
  gg <- ggplot( corrs.dt.mu, 
                aes( x = base.model,
                     y = variable,
                     size = type,
                     color = value
                )) +
    facet_grid( rows = vars( state.factor),
                cols = vars( year)) +
    geom_point() + #position = position_dodge(width = 1)
    scale_size_discrete( range = c(15,8),
                         guide = guide_legend(override.aes = list(
                           size = c(10,5)))) +
    scale_y_discrete("Model 1") +
    scale_x_discrete("Model 2") +
    scale_color_viridis( limits = c(0.8,1), 
                         oob = scales::squish,
                         option = 'C') +
    theme_bw() +
    theme(plot.title = element_text(size = 16, hjust = 0.5),
          axis.text = element_text( size = 8),
          axis.text.x = element_text( angle = 0),
          axis.text.y = element_blank(),
          axis.ticks = element_blank(),
          axis.title = element_blank(),
          legend.position = legend.pos,	
          legend.title = element_blank(),
          legend.text = element_text(size = 8),
          legend.background = element_rect(fill = 'transparent'),
          legend.key = element_rect(fill = 'transparent'),
          legend.key.width = unit(.5, 'cm'),
          legend.key.height = unit(.04, 'npc'),
          legend.box = 'vertical',
          legend.box.just = 'top',
          legend.box.background = element_rect(),
          legend.box.margin = margin(0, 0, 0, 0),
          legend.direction = 'horizontal',
          legend.spacing = unit(-.1,'npc'),
          legend.margin = margin(t = -0.02, r = .25, b = 0.05, l = .01, unit = "cm"),
          panel.grid = element_blank(),
          strip.background = element_rect(fill='white'),
          strip.text = element_text(size=10,
                                    margin = margin( 0, 0, 0, 0))
    )
  
  print(gg)
  return( gg)
}


## =========================================================== ##
## Load adjoint runs, perform rankings, merge
## =========================================================== ##

read.adj <- function( base.dir = '~/Dropbox/Harvard/RFMeval_Local/Adjoint',
                      adj.name = 'other_states'){
  
  # ID directory, list files
  adj.dir <- file.path( base.dir,
                        adj.name)
  adj.files <- data.table( file = list.files( adj.dir,
                                              full.names = T))
  
  # ID states and years
  adj.files[, `:=` ( state = gsub( '.*_adjoint_|.csv', 
                                   '', 
                                   file),
                     year = gsub( '.*_all_units_|_adjoint.*.csv', 
                                  '', 
                                  file))]
  adj.files[, state := gsub( 'USA', 'US', state)]
  
  # read in file that includes names, assign these names to read files
  adj.withnames <- fread( '~/Dropbox/Harvard/RFMeval_Local/Adjoint/final_merge_nei_ampd_all_units_2005_adjoint.csv')
  adj.names <- names( adj.withnames)
  
  # Read in, listify
  adj.in <- lapply( 1:nrow( adj.files),
                    function( n, dt, adj.names.in) { 
                      x <- dt[n]
                      adj.in.n <- fread( x$file)
                      names( adj.in.n) <- adj.names.in
                      
                      adj.in.n[, `:=` (state = x$state,
                                       year = as( x$year,
                                                  'numeric'))]
                      return( adj.in.n)
                    },
                    adj.files,
                    adj.names)
  adj.in.dt <- rbindlist( adj.in)
  
  # add in adjoint name
  adj.in.dt[, `:=` (adj.name = adj.name)]
  
  # minor adjustments, rank by SOx exposure
  adj.in.dt[, uID := gsub( ' |"', '', ID)]
  adj.in.dt[, uID := gsub('_|-|\\*', '.', uID)]
  
  # minor adjustments, rank by SOx exposure
  adj.in.dt <- adj.in.dt[ SOx > 0]
  adj.in.dt[, SOx.rank.adj  := frankv( `SOx exposure`, order = -1),
            by = .( year, state)]
  
  return( adj.in.dt[, .( uID, adj.name, year, 
                         state, SOx,
                         `SOx exposure`,
                         SOx.rank.adj)])
}

## =========================================================== ##
## Take the metrics
## =========================================================== ##
eval.fn <- function( Yhat, Yact, mod.name){
  num.diff <- sum( Yhat - Yact, na.rm = T)
  abs.diff <- sum( abs( Yhat - Yact), na.rm = T)
  denom <- sum( Yact, na.rm = T)
  metrics <- data.table( mod.name = mod.name,
                         NMB = num.diff / denom,
                         NME = abs.diff / denom,
                         MB   = num.diff / length( Yhat),
                         RMSE = sqrt( sum( ( Yhat - Yact) ^ 2, na.rm = T) / length( Yhat)),
                         R.p = cor( Yhat, Yact, use = 'complete.obs') ^ 2,
                         R.s = cor( Yhat, Yact, use = 'complete.obs', method = 'spearman') ^ 2)
  return( metrics)
}


## =========================================================== ##
## Read in population data, allocate to states
## =========================================================== ##
grid_popwgt.xyz <- fread( '~/Dropbox/Harvard/RFMeval_Local/HyADS_to_pm25/HyADS_grid/population/hyads_grid_population.csv',
                          drop = 'V1')

# rasterize
p4s <- "+proj=lcc +lat_1=33 +lat_2=45 +lat_0=40 +lon_0=-97 +a=6370000 +b=6370000"
grid_popwgt.r <- rasterFromXYZ( grid_popwgt.xyz, crs = p4s)

# as sf
grid_popwgt.sf <- st_as_sf( rasterToPolygons( grid_popwgt.r))

# get total state populations
us_states <- st_transform( USAboundaries::us_states(), p4s)
us_states.p <- st_interpolate_aw( grid_popwgt.sf, us_states, extensive = T)

# merge back to state information
us_states.pop <- cbind( us_states.p, us_states[us_states.p$Group.1,])
us_states.pop$geometry.1 <- NULL
us_states.pop.dt <- data.table( us_states.pop)

# melt to long format
us_states.pop.dt.m <- melt( us_states.pop.dt[,.( X2006, X2011, state_abbr)],
                            id.vars = 'state_abbr', value.name = 'pop', variable.name = 'year')
us_states.pop.dt.m[, year := as( gsub( 'X', '', year), 'numeric')]

# sum to US
us_states.pop.dt.us <- us_states.pop.dt.m[, .( pop = sum( pop)), by = year][, state_abbr := 'US']
us_states.pop.dt.m <- rbind( us_states.pop.dt.m, us_states.pop.dt.us)


## =========================================================== ##
## Read in adjoint
## =========================================================== ##
adj.other_states <-       read.adj( base.dir = '~/Dropbox/Harvard/RFMeval_Local/Adjoint/redraftindividualsourceevaluation',
                                    adj.name = 'initial')
adj.layers_2_5 <-         read.adj( base.dir = '~/Dropbox/Harvard/RFMeval_Local/Adjoint/redraftindividualsourceevaluation',
                                    adj.name = 'layers_2-5')
adj.stack_height <-       read.adj( base.dir = '~/Dropbox/Harvard/RFMeval_Local/Adjoint/redraftindividualsourceevaluation',
                                    adj.name = 'stack_height')
adj.stack_height_plus1 <- read.adj( base.dir = '~/Dropbox/Harvard/RFMeval_Local/Adjoint/redraftindividualsourceevaluation',
                                    adj.name = 'stack_height_plus1')
adj.stack_height_plus2 <- read.adj( base.dir = '~/Dropbox/Harvard/RFMeval_Local/Adjoint/redraftindividualsourceevaluation',
                                    adj.name = 'stack_height_plus2')

adjoint_all <- rbind( adj.other_states,
                      adj.layers_2_5,
                      adj.stack_height,
                      adj.stack_height_plus1,
                      adj.stack_height_plus2)

#======================================================================#
## area weight over states
#======================================================================#
# combine exposure and population
adjoint_all_pop <- merge( us_states.pop.dt.m, adjoint_all, 
                          by.x = c( 'state_abbr', 'year'), by.y = c( 'state', 'year'))
setnames( adjoint_all_pop, 'state_abbr', 'state')
adjoint_all_pop[, adj_popwgt := `SOx exposure` / pop] 

## =========================================================== ##
## adjoint variability between plume assumptions
## =========================================================== ##
adj.variability.name <- adjoint_all_pop[ , sum( adj_popwgt), 
                                         by = .( adj.name, year, state)]#[state == 'US' & year == 2006]
adj.variability.name[, fract.V1 := V1 / min(V1)]

adj.variability.year <- adjoint_all_pop[ , sum( `adj_popwgt`), 
                                         by = .( adj.name, year, state)][state == 'US' & 
                                                                           year %in% c(2006,2011) & 
                                                                           adj.name == 'initial']
adj.variability.year[, fract.V1 := V1 / min(V1[1])]

## =========================================================== ##
## Load multiple states - HyADS and IDWE
## =========================================================== ##
exp_home <- '~/Dropbox/Harvard/RFMeval_Local/HyADS_to_pm25/PopWeighted'
exp_home_annual <- '~/Dropbox/Harvard/RFMeval_Local/HyADS_to_pm25/RData'
dropcols <- c( 'V1', 'mean_popwgt', 'pop_amnt')

# IDWE
idwe_2006 <- fread( file.path( exp_home_annual, 'popwgt_idwe_annual_diff2006_2.csv'), drop = dropcols)
idwe_2011 <- fread( file.path( exp_home_annual, 'popwgt_idwe_annual_diff2011_2.csv'), drop = dropcols)
idwe_all <- rbind( idwe_2006, idwe_2011)[uID != 'Xtot.sum']
setnames( idwe_all, c( 'popwgt', 'mean_pm'), c( 'idwe.pw', 'idwe.mean'))

# HyADS
hyads_2006 <- fread( file.path( exp_home_annual, 'popwgt_hyads_annual_diff2006_2.csv'), drop = dropcols)
hyads_2011 <- fread( file.path( exp_home_annual, 'popwgt_hyads_annual_diff2011_2.csv'), drop = dropcols)
hyads_all <- rbind( hyads_2006, hyads_2011)
setnames( hyads_all, c( 'popwgt', 'mean_pm'), c( 'hyads.pw', 'hyads.mean'))

# hyads xgboost
hyadsx_2006 <- fread( file.path( exp_home_annual, 'popwgt_hyads_annual2006_xb.csv'), drop = dropcols)
hyadsx_2011 <- fread( file.path( exp_home_annual, 'popwgt_hyads_annual2011_xb.csv'), drop = dropcols)
hyadsx_all <- rbind( hyadsx_2006, hyadsx_2011)
setnames( hyadsx_all, 'popwgt', 'hyadsx')

# hyads xgboost
idwex_2006 <- fread( file.path( exp_home_annual, 'popwgt_idwe_annual2006_xb.csv'), drop = dropcols)
idwex_2011 <- fread( file.path( exp_home_annual, 'popwgt_idwe_annual2011_xb.csv'), drop = dropcols)
idwex_all <- rbind( idwex_2006, idwex_2011)
setnames( idwex_all, 'popwgt', 'idwex')

# combine, limit to input states
states.use <- c( 'PA', 'KY', 'GA', 'WI', 'TX', 'CO', 'CA', 'US')
rcm_exp <- Reduce(function(...) merge(..., all = TRUE, by = c( 'state_abbr', 'uID', 'year', 'ID')), 
                  list( idwe_all, hyads_all))[ state_abbr %in% states.use]
rcm_exp[, uID := gsub( '^X', '', uID)]
setnames( rcm_exp, 'state_abbr', 'state')


ggplot( data = rcm_exp[ state %in% states.use]) +
  geom_point( aes( x = hyads, y = idwe, color = year)) + 
  # scale_y_continuous( limits = c( 0,1)) + 
  facet_wrap( . ~ state)
## =========================================================== ##
## Merge adjoint and HyADS/dist
## =========================================================== ##
ranks_adj_all <- merge( adjoint_all_pop, rcm_exp,
                        by = c( "uID", "year", "state"),
                        all = T, allow.cartesian = TRUE)
ranks_adj_all[, state.factor := factor(`state`, levels =  c( 'US', 'PA', 'KY', 'GA', 'NY', 'WI', 'TX', 'CO', 'CA'))]

# not looking at the ones with zero
ranks_adj_all <- ranks_adj_all[ uID != 'tot.sum']

# rank hyads and idwe
# ranks_adj_all[, `:=` ( hyads.rank  = frankv( `hyads`, order = -1),
#                        idwe.rank   = frankv( `idwe`, order = -1),
#                        hyadsx.rank  = frankv( `hyadsx`, order = -1),
#                        idwex.rank   = frankv( `idwex`, order = -1)),
#               by = .( year, state, adj.name)]

## =========================================================== ##
## Total contributions to each state - how do they compare to CMAQ-DDM?
## =========================================================== ##
# sum by states
rcm_statesums <- na.omit( ranks_adj_all[ adj.name == 'initial'])[, 
                                         .( adj.pw = sum( adj_popwgt),
                                              idwe.pw    = sum( idwe.pw),
                                              idwe.mean  = sum( idwe.mean),
                                              hyads.pw   = sum( hyads.pw),
                                              hyads.mean = sum( hyads.mean)),
                                         by = c( 'state', 'year')]

# load CMAQ data
load( '~/Dropbox/Harvard/RFMeval_Local/HyADS_to_pm25/RData/hyads_to_cmaq_models2.RData')

# grid-weighted population as raster
ddm2006.sf <- st_as_sf( rasterToPolygons( dats2006.a$cmaq.ddm))
ddm2011.sf <- st_as_sf( rasterToPolygons( dats2011.a$cmaq.ddm))

# get total state populations
us_states <- st_transform( USAboundaries::us_states(), p4s)
ddm2006.states <- data.table( st_interpolate_aw( ddm2006.sf, us_states, extensive = F))
ddm2011.states <- data.table( st_interpolate_aw( ddm2011.sf, us_states, extensive = F))

# merge together
ddm2006.states[, `:=` ( year = 2006, state = us_states$state_abbr[Group.1], Group.1 = NULL, geometry = NULL)]
ddm2011.states[, `:=` ( year = 2011, state = us_states$state_abbr[Group.1], Group.1 = NULL, geometry = NULL)]
ddm.states <- rbind( ddm2006.states, ddm2011.states)
rcm_statesums <- merge( rcm_statesums, ddm.states, by = c( 'state', 'year'))

# load in model fits
modfits.annual <- fread( "~/Dropbox/Harvard/RFMeval_Local/HyADS_to_pm25/RData/annual_fields_hyads_idwe.csv", drop = 'V1')
modfits.annual.c <- dcast( modfits.annual, x + y + year ~ field, value.var = 'layer')


## =========================================================== ##
## Mkae plots
## =========================================================== ##

gghyads.lin <- ggplot( data = ranks_adj_all[adj.name == 'initial' & year %in% c( 2006, 2012)],
                       aes( x = adj_popwgt,
                            y = hyads)) +
  geom_point() +
  geom_smooth( method = 'lm') +
  facet_grid( rows = vars( state.factor),
              cols = vars( year)) +
  ylab( "HyADS unit PWE") +
  xlab( "Adjoint unit PWE") +
  theme_bw() + 
  theme( axis.text.x = element_text( angle = 30,
                                     vjust = .75))

gghyads.rank <- ggplot( data = ranks_adj_all[adj.name == 'initial' & year %in% c( 2006, 2012)],
                        aes( x = SOx.rank.adj,
                             y = hyads.rank)) +
  geom_point() +
  # scale_x_continuous( limits = c( 0, 1060)) +
  # scale_y_continuous( limits = c( 0, 1060)) +
  geom_smooth( method = 'lm') +
  facet_grid( rows = vars( state.factor),
              cols = vars( year)) +
  ylab( "HyADS unit rank") +
  xlab( "Adjoint unit rank") +
  theme_bw() 

ggIDWE.lin <- ggplot( data = ranks_adj_all[adj.name == 'initial' & year %in% c( 2006, 2012)],
                      aes( x = adj_popwgt,
                           y = idwe)) +
  geom_point() +
  geom_smooth( method = 'lm') +
  facet_grid( rows = vars( state.factor),
              cols = vars( year)) +
  ylab( "IDWE unit ranks") +
  xlab( "Adjoint unit ranks") +
  theme_bw() + 
  theme( axis.text.x = element_text( angle = 30,
                                     vjust = .75))

ggIDWE.rank <- ggplot( data = ranks_adj_all[adj.name == 'initial' & year %in% c( 2006, 2012)],
                       aes( x = SOx.rank.adj,
                            y = idwe.rank)) +
  geom_point() +
  # scale_x_continuous( limits = c( 0, 1060)) +
  # scale_y_continuous( limits = c( 0, 1060)) +
  geom_smooth( method = 'lm') +
  facet_grid( rows = vars( state.factor),
              cols = vars( year)) +
  ylab( "HyADS unit rank") +
  xlab( "Adjoint unit rank") +
  theme_bw() 

## =========================================================== ##
## Plot fraction of emissions vs. exposure fraction
## =========================================================== ##
eval.hyads <- ranks_adj_all[, eval.fn( hyads, adj_popwgt, 'HyADS'), by = .( state, year, adj.name)]
eval.idwe  <- ranks_adj_all[, eval.fn( idwe, adj_popwgt, 'IDWE'), by = .( state, year, adj.name)]
eval.all <- rbind( eval.hyads, eval.idwe)

# melt
eval.all.m <- melt( eval.all, id.vars = c( 'state', 'year', 'adj.name', 'mod.name'),
                    variable.name = 'metric')

## define adjoint names, set up for plotting
adj.names <- data.table( adj.name = c( "initial", "layers_2-5", "stack_height", "stack_height_plus1", "stack_height_plus2"),
                         adj.name.adjust = c("Average", "Layers 2-5", "Stack Height", "Stack Height +1", "Stack Height +2"))
corrs.all <- merge( eval.all.m, adj.names, by = 'adj.name')
corrs.all[, state.factor := factor(`state`, levels = c( 'US', 'CA', 'CO', 'TX', 'WI', 'NY', 'GA', 'KY', 'PA'))]

ggcorrs <- ggplot( data = corrs.all,
                   aes( x = state.factor,
                        y = value,
                        color = mod.name,
                        shape = adj.name.adjust
                   )) +
  geom_point( size = 4,
              position = position_dodge( width = .25)) + 
  scale_shape_manual( values = c( 0:2, 5:6)) +
  scale_color_manual( values = c( '#479ddd', '#dd8747')) +
  geom_vline( xintercept = 1.5,
              lty = 2) +
  facet_grid( rows = vars( metric),
              cols = vars( year), scales = 'free_y') +
  expand_limits( y = 0) +
  # scale_y_continuous( limits = c( 0.00, 1)) + 
  guides(
    color = guide_legend(order = 1,
                         keywidth = 1.5),
    shape = guide_legend(order = 0,
                         keywidth = 1.5)
  ) +
  theme_bw() +
  theme( axis.text = element_text( size = 12),
         axis.title = element_text( size = 18),
         axis.title.x = element_blank(),
         legend.direction = 'horizontal',
         legend.position = 'bottom',
         legend.text = element_text( size = 12),
         legend.title = element_blank(),
         panel.grid.major.x = element_line( size = 10,
                                            color = 'grey90'),
         panel.grid.minor = element_blank(),
         strip.background = element_blank(),
         strip.text = element_text( size = 18))

ggsave( file = '~/Dropbox/Harvard/Manuscripts/eval_local_exposure/figures/corrs_adjoints_heights_2006_2011.png',
        ggcorrs,
        width  = 8 * 1.1,
        height = 5 )


corrs.all[adj.name == 'initial' & state == 'US']
corrs.all[adj.name == 'initial' & state %in% c( 'KY', 'PA', 'GA')]

## =========================================================== ##
## Check out average distance metric
## =========================================================== ##
ave_distances <- fread( '~/Dropbox/Harvard/RFMeval_Local/HyADS/unit_rankings/ave_distances.csv')
ave_distances[, state.factor := factor(`state`, levels = c( 'US', 'CA', 'CO', 'TX', 'WI', 'NY', 'GA', 'KY', 'PA'))]
ave_distances[, `:=` (avedist.pew.relative = avedist.pew / ave_distances[state == 'US']$avedist.pew,
                      avehyads.pew.relative = avehyads.pew / ave_distances[state == 'US']$avehyads.pew)]


gg_avedist <- ggplot( data = ave_distances,
                      aes( x = state.factor,
                           y = avedist.pew)) + 
  geom_col( width = .6,
            fill = 'grey50') + 
  geom_text( aes( label = round( avedist.pew)),
             y = 250,
             color = 'white',
             fontface = 'bold') +
  labs( y = expression( paste( D^{pew}, ' [km]'))) +
  geom_vline( xintercept = 1.5,
              lty = 2) +
  theme_bw( ) + 
  theme( axis.title = element_text( size = 18),
         axis.title.x = element_blank(),
         axis.text = element_text( size = 12),
         panel.grid.major.x = element_blank())


ave_distances.m <- melt( ave_distances,
                         id.vars = 'state.factor',
                         measure.vars = c( 'avedist.pew.relative',
                                           'avehyads.pew.relative'))
gg_avedist.rel <- ggplot( data = ave_distances.m,
                          aes( x = state.factor,
                               y = value,
                               fill = variable)) + 
  geom_col( width = .6,
            position = position_dodge2()) + 
  labs( y = expression( paste( D^{pew}, ' & ', H^{pew}, ' relative to US'))) +
  geom_vline( xintercept = 1.5,
              lty = 2) +
  theme_bw( ) + 
  theme( axis.title = element_text( size = 18),
         axis.title.x = element_blank(),
         axis.text = element_text( size = 12),
         panel.grid.major.x = element_blank())



gg_combine <- plot_grid(ggcorrs, gg_avedist,
                        rel_heights = c( 3, 1),
                        labels = NULL, ncol = 1, 
                        align = 'hv', axis = 'lr') 

ggsave( file = '~/Dropbox/Harvard/Manuscripts/eval_local_exposure/figures/corrs_adjoints_heights_2006_2011_dists.png',
        gg_combine,
        width  = 8 * 1.1,
        height = 5 * 1.3)

## =========================================================== ##
## Plot fraction of emissions vs. exposure fraction
## =========================================================== ##
ranks_adj_fracs.corrs <- ranks_adj_all[, .( hya = cor( hya.sox.fract, adj.sox.fract),
                                            iwd = cor( iwd.sox.fract, adj.sox.fract)), 
                                       by = .( year, state.factor)]


ranks_adj_all.m <- melt( ranks_adj_all[, .( uID, year, state.factor, 
                                            adj.sox.fract, hya.sox.fract, 
                                            iwd.sox.fract, emiss.sox.fract)],
                         id.vars = c( 'uID', 'year', 
                                      'state.factor', 'emiss.sox.fract'))

'%ni%' <- function(x,y)!('%in%'(x,y))
ggplot( data = ranks_adj_all.m[state.factor %ni% c('NY', 'CA') & 
                                 year == 2005],
        aes( x = emiss.sox.fract,
             y = value,
             color = variable)) +
  geom_point( size = 0.1) +
  geom_smooth( method = 'lm') +
  facet_grid( rows = vars( state.factor),
              cols = vars( year)) +
  scale_x_continuous( trans = 'log',
                      breaks = c(1e-1, 1e-3, 1e-5, 1e-7, 1e-9, 1e-11),
                      labels = scales::scientific) +
  scale_y_continuous( trans = 'log',
                      breaks = c(1e-1, 1e-3, 1e-5, 1e-7, 1e-9, 1e-11),
                      labels = scales::scientific) +
  ylab( "Fraction of exposure") +
  xlab( "Fraction of total emissions") +
  theme_bw() + 
  theme( axis.title = element_text( size = 10),
         axis.text = element_text( size = 8),
         strip.background = element_rect( fill = 'white'))

## =========================================================== ##
## Check out monthly rankings
## =========================================================== ##
library( hyspdisp)
units.use <- units2011$uID

# read, create monthly factor
ranks_allstates.mo2011 <- fread( "~/Dropbox/Harvard/RFMeval_Local/HyADS/unit_rankings//ranks_all_mo_2006_2011.csv")
ranks_allstates.mo <- fread( "~/Dropbox/Harvard/RFMeval_Local/HyADS/unit_rankings//ranks_all_mo_2006_2011.csv")
#ranks_allstates.mo <- rbind( ranks_allstates.mo, ranks_allstates.mo2011)[uID %in% units.use]
ranks_allstates.mo[, `:=` ( month =  factor( month.name[ as.integer( gsub( '.*_', '', year_month))],
                                             levels = month.name),
                            state.factor = factor( state,
                                                   levels = c( 'US', 'PA', 'KY', 'GA', 'NY', 'WI', 'TX', 'CO', 'CA')))]



# calculate correlations, slopes, etc by year, state, month, etc
ranks_all.mo.cors <- ranks_allstates.mo[, .( slope = coef( lm( hyspdisp.rank ~ sox.inverse.dist.rank))[2],
                                             inter = coef( lm( hyspdisp.rank ~ sox.inverse.dist.rank))[1],
                                             Pearson= cor( hyspdisp.py.mean, sox.inverse.dist.py.mean),
                                             Spearman= cor( hyspdisp.py.mean, sox.inverse.dist.py.mean, method = 'spearman')), 
                                        by = .( year, month, state.factor)]
ranks_all.mo.cors.m <- melt( ranks_all.mo.cors,
                             id.vars = c( 'year', 'month', 'state.factor'))

ranks_allstates.mo[state == 'GA' & year_month == '2012_07' & hyspdisp.rank > 800]
ranks_all.mo.cors[state.factor == 'US', min(Pearson), by = .( year)]
ranks_all.mo.cors[state.factor == 'US', max(Pearson), by = .( year)]
ranks_all.mo.cors[state.factor == 'US', min(Spearman), by = .( year)]
ranks_all.mo.cors[state.factor == 'US', max(Spearman), by = .( year)]
ranks_all.mo.cors[state.factor == 'CA', max(Spearman), by = .( year)]

ggplot( data = ranks_allstates.mo[state == 'GA'],
        aes( x = hyspdisp.rank,
             y = sox.inverse.dist.rank)) +
  geom_point() +
  facet_wrap( . ~ year_month, ncol = 2)

# create rectangles for shading

rect_bot <- c( 0.25, 0.75)
rectangles <- data.frame(
  xmin = rect_bot,
  xmax = rect_bot + .25,
  ymin = -Inf,
  ymax = +Inf,
  fill = 'grey50'
)

# plot correlations across the year
ggcors <- ggplot( data = ranks_all.mo.cors.m[variable %in% c( 'Spearman', 'Pearson') &
                                               year %in% c( 2006, 2011) &
                                               state.factor != 'NY'],
                  aes( x = month,
                       y =  value, #substring( value, 2)
                       color = variable,
                       group = variable)) +
  geom_rect( data = rectangles,
             aes( ymin = xmin,
                  ymax = xmax),
             fill = 'grey95',
             xmin = -Inf,
             xmax = +Inf,
             inherit.aes = FALSE) +
  geom_hline( yintercept = 1,
              size = 0.5,
              color = 'grey50') + 
  geom_hline( yintercept = 0,
              size = 0.5,
              color = 'grey50') + 
  geom_line( size = 1.5) +
  scale_y_continuous( limits = c( NA,1)) +
  facet_grid( rows = vars( state.factor),
              cols = vars( year)) +
  scale_color_viridis( begin = 0.4, 
                       end = 0.8,
                       discrete = T) +
  ylab( "Correlations") +
  xlab( "Month") +
  theme_bw() +
  theme(plot.title = element_text(size = 16, hjust = 0.5),
        axis.text = element_text( size = 14),
        axis.text.x = element_text( angle = 45,
                                    vjust = 0.6),
        axis.title = element_text( size = 20),
        axis.title.x = element_blank(),
        legend.position = c( .19, .92),	
        legend.title = element_blank(),
        legend.text = element_text(size = 16),
        legend.box.background = element_rect(),
        legend.box.margin = margin(0, 0, 0, 0),
        legend.direction = 'horizontal',
        panel.grid = element_blank(),
        strip.background = element_blank(), #rect(fill='white'),
        strip.text = element_text(size=24))


ggsave( file = '~/Dropbox/Harvard/Manuscripts/eval_local_exposure/figures/correlations_monthly_2006_2011.png',
        plot = ggcors,
        width =  7 * 1.3,
        height = 7 * 1.7)


ggplot( data = ranks_allstates.mo,
        aes( x = hyspdisp.py.mean,
             y = sox.inverse.dist.py.mean,
             color = state.factor
        )) +
  scale_color_viridis( discrete = T) +
  geom_point( alpha = .3,
              size = .05) +
  #  geom_smooth( method = 'lm') +
  facet_grid( rows = vars( month),
              cols = vars( year)) +
  theme_bw()

ggplot( data = ranks_allstates.mo,
        aes( x = hyspdisp.rank,
             y = sox.inverse.dist.rank,
             color = state.factor)) +
  scale_color_viridis( discrete = T) +
  geom_point( alpha = .3,
              size = .05) +
  geom_abline( intercept = 0, 
               slope = 1,
               color = 'red') +
  geom_smooth( method = 'lm') +
  facet_grid( rows = vars( month),
              cols = vars( year)) +
  # ylab( "HyADS unit ranks") +
  # xlab( "SOx inverse dist unit ranks") +
  theme_bw()



## =========================================================== ##
## Average distance by state
## =========================================================== ##
zips_inv_dist <- fread( "~/Dropbox/Harvard/RFMeval_Local/Comparisons_Intermodel/evaluate_RFMs_intermediates/ampd_zips_links_annual_5km.csv")

## =========================================================== ##
## Check decreasing exposure by distance
## =========================================================== ##
ziplinks_sub <- fread( "~/Dropbox/Harvard/RFMeval_Local/HyADS/unit_rankings/ziplinks_hyads_ampd_sub.csv")

ziplinks_sub.m <- melt( ziplinks_sub[, .( uID, year, hyspdisp, dist.km,
                                          sox.inverse.dist)],
                        id.vars = c( 'uID', 'year', 'dist.km'))

ziplinks_sub.m[, val.scale := value / value[which.min( dist.km)],
               by = .( uID, year, variable)]

ggplot( data = ziplinks_sub.m[year == 2005],
        aes( x = dist.km,
             y = val.scale)) +
  scale_color_viridis( discrete = F) +
  geom_point( ) +
  scale_x_continuous( limits = c( 0, 500)) +
  facet_grid( rows = vars( uID),
              cols = vars( variable)) +
  # ylab( "HyADS unit ranks") +
  # xlab( "SOx inverse dist unit ranks") +
  theme_bw()

## =========================================================== ##
## What's going on in WI?
## =========================================================== ##
ggplot( data = adjoint_all[ state %in% c('WI', 'KY', 'GA')],
        aes( x = SOx,
             y = `SOx exposure`,
             color = adj.name)) + 
  geom_point() + 
  facet_grid( rows = vars( state),
              cols = vars( year))


## =========================================================== ##
## Presentation plots
## =========================================================== ##

ranks.melt <- melt( ranks_KY_adj,
                    id.vars = c( 'uID', 'year', 'SOx.rank.adj'))

r2.1 <- format( cor( ranks_KY_adj[year == 2005, 
                                  SOx.rank.adj], 
                     ranks_KY_adj[year == 2005, 
                                  hyspdisp.rank])^2, digits = 2)
r2.2 <- format( cor( ranks_KY_adj[year == 2005, 
                                  SOx.rank.adj], 
                     ranks_KY_adj[year == 2005, 
                                  sox.inverse.dist.rank])^2, digits = 2)

hospital_names <- c(
  `HyADS~Ranks` = "hyspdisp.rank",
  `SO[x] * distance^{'-1'}` = "sox.inverse.dist.rank",
  `Hospital#3` = "Hospital Number 3",
  `Hospital#4` = "The Other Hospital"
)

ggplot( data = ranks.melt[ year == 2005 & variable %in% c( 'hyspdisp.rank',
                                                           'sox.inverse.dist.rank')],
        aes( x = SOx.rank.adj,
             y = value)) +
  geom_point() +
  geom_abline( intercept = 0, 
               slope = 1,
               color = 'red') +
  facet_wrap( . ~ variable, ncol = 2) +#, labeller = as_labeller(hospital_names)) +
  #  ylab( "HyADS unit ranks") +
  xlab( "Adjoint unit ranks") +
  theme_bw() +
  theme( axis.title.y = element_blank(),
         #   strip.text = element_blank(),
         strip.background = element_blank())

gg.annual <- ggplot( data = ranks.melt[ variable %in% c( 'hyspdisp.rank',
                                                         'sox.inverse.dist.rank')],
                     aes( x = SOx.rank.adj,
                          y = value)) +
  geom_point() +
  geom_abline( intercept = 0, 
               slope = 1,
               color = 'red') +
  facet_wrap( . ~ variable, ncol = 1) +#, labeller = as_labeller(hospital_names)) +
  #  ylab( "HyADS unit ranks") +
  xlab( "Adjoint unit ranks") +
  theme_bw() +
  theme( axis.title.y = element_blank(),
         #   strip.text = element_blank(),
         strip.background = element_blank())

ranks.hyspdisp.rank <- ranks.melt[ variable %in% c( 'hyspdisp.rank')]
ranks.hyspdisp.rank <- ranks.hyspdisp.rank[ order(value)]


















