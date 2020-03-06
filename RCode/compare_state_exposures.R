# copied from evaluate_RFMs/plot_unit_ranks

library( data.table)
library( ggplot2)
library( viridis)
library( grid)
library( ggpubr)
library( cowplot)
library( raster)
library( sf)
source( '~/Dropbox/Harvard/RFMeval_Local/HyADS_to_pm25/RCode/hyads_to_pm25_functions.R')
#coordinate reference system projection string for spatial data
p4s <- "+proj=lcc +lat_1=33 +lat_2=45 +lat_0=40 +lon_0=-97 +a=6370000 +b=6370000"


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
  # if( state.abbrev == 'US')
  #   state.abbrev <- state.abb
  state <- tolower( state.name[ match( state.abbrev, state.abb)])
  state.poly <- states[region %in% state]
  
  # -- make the plot  -- #
  if( is.null( geo.name)){
    plot.title <- NULL
  } else
    plot.title <- paste("Top facilities by pop-weighted exposure on", geo.name)
  
  # -- define columns for scalebar  -- #
  ranks.dt.units[, `:=` ( long = Longitude,
                          lat = Latitude)]
  
  # -- make the plot  -- #
  gg_coal <- ggplot() + 
    theme_bw() + 
    labs(title = plot.title) +
    facet_wrap( . ~ year, ncol = 2) +
    geom_polygon(data = states, 
                 aes(x = long, y = lat, group = group),
                 fill = 'white', 
                 color = "grey70",
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
                   color = size.var),
               stroke = 1,
               size = 3,
               shape = 21) +#,
    #      color = '#479ddd') + 
    geom_rect( data = latlonrange,
               aes(xmin = xlim[1] - 5, 
                   xmax = xlim[1] + (xlim[2] - xlim[1]) / 2, 
                   ymin = ylim[1] - 5, 
                   ymax = ylim[1] + .5), 
               fill = 'white', 
               color = NA) +
    scale_color_gradient( low = '#dd8747', 
                          high = '#4C061D', #'#6d213c',
                          guide = guide_colorbar(title.position = "left"),
                          name = bquote( atop( .( size.name), 'source impacts, µg'~m^{'-3'}))#,
                          # expression(paste('2012 total ', SO[2], ' emissions [tons]'))
                          # ylab( parse( text = 'list(PWSI["i,P"]^{m},~mu*"g"~m^{-3})')) +
                          # range = c(.5,5.0)
    ) +
    expand_limits(color = 0) +
    theme(
      plot.title = element_text(size = 16, hjust = 0.5), #element_blank(), #
      axis.title = element_text(size = 24),
      axis.text = element_blank(),
      axis.title.x = element_blank(),
      axis.title.y = element_blank(),
      legend.title = element_text(size = 16, vjust = 1, face = 'bold'),
      legend.title.align = 0,
      legend.position = "bottom", #c(.22, .15),
      legend.text = element_text(size = 14, angle = 30, vjust = .5),
      legend.background = element_rect(fill = 'transparent'),
      legend.key.size = unit(.04, 'npc'),
      legend.direction = 'horizontal',
      strip.text = element_text( size = 20, face = 'bold'),
      strip.background = element_blank( )
    )  +
    ggsn::scalebar( data = ranks.dt.units,
                    location = 'bottomleft',
                    anchor = c( x = latlonrange$xlim[1] - 1, y = latlonrange$ylim[1] + 2),
                    x.min = latlonrange$xlim[1],
                    y.min = latlonrange$ylim[1],
                    x.max = latlonrange$xlim[2],
                    y.max = latlonrange$ylim[2],
                    dist = dist.scalebar / 2, 
                    dist_unit = 'km',
                    height = 0.05, 
                    st.dist = 0.1, 
                    st.size = 4, 
                    box.fill = c( 'white', 'grey50'),
                    box.color = 'grey50',
                    transform = T,
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

# convert from ton/yr to kg / hr
adjoint_all[, sox_exp.adj := `SOx exposure` * 907.185 / 8760]

#======================================================================#
## area weight over states
#======================================================================#
# combine exposure and population
adjoint_all_pop <- merge( us_states.pop.dt.m, adjoint_all, 
                          by.x = c( 'state_abbr', 'year'), by.y = c( 'state', 'year'))
setnames( adjoint_all_pop, 'state_abbr', 'state')
adjoint_all_pop[, adj_popwgt := sox_exp.adj / pop] 

## =========================================================== ##
## adjoint variability between plume assumptions
## =========================================================== ##
adjoint_all_pop[, adj_popwgt.rank := frank( adj_popwgt), by = .( year, adj.name, state)]
ranks_US150.hys <- plot_ranked_facs( ranks.dt = adjoint_all_pop,
                                     rank.name = 'adj_popwgt.rank',
                                     size.var = 'adj_popwgt',
                                     toprank = 150,
                                     geo.name = "United States",
                                     dist.scalebar = 1000,
                                     latlonrange.year = 2006
)
ranks_US50.adj <- plot_ranked_facs( ranks.dt = adjoint_all_pop,
                                    state.abbrev = 'US',
                                    rank.name = 'adj_popwgt.rank',
                                    size.var = 'adj_popwgt',
                                    size.name = 'Adjoint population-weighted',
                                    toprank = 50,
                                    geo.name = NULL,
                                    dist.scalebar = 1000,
                                    latlonrange.year = 2006,
                                    xlims = c(-123, -68), #ranks_US150.hys$latlonrange$xlim,
                                    ylims = c(25, 49) #ranks_US150.hys$latlonrange$ylim
)

save.loc <- '~/Dropbox/Harvard/Manuscripts/eval_local_exposure/Revision_JESEE/figures/'
# ggsave( file = paste0( save.loc, 'PP_locs_2006_2011.png'),
#         plot = ranks_US50.adj$plot,
#         width =  11.5,
#         height = 4)

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
idwe_2006 <- fread( file.path( exp_home_annual, 'popwgt_idwe_annual_diff2006_3.csv'), drop = dropcols)
idwe_2011 <- fread( file.path( exp_home_annual, 'popwgt_idwe_annual_diff2011_3.csv'), drop = dropcols)
idwe_all <- rbind( idwe_2006, idwe_2011)[uID != 'Xtot.sum']
setnames( idwe_all, c( 'popwgt', 'mean_pm'), c( 'idwe.pw', 'idwe.mean'))

# HyADS
hyads_2006 <- fread( file.path( exp_home_annual, 'popwgt_hyads_annual_diff2006_3.csv'), drop = dropcols)
hyads_2011 <- fread( file.path( exp_home_annual, 'popwgt_hyads_annual_diff2011_3.csv'), drop = dropcols)
hyads_all <- rbind( hyads_2006, hyads_2011)
setnames( hyads_all, c( 'popwgt', 'mean_pm'), c( 'hyads.pw', 'hyads.mean'))

# raw IDWE
idwe_2006raw <- fread( file.path( exp_home_annual, 'popwgt_idwe_annual_raw2006_3.csv'), drop = dropcols)
idwe_2011raw <- fread( file.path( exp_home_annual, 'popwgt_idwe_annual_raw2011_3.csv'), drop = dropcols)
idwe_allraw <- rbind( idwe_2006raw, idwe_2011raw)[uID != 'Xtot.sum']
setnames( idwe_allraw, c( 'popwgt', 'mean_pm'), c( 'idwe.raw.pw', 'idwe.raw.mean'))

# raw HyADS
hyads_2006raw <- fread( file.path( exp_home_annual, 'popwgt_hyads_annual_raw2006_3.csv'), drop = dropcols)
hyads_2011raw <- fread( file.path( exp_home_annual, 'popwgt_hyads_annual_raw2011_3.csv'), drop = dropcols)
hyads_allraw <- rbind( hyads_2006raw, hyads_2011raw)
setnames( hyads_allraw, c( 'popwgt', 'mean_pm'), c( 'hyads.raw.pw', 'hyads.raw.mean'))

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
states.use <- c( 'US', 'CA', 'CO', 'TX', 'WI', 'KY', 'GA', 'PA')
rcm_exp <- Reduce(function(...) merge(..., all = TRUE, by = c( 'state_abbr', 'uID', 'year', 'ID')), 
                  list( idwe_all, hyads_all))[ state_abbr %in% states.use]
rcm_exp[, uID := gsub( '^X', '', uID)]
setnames( rcm_exp, 'state_abbr', 'state')

# combine raw, limit to input states
rcm_expraw <- Reduce(function(...) merge(..., all = TRUE, by = c( 'state_abbr', 'uID', 'year', 'ID')), 
                  list( idwe_allraw, hyads_allraw))[ state_abbr %in% states.use]
rcm_expraw[, uID := gsub( '^X', '', uID)]
setnames( rcm_expraw, 'state_abbr', 'state')


ggplot( data = rcm_exp[ state %in% states.use]) +
  geom_point( aes( x = hyads.pw, y = idwe.pw, color = year)) + 
  # scale_y_continuous( limits = c( 0,1)) + 
  facet_wrap( . ~ state)
## =========================================================== ##
## Merge adjoint and HyADS/dist
## =========================================================== ##
ranks_adj_all <- merge( adjoint_all_pop, rcm_exp,
                        by = c( "uID", "year", "state"),
                        all = T, allow.cartesian = TRUE)
ranks_adj_all[, state.factor := factor(`state`, levels =  states.use)]

# not looking at the ones with zero
ranks_adj_all <- ranks_adj_all[ uID != 'tot.sum'][!is.na( adj.name)]

# create dataset for sharing
ranks_adj_all_share.m <- ranks_adj_all[, .( uID, year, state, pop, 
                                          adj.name, adj_popwgt, idwe.pw, hyads.pw)]
ranks_adj_all_share.c <- dcast( ranks_adj_all_share.m, 
                                uID + year + state + pop + idwe.pw + hyads.pw ~ adj.name,
                                value.var = 'adj_popwgt')
fwrite( ranks_adj_all_share.c,
        '~/Dropbox/Harvard/Manuscripts/eval_local_exposure/CopyEdits_JESEE/Repo/dataset.csv')
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
                                                                    hyads.pw   = sum( hyads.pw)),
                                                                 by = c( 'state', 'year')]
rcm_statesums.m <- melt( rcm_statesums, id.vars = c( 'state', 'year', 'adj.pw'))
rcm_statesums.m[, bias := ( value - adj.pw) / adj.pw]

# load CMAQ data
# load( '~/Dropbox/Harvard/RFMeval_Local/HyADS_to_pm25/RData/hyads_to_cmaq_models3.RData')

# grid-weighted population as raster
# ddm2006.sf <- st_as_sf( rasterToPolygons( dats2006.a$cmaq.ddm))

# get total state populations
# us_states <- st_transform( USAboundaries::us_states(), p4s)
# ddm2006.states <- data.table( st_interpolate_aw( ddm2006.sf, us_states, extensive = F))

# merge together
# ddm2006.states[, `:=` ( year = 2006, state = us_states$state_abbr[Group.1], Group.1 = NULL, geometry = NULL)]
# rcm_statesums <- merge( rcm_statesums, ddm2006.states, by = c( 'state', 'year'))

# load in model fits
# modfits.annual <- fread( "~/Dropbox/Harvard/RFMeval_Local/HyADS_to_pm25/RData/annual_fields_hyads_idwe.csv", drop = 'V1')
# modfits.annual.c <- dcast( modfits.annual, x + y + year + model ~ field, value.var = 'layer')


## =========================================================== ##
## Mkae plots
## =========================================================== ##
ranks_adj_all.m <- melt( ranks_adj_all[, .( uID, year, state.factor, adj.name, adj_popwgt, hyads.pw, idwe.pw)],
                         id.vars = c( 'uID', 'year', 'state.factor', 'adj.name', 'adj_popwgt'))
ranks_adj_all.m[ variable == 'hyads.pw', var.name := 'HyADS']
ranks_adj_all.m[ variable == 'idwe.pw', var.name := 'IDWE']
states.use.vert <- c( 'US', 'PA', 'GA', 'KY', 'WI', 'TX', 'CO', 'CA')
ranks_adj_all.m[, state.factor := factor(`state.factor`, levels =  states.use.vert)]

gghyads.lin <- ggplot( data = ranks_adj_all.m[adj.name == 'initial'],
                       aes( x = adj_popwgt,
                            y = value,
                            color = var.name)) +
  geom_point( size = .5) +
  coord_cartesian( ylim = c( 0, 0.1)) +
  geom_smooth( method = 'lm', size = .5, fullrange = T, se = F) +
  scale_color_manual( values = c( '#479ddd', '#dd8747')) +
  facet_grid( state.factor ~ year) +
  ylab( parse( text = 'list(PWSI["i,P"]^{m},~mu*"g"~m^{-3})')) +
  xlab( parse( text = 'list(PWSI["i,P"]^{Adjoint},~mu*"g"~m^{-3})')) +
  theme_bw() + 
  theme( axis.text = element_text( size = 12),
         axis.title = element_text( size = 16),
         legend.direction = 'horizontal',
         legend.position = 'bottom',
         legend.text = element_text( size = 14),
         legend.title = element_blank(),
         strip.background = element_blank(),
         strip.text = element_text( size = 18))

gghyads.lin
save.loc <- '~/Dropbox/Harvard/Manuscripts/eval_local_exposure/Review_JESEE/figures/'
# ggsave( file = paste0( save.loc, 'scatter_annual_2006_2011.png'),
#         plot = gghyads.lin,
#         width =  7 * 1.3,
#         height = 7 * 1.7)


## =========================================================== ##
## Plot in the paper
## =========================================================== ##
eval.hyads <- ranks_adj_all[, eval.fn( hyads.pw, adj_popwgt, 'HyADS'), by = .( state.factor, year, adj.name)]
eval.idwe  <- ranks_adj_all[, eval.fn( idwe.pw, adj_popwgt, 'IDWE'), by = .( state.factor, year, adj.name)]
eval.all <- rbind( eval.hyads, eval.idwe)

# melt
eval.all.m <- melt( eval.all, id.vars = c( 'state.factor', 'year', 'adj.name', 'mod.name'),
                    variable.name = 'metric')

# redefine metric names
eval.all.m[metric == 'R.p', metric := 'Pearson~R']
eval.all.m[metric == 'R.s', metric := 'Spearman~R']

## define adjoint names, set up for plotting
adj.names <- data.table( adj.name = c( "initial", "layers_2-5", "stack_height", "stack_height_plus1", "stack_height_plus2"),
                         adj.name.adjust = c("Average", "Layers 2-5", "Stack Height", "Stack Height +1", "Stack Height +2"))
corrs.all <- merge( eval.all.m, adj.names, by = 'adj.name')
corrs.all[, state.factor := factor( state.factor, levels = states.use)]

# don't plot very high NME for CA and CO
cors.all.use <- corrs.all[ metric %in% c( 'Pearson~R', 'NMB', 'RMSE')]
corrs.removed <- cors.all.use[ state.factor %in% c( 'CA', 'CO') & value > 2]
cors.all.use[ state.factor %in% c( 'CA', 'CO') & value > 2, value := NA]
cors.all.use[ metric == 'NMB', value := value * 100]
cors.all.use[ metric == 'NMB', metric := 'NMB~"%"']
cors.all.use[ metric == 'RMSE', metric := 'RMSE~mu*"g"~m^{-3}']
gguse <- ggplot( data = cors.all.use,
                 aes( x = state.factor,
                      y = value,
                      color = mod.name,
                      shape = adj.name.adjust
                 )) +
  geom_hline( yintercept = 0) +
  geom_point( size = 4,
              position = position_dodge( width = .25)) + 
  scale_shape_manual( values = c( 0:2, 5:6)) +
  scale_color_manual( values = c( '#479ddd', '#dd8747')) +
  geom_vline( xintercept = 1.5,
              lty = 2) +
  facet_grid( metric ~ year, scales = 'free_y', labeller = label_parsed) +
  expand_limits( y = 0) +
  scale_y_continuous( expand = expansion( .1)) +
  guides(
    color = guide_legend(order = 1,
                         keywidth = 1.5),
    shape = guide_legend(order = 0,
                         keywidth = 1.5)
  ) +
  theme_bw() +
  theme( axis.text = element_text( size = 12),
         axis.title = element_blank( ),
         legend.direction = 'horizontal',
         legend.position = 'bottom',
         legend.text = element_text( size = 12),
         legend.title = element_blank(),
         panel.grid.major.x = element_line( size = 10,
                                            color = 'grey90'),
         panel.grid.minor = element_blank(),
         strip.background = element_blank(),
         strip.text = element_text( size = 18))

# ggsave( file = '~/Dropbox/Harvard/Manuscripts/eval_local_exposure/Revision_JESEE/figures/mixed_adjoints_heights_2006_2011.png',
#         gguse,
#         width  = 8 * 1.2,
#         height = 7 )


## =========================================================== ##
## Check out average distance metric
## =========================================================== ##
# read sox inverse distance and population, housekeeping
sox_inv_dist06 <- fread('~/Dropbox/Harvard/RFMeval_Local/HyADS_to_pm25/RData/ampd_dists_sox_weighted_2006_unit.csv', drop = 'V1')
pop06 <- fread( '~/Dropbox/Harvard/RFMeval_Local/HyADS_to_pm25/HyADS_grid/population/hyads_grid_population.csv', drop = c( 'V1', '2011'))
sox_inv_dist06.m <- melt( sox_inv_dist06[,!(names(sox_inv_dist06) %in% 'tot.sum'), with = F], id.vars = c( 'x', 'y'), variable.name = 'uID')
sox_inv_dist06.m[, uID := gsub( '^X', '', uID)]

#merge with units and population
units <- hyspdisp::units2006
sox_inv_dist06.m2 <- merge( units[, .(uID, SOx)], sox_inv_dist06.m, by = 'uID')
sox_inv_dist06.m2 <- merge( sox_inv_dist06.m2, pop06, by = c( 'x', 'y'))

#lculate distance, "Inf" dist interpreted as zero
sox_inv_dist06.m2[, dist := SOx / ( value)]
sox_inv_dist06.m2[ is.nan( dist), dist := 0]
sox_inv_dist06.m2[, numerator := SOx * dist * `2006`]
sox_inv_dist06.m2[, denominator := SOx * `2006`]

# sum by locations
sox_inv_dist06.dt <- sox_inv_dist06.m2[, .( Numerator = sum( numerator),
                                            Denominator = sum( denominator)), by = .( x, y)]

# assign states
# get usa mask for masking
# download USA polygon from rnaturalearth
us_states.names <- state.abb[!(state.abb %in% c( 'HI', 'AK'))]
us_states <- st_transform( USAboundaries::us_states(), p4s)
mask.usa <- sf::as_Spatial(us_states)[ us_states$state_abbr %in% us_states.names,]
mask.r <- rasterize( mask.usa[,'state_abbr'], Dpew.all)


# put back to dt and average by states
Dpew.all$ID <- mask.r$layer
mask.a <- levels( mask.r)[[1]]
Dpew.dt <- data.table( values( Dpew.all))[! is.na( ID)]
Dpew.dt <- merge( Dpew.dt, mask.a, by = 'ID')

# format for plotting, calculate Dpew
states.dpew <- c( 'US', 'CA', 'CO', 'TX', 'WI', 'GA', 'KY', 'PA')
ave_distances <- Dpew.dt[, .( avedist.pew = mean( layer, rm.na = T) / 1000), by = state_abbr]
ave_distances.all <- Dpew.dt[, .( avedist.pew = mean( layer, na.rm = T) / 1000, state_abbr = 'US')]
ave_distances <- rbind( ave_distances, ave_distances.all)[state_abbr %in% states.dpew]
ave_distances[, state.factor := factor(`state_abbr`, levels = states.dpew)]

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




gg_combine <- plot_grid(gguse, gg_avedist,
                        rel_heights = c( 4, 1),
                        labels = NULL, ncol = 1, 
                        align = 'hv', axis = 'lr') 

ggsave( file = '~/Dropbox/Harvard/Manuscripts/eval_local_exposure/Revision_JESEE/figures/eval_dists_2006_2011.png',
        gg_combine,
        width  = 8 * 1.3,
        height = 5 * 1.5)


## =========================================================== ##
# plot Dpew
## =========================================================== ##
sox_inv_dist06.r <- rasterFromXYZ( sox_inv_dist06.dt)
plot( sox_inv_dist06.r$Numerator / sox_inv_dist06.r$Denominator)
Dpew.all <- sox_inv_dist06.r$Numerator / sox_inv_dist06.r$Denominator
Dpes.all <- focal( Dpew.all, w = matrix( 1, 3, 3), fun = mean, NAonly = T, na.rm = T)

# state shapefile
mask.usa.st <- data.table( st_as_sf(mask.usa))


# make the plot
gg.dpew <- ggplot.a.raster( Dpes.all / 1000, facet.names = "layer", mask.raster = mask.usa,
                 legend.name = expression( D^{'pew'}~'km'),
                 theme.obj = theme( legend.key.width = unit( 0.1, 'npc'),
                                    strip.text = element_blank())) + 
  scale_fill_viridis( direction = -1,
                      name = expression( D^{'pew'}~'[km]')) + 
  geom_sf( data = mask.usa.st, aes( geometry = geometry), fill = NA, color = 'grey70') +
  geom_sf_label( data = mask.usa.st[state_abbr %in% states.dpew], aes( label = state_abbr,
                                                                       geometry = geometry))

  
ggsave( file = '~/Dropbox/Harvard/Manuscripts/eval_local_exposure/Revision_JESEE/figures/Dpew.png',
        gg.dpew,
        width  = 7,
        height = 5)



## =========================================================== ##
## Plot in the SI
## =========================================================== ##
cors.all.SI <- corrs.all[ metric %in% c( 'Spearman~R', 'NME', 'MB')]
corrs.removed.SI <- cors.all.SI[ state.factor %in% c( 'CA', 'CO') & value > 2]
cors.all.SI[ state.factor %in% c( 'CA', 'CO') & value > 2, value := NA]
cors.all.SI[ metric == 'NME', value := value * 100]
cors.all.SI[ metric == 'NME', metric := 'NME~"%"']
cors.all.SI[ metric == 'MB', metric := 'MB~mu*"g"~m^{-3}']

ggSI <- ggplot( data = cors.all.SI,
                aes( x = state.factor,
                     y = value,
                     color = mod.name,
                     shape = adj.name.adjust
                )) +
  geom_hline( yintercept = 0) +
  geom_point( size = 4,
              position = position_dodge( width = .25)) + 
  scale_shape_manual( values = c( 0:2, 5:6)) +
  scale_color_manual( values = c( '#479ddd', '#dd8747')) +
  geom_vline( xintercept = 1.5,
              lty = 2) +
  facet_grid( metric ~ year, scales = 'free_y', labeller = label_parsed) +
  expand_limits( y = 0) +
  scale_y_continuous( expand = expand_scale( mult = c( .1))) +
  guides(
    color = guide_legend(order = 1,
                         keywidth = 1.5),
    shape = guide_legend(order = 0,
                         keywidth = 1.5)
  ) +
  theme_bw() +
  theme( axis.text = element_text( size = 12),
         axis.title = element_blank( ),
         legend.direction = 'horizontal',
         legend.position = 'bottom',
         legend.text = element_text( size = 12),
         legend.title = element_blank(),
         panel.grid.major.x = element_line( size = 10,
                                            color = 'grey90'),
         panel.grid.minor = element_blank(),
         strip.background = element_blank(),
         strip.text = element_text( size = 18))
ggSI
ggsave( file = '~/Dropbox/Harvard/Manuscripts/eval_local_exposure/Revision_JESEE/figures/SI_adjoints_heights_2006_2011.png',
        ggSI,
        width  = 8 * 1.2,
        height = 6 )

## =========================================================== ##
## Old figures
## =========================================================== ##
ggNM <- ggplot( data = corrs.all[ metric %in% c( 'NME', 'NMB')],
                aes( x = state.factor,
                     y = value,
                     color = mod.name,
                     shape = adj.name.adjust
                )) +
  scale_y_continuous( limits = c( NA, 1), labels = scales::percent_format(accuracy = 1),
                      expand = expand_scale( mult = c( .1))) + 
  geom_hline( yintercept = 0) +
  geom_point( size = 4,
              position = position_dodge( width = .25)) + 
  scale_shape_manual( values = c( 0:2, 5:6)) +
  scale_color_manual( values = c( '#479ddd', '#dd8747')) +
  geom_vline( xintercept = 1.5,
              lty = 2) +
  facet_grid( metric ~ year, scales = 'free_y') +
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
         axis.title = element_blank( ),
         legend.direction = 'horizontal',
         legend.position = 'bottom',
         legend.text = element_text( size = 12),
         legend.title = element_blank(),
         panel.grid.major.x = element_line( size = 10,
                                            color = 'grey90'),
         panel.grid.minor = element_blank(),
         strip.background = element_blank(),
         strip.text = element_text( size = 18))

ggsave( file = '~/Dropbox/Harvard/Manuscripts/eval_local_exposure/Review_JESEE/figures/NM_adjoints_heights_2006_2011.png',
        ggNM,
        width  = 8 * 1.2,
        height = 5 )

ggMB <- ggplot( data = corrs.all[ metric %in% c( 'MB', 'RMSE')],
                aes( x = state.factor,
                     y = value,
                     color = mod.name,
                     shape = adj.name.adjust
                )) +
  scale_y_continuous( limits = c( NA, NA), expand = expand_scale( mult = c( .1))) + 
  geom_hline( yintercept = 0) +
  geom_point( size = 4,
              position = position_dodge( width = .25)) + 
  scale_shape_manual( values = c( 0:2, 5:6)) +
  scale_color_manual( values = c( '#479ddd', '#dd8747')) +
  geom_vline( xintercept = 1.5,
              lty = 2) +
  facet_grid( metric ~ year, scales = 'free_y') +
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
         axis.title = element_blank( ),
         legend.direction = 'horizontal',
         legend.position = 'bottom',
         legend.text = element_text( size = 12),
         legend.title = element_blank(),
         panel.grid.major.x = element_line( size = 10,
                                            color = 'grey90'),
         panel.grid.minor = element_blank(),
         strip.background = element_blank(),
         strip.text = element_text( size = 18))

ggsave( file = '~/Dropbox/Harvard/Manuscripts/eval_local_exposure/Review_JESEE/figures/MB_adjoints_heights_2006_2011.png',
        ggMB,
        width  = 8 * 1.2,
        height = 5 )


corrs.all[adj.name == 'initial' & state.factor == 'US']
corrs.all[adj.name == 'initial' & state %in% c( 'KY', 'PA', 'GA')]

summary( corrs.all[state.factor == 'US' & metric == 'NMB' & mod.name == 'HyADS'])
summary( corrs.all[state.factor == 'US' & metric == 'RMSE' & mod.name == 'HyADS'])
summary( corrs.all[state.factor == 'US' & metric == 'NMB' & mod.name == 'IDWE'])
summary( corrs.all[state.factor == 'US' & metric == 'RMSE' & mod.name == 'IDWE'])
summary( corrs.all[state.factor %in% c( 'KY', 'PA', 'GA') & metric == 'NMB' & mod.name == 'IDWE'])
summary( corrs.all[state.factor %in% c( 'KY', 'PA', 'GA') & metric == 'NMB' & mod.name == 'HyADS'])
summary( corrs.all[state.factor %in% c( 'KY', 'PA', 'GA') & metric == 'RMSE' & mod.name == 'HyADS'])
summary( corrs.all[state.factor %in% c( 'KY', 'PA', 'GA') & metric == 'RMSE' & mod.name == 'IDWE'])


## =========================================================== ##
## correlations of raw hyads/idwe
## =========================================================== ##
ranks_adj_allraw <- merge( adjoint_all_pop, rcm_expraw,
                        by = c( "uID", "year", "state"),
                        all = T, allow.cartesian = TRUE)
ranks_adj_allraw[, state.factor := factor(`state`, levels =  states.use)]

# not looking at the ones with zero
ranks_adj_allraw <- ranks_adj_allraw[ uID != 'tot.sum'][!is.na( adj.name)]

eval.hyadsraw <- ranks_adj_allraw[, eval.fn( hyads.raw.pw, adj_popwgt, 'HyADS'), by = .( state.factor, year, adj.name)]
eval.idweraw  <- ranks_adj_allraw[, eval.fn( idwe.raw.pw, adj_popwgt, 'IDWE'), by = .( state.factor, year, adj.name)]
eval.allraw <- rbind( eval.hyadsraw, eval.idweraw)

# melt
eval.all.m <- melt( eval.allraw, id.vars = c( 'state.factor', 'year', 'adj.name', 'mod.name'),
                    variable.name = 'metric')

# redefine metric names
eval.all.m[metric == 'R.p', metric := 'Pearson~R']
eval.all.m[metric == 'R.s', metric := 'Spearman~R']

## define adjoint names, set up for plotting
adj.names <- data.table( adj.name = c( "initial", "layers_2-5", "stack_height", "stack_height_plus1", "stack_height_plus2"),
                         adj.name.adjust = c("Average", "Layers 2-5", "Stack Height", "Stack Height +1", "Stack Height +2"))
corrs.all <- merge( eval.all.m, adj.names, by = 'adj.name')
corrs.all[, state.factor := factor( state.factor, levels = states.use)]

# don't plot very high NME for CA and CO
cors.all.use <- corrs.all[ metric %in% c( 'Pearson~R', 'Spearman~R')]
gguseraw <- ggplot( data = cors.all.use,
                 aes( x = state.factor,
                      y = value,
                      color = mod.name,
                      shape = adj.name.adjust
                 )) +
  geom_hline( yintercept = 0) +
  geom_point( size = 4,
              position = position_dodge( width = .25)) + 
  scale_shape_manual( values = c( 0:2, 5:6)) +
  scale_color_manual( values = c( '#479ddd', '#dd8747')) +
  geom_vline( xintercept = 1.5,
              lty = 2) +
  facet_grid( metric ~ year, scales = 'free_y', labeller = label_parsed) +
  expand_limits( y = 0) +
  scale_y_continuous( expand = expansion( .1)) +
  guides(
    color = guide_legend(order = 1,
                         keywidth = 1.5),
    shape = guide_legend(order = 0,
                         keywidth = 1.5)
  ) +
  theme_bw() +
  theme( axis.text = element_text( size = 12),
         axis.title = element_blank( ),
         legend.direction = 'horizontal',
         legend.position = 'bottom',
         legend.text = element_text( size = 12),
         legend.title = element_blank(),
         panel.grid.major.x = element_line( size = 10,
                                            color = 'grey90'),
         panel.grid.minor = element_blank(),
         strip.background = element_blank(),
         strip.text = element_text( size = 18))

ggsave( file = '~/Dropbox/Harvard/Manuscripts/eval_local_exposure/Revision_JESEE/figures/raw_correlations_2006_2011.png',
        gguseraw,
        width  = 8 * 1.2,
        height = 7 * .8)




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


















