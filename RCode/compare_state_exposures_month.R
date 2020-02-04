source( '~/Dropbox/Harvard/RFMeval_Local/HyADS_to_pm25/RCode/hyads_to_pm25_functions.R')
save.loc <- '~/Dropbox/Harvard/Manuscripts/eval_local_exposure/Review_JESEE/figures/'

## =================================================== ##
#   Monthly evaluations of original predictions
## =================================================== ##
load( '~/Dropbox/Harvard/RFMeval_Local/HyADS_to_pm25/RData/hyads_to_cmaq_models3.RData')

preds.metrics.hyads <- preds.mon.hyads06w05[ 'metrics',]
preds.metrics.idwe  <- preds.mon.idwe06w05[ 'metrics',]

metrics <- data.table( month = c( as.Date( gsub( '\\.', '-', gsub( 'X', '', names( preds.metrics.hyads)))),
                                  as.Date( gsub( '\\.', '-', gsub( 'X', '', names( preds.metrics.idwe))))),
                       model = c( rep( 'HyADS', length( names( preds.metrics.hyads))),
                                  rep( 'IDWE', length( names( preds.metrics.idwe)))),
                       class = c( rep( 'GAM', 2 * length( names( preds.metrics.hyads))),
                                  rep( 'Linear', 2 * length( names( preds.metrics.idwe)))),
                       'Pearson~R' = c( sapply( preds.metrics.hyads, function( dt) dt[ mod.name == 'gam.cv']$R),
                                        sapply( preds.metrics.idwe, function( dt) dt[ mod.name == 'gam.cv']$R),
                                        sapply( preds.metrics.hyads, function( dt) dt[ mod.name == 'lm.cv']$R),
                                        sapply( preds.metrics.idwe, function( dt) dt[ mod.name == 'lm.cv']$R)),
                       'NMB~"%"' = c( sapply( preds.metrics.hyads, function( dt) dt[ mod.name == 'gam.cv']$NMB*100),
                                      sapply( preds.metrics.idwe, function( dt) dt[ mod.name == 'gam.cv']$NMB*100),
                                      sapply( preds.metrics.hyads, function( dt) dt[ mod.name == 'lm.cv']$NMB*100),
                                      sapply( preds.metrics.idwe, function( dt) dt[ mod.name == 'lm.cv']$NMB*100)),
                       'NME~"%"' = c( sapply( preds.metrics.hyads, function( dt) dt[ mod.name == 'gam.cv']$NME*100),
                                      sapply( preds.metrics.idwe, function( dt) dt[ mod.name == 'gam.cv']$NME*100),
                                      sapply( preds.metrics.hyads, function( dt) dt[ mod.name == 'lm.cv']$NME*100),
                                      sapply( preds.metrics.idwe, function( dt) dt[ mod.name == 'lm.cv']$NME*100)),
                       'MB~mu*"g"~m^{-3}' = c( sapply( preds.metrics.hyads, function( dt) dt[ mod.name == 'gam.cv']$MB),
                                               sapply( preds.metrics.idwe, function( dt) dt[ mod.name == 'gam.cv']$MB),
                                               sapply( preds.metrics.hyads, function( dt) dt[ mod.name == 'lm.cv']$MB),
                                               sapply( preds.metrics.idwe, function( dt) dt[ mod.name == 'lm.cv']$MB)),
                       'RMSE~mu*"g"~m^{-3}' = c( sapply( preds.metrics.hyads, function( dt) dt[ mod.name == 'gam.cv']$RMSE),
                                                 sapply( preds.metrics.idwe, function( dt) dt[ mod.name == 'gam.cv']$RMSE),
                                                 sapply( preds.metrics.hyads, function( dt) dt[ mod.name == 'lm.cv']$RMSE),
                                                 sapply( preds.metrics.idwe, function( dt) dt[ mod.name == 'lm.cv']$RMSE)))
metrics.m <- melt( metrics, id.vars = c( 'model', 'month', 'class'), variable.name = 'metric')

gg_monthlyevals <- ggplot( data = metrics.m,
                           aes( x = month, y = value, lty = class, color = model)) + 
  geom_hline( yintercept = 0) +
  geom_line() + geom_point() + 
  facet_wrap( . ~ metric, scales = 'free_y', ncol = 1, labeller = label_parsed) + 
  scale_x_date( date_labels = '%B', breaks = unique( metrics.m$month)) +
  scale_color_viridis( begin = 0.4, end = 0.8, discrete = T, option = 'A') +
  expand_limits( y = 0) + 
  theme_bw() +
  theme(plot.title = element_text(size = 16, hjust = 0.5),
        axis.text = element_text( size = 14),
        axis.text.x = element_text( angle = 45,
                                    vjust = 1,
                                    hjust = 1),
        axis.title = element_blank(),
        legend.position = 'bottom',	
        legend.title = element_blank(),
        legend.text = element_text(size = 16),
        legend.box.background = element_rect(),
        legend.box.margin = margin(0, 0, 0, 0),
        legend.direction = 'horizontal',
        panel.grid.minor.x = element_blank(),
        strip.background = element_blank(), #rect(fill='white'),
        strip.text = element_text(size=24))
ggsave( file = paste0( save.loc, 'monthly_modelevals_2006_2011.png'),
        plot = gg_monthlyevals,
        width =  7 * 1.3,
        height = 7 * 1.7)


## =================================================== ##
#   Read in monthly unit links
## =================================================== ##
loc <- '~/Dropbox/Harvard/RFMeval_Local/HyADS_to_pm25/RData/units_monthly/'
files.idwe.gam  <- paste0( loc, 'popwgt_idwe', rep( c( 2006, 2011), each = 12), '_', 1:12, '_3g.csv')
files.hyads.gam <- paste0( loc, 'popwgt_hyads', rep( c( 2006, 2011), each = 12), '_', 1:12, '_3g.csv')
files.idwe.lm  <- paste0( loc, 'popwgt_idwe', rep( c( 2006, 2011), each = 12), '_', 1:12, '_3.csv')
files.hyads.lm <- paste0( loc, 'popwgt_hyads', rep( c( 2006, 2011), each = 12), '_', 1:12, '_3.csv')

idwe_month_gam.dt <-  rbindlist( lapply( files.idwe.gam,  fread, drop = c( 'V1', 'ID')))[, `:=` (model = 'Model 6', name = 'IDWE')]
hyads_month_gam.dt <- rbindlist( lapply( files.hyads.gam, fread, drop = c( 'V1', 'ID')))[, `:=` (model = 'Model 6', name = 'HyADS')]
idwe_month_lm.dt <-  rbindlist( lapply( files.idwe.lm,  fread, drop = c( 'V1', 'ID')))[, `:=` (model = 'Model 5', name = 'IDWE')]
hyads_month_lm.dt <- rbindlist( lapply( files.hyads.lm, fread, drop = c( 'V1', 'ID')))[, `:=` (model = 'Model 5', name = 'HyADS')]

# syntax/naming conventions
units_month.dt <- rbind( idwe_month_gam.dt, hyads_month_gam.dt, idwe_month_lm.dt, hyads_month_lm.dt)
units_month.dt.c <- na.omit( dcast( units_month.dt, state_abbr + uID + month + model ~ name, value.var = 'popwgt'))

## =================================================== ##
#   Do the evaluations
## =================================================== ##
# apply evaluation function
units_month.state <- units_month.dt.c[, evals.fn( IDWE, HyADS), 
                                      by = .( state_abbr, month, model)]

# melt
units_month.state.m <- melt( units_month.state, id.vars = c( 'state_abbr', 'month', 'model'),
                             variable.name = 'metric')

# pick out states, format months
states.use <- c( 'US', 'PA', 'KY', 'GA', 'WI', 'TX', 'CO', 'CA')
units_month.state.m <- units_month.state.m[ state_abbr %in% states.use]
units_month.state.m[, `:=` ( year.n = year( month),
                             month.n =  factor( month.name[ month( month)], levels = month.name),
                             state.factor = factor( state_abbr, levels = states.use))]

# negative correlations in CO in 3 months come from negative 
# results predictions in GAM, which lead to inverted correlation
# == higher IDWE means less PM2.5
## =================================================== ##
#   plot
## =================================================== ##
# ggplot( data = units_month.state.m[ state_abbr %in% states.use & 
#                                       year( month) == 2011]) + 
#   geom_line( aes( x = as.Date( month),
#                   y = value)) + 
#   facet_grid( metric ~ state_abbr, scales = 'free')

## split up into three plots
# 1. correlations (same plot for both)
# 2. NMB/NME
# 3. MB, RMSE


rect_bot <- seq( -10.75, 10.75, .5)
rect_bot2 <- seq( -30, 30, 1)
rect_bot3 <- seq( -30, 30, .04)
rect_bot4 <- seq( -30, 30, .02)
rectangles <- data.frame(
  xmin = rect_bot,
  xmax = rect_bot + .25,
  ymin = -Inf,
  ymax = +Inf,
  fill = 'grey50'
)
rectangles2 <- data.frame(
  xmin = rect_bot2,
  xmax = rect_bot2 + .5,
  ymin = -Inf,
  ymax = +Inf,
  fill = 'grey50'
)
rectangles3 <- data.frame(
  xmin = rect_bot3,
  xmax = rect_bot3 + .02,
  ymin = -Inf,
  ymax = +Inf,
  fill = 'grey50'
)
rectangles4 <- data.frame(
  xmin = rect_bot4,
  xmax = rect_bot4 + .01,
  ymin = -Inf,
  ymax = +Inf,
  fill = 'grey50'
)

# plot correlations across the year
ggcors <- ggplot( data = units_month.state.m[metric %in% c( 'R.p', 'R.s') & model == 'Model 5'],
                  aes( x = as.Date( month),
                       y =  value, #substring( value, 2)
                       color = metric)) +
  coord_cartesian( ylim = c( 0,1)) +
  scale_x_date( date_labels = '%B', breaks = unique( as.Date( units_month.state.m$month))) +
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
  facet_grid( state.factor ~ year.n, scales = 'free_x') +
  scale_color_viridis( begin = 0.4, end = 0.8, discrete = T,
                       breaks = c( 'R.p', 'R.s'),
                       labels = c( 'Pearson', 'Spearman')) +
  ylab( "Correlations") +
  xlab( "Month") +
  theme_bw() +
  theme(plot.title = element_text(size = 16, hjust = 0.5),
        axis.text = element_text( size = 14),
        axis.text.x = element_text( angle = 45,
                                    vjust = 1,
                                    hjust = 1),
        axis.title = element_text( size = 20),
        axis.title.x = element_blank(),
        legend.position = 'bottom', #c( .185, .025),	
        legend.title = element_blank(),
        legend.text = element_text(size = 16),
        legend.box.background = element_rect(),
        legend.box.margin = margin(0, 0, 0, 0),
        legend.direction = 'horizontal',
        panel.grid = element_blank(),
        strip.background = element_blank(), #rect(fill='white'),
        strip.text = element_text(size=24))
ggcors
ggsave( file = paste0( save.loc, 'correlations_monthly_2006_2011.png'),
        plot = ggcors,
        width =  7 * 1.3,
        height = 7 * 1.7)

# plot normalized mean error/bias
ggNMB <- ggplot( data = units_month.state.m[metric %in% c( 'NMB') & model == 'Model 5'],
                 aes( x = as.Date( month),
                      y =  value)) +
  coord_cartesian( ylim = c( -1,3.5)) +
  scale_x_date( date_labels = '%B', breaks = unique( as.Date( units_month.state.m$month))) +
  scale_y_continuous( labels = scales::percent_format(accuracy = 1)) +
  geom_rect( data = rectangles2,
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
  facet_grid( state.factor ~ year.n, scales = 'free') +
  ylab( "NMB") +
  theme_bw() +
  theme(plot.title = element_text(size = 16, hjust = 0.5),
        axis.text = element_text( size = 14),
        axis.text.x = element_text( angle = 45,
                                    vjust = 1,
                                    hjust = 1),
        axis.title = element_text( size = 20),
        axis.title.x = element_blank(),
        legend.position = 'bottom', #c( .185, .025),	
        legend.title = element_blank(),
        legend.text = element_text(size = 16),
        legend.box.background = element_rect(),
        legend.box.margin = margin(0, 0, 0, 0),
        legend.direction = 'horizontal',
        panel.grid.major.x = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(), #rect(fill='white'),
        strip.text = element_text(size=24))

ggsave( file = paste0( save.loc, 'NMB_monthly_2006_2011.png'),
        plot = ggNMB,
        width =  7 * 1.3,
        height = 7 * 1.7)
summary( units_month.state.m[metric %in% c( 'NMB') & model == 'Model 5' & state_abbr == 'CA'])
summary( units_month.state.m[metric %in% c( 'NMB') & model == 'Model 5' & state_abbr == 'CO'])
summary( units_month.state.m[metric %in% c( 'NMB') & model == 'Model 5' & state_abbr == 'US'])
summary( units_month.state.m[metric %in% c( 'NMB') & model == 'Model 5' & state_abbr == 'PA'])
summary( units_month.state.m[metric %in% c( 'NMB') & model == 'Model 5' & state_abbr == 'GA'])
summary( units_month.state.m[metric %in% c( 'NMB') & model == 'Model 5' & state_abbr == 'KY'])

# plot normalized mean error/bias
ggNME <- ggplot( data = units_month.state.m[metric %in% c( 'NME') & model == 'Model 5'],
                 aes( x = as.Date( month),
                      y =  value)) +
  coord_cartesian( ylim = c( 0,4)) +
  scale_x_date( date_labels = '%B', breaks = unique( as.Date( units_month.state.m$month))) +
  scale_y_continuous( labels = scales::percent_format(accuracy = 1)) +
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
  facet_grid( state.factor ~ year.n, scales = 'free_x') +
  ylab( "NME") +
  theme_bw() +
  theme(plot.title = element_text(size = 16, hjust = 0.5),
        axis.text = element_text( size = 14),
        axis.text.x = element_text( angle = 45,
                                    vjust = 1,
                                    hjust = 1),
        axis.title = element_text( size = 20),
        axis.title.x = element_blank(),
        legend.position = 'bottom', #c( .185, .025),	
        legend.title = element_blank(),
        legend.text = element_text(size = 16),
        legend.box.background = element_rect(),
        legend.box.margin = margin(0, 0, 0, 0),
        legend.direction = 'horizontal',
        panel.grid = element_blank(),
        strip.background = element_blank(), #rect(fill='white'),
        strip.text = element_text(size=24))

ggsave( file = paste0( save.loc, 'NME_monthly_2006_2011.png'),
        plot = ggNME,
        width =  7 * 1.3,
        height = 7 * 1.7)
summary( units_month.state.m[metric %in% c( 'NME') & model == 'Model 5' & state_abbr == 'CA'])
summary( units_month.state.m[metric %in% c( 'NME') & model == 'Model 5' & state_abbr == 'CO'])
# plot mean error and EMSE
ggME <- ggplot( data = units_month.state.m[metric %in% c( 'ME', 'RMSE') & model == 'Model 5'],
                aes( x = as.Date( month),
                     y =  value#, #substring( value, 2)
                )) +
  coord_cartesian( ylim = c( 0,.099)) +
  scale_x_date( date_labels = '%B', breaks = unique( as.Date( units_month.state.m$month))) +
  scale_y_continuous( breaks = seq( 0,.08,.04)) +
  geom_rect( data = rectangles3,
             aes( ymin = xmin,
                  ymax = xmax),
             fill = 'grey95',
             xmin = -Inf,
             xmax = +Inf,
             inherit.aes = FALSE) +
  geom_hline( yintercept = 0,
              size = 0.5,
              color = 'grey50') + 
  geom_line( aes( color = metric), size = 1.5) +
  facet_grid( state.factor ~ year.n, scales = 'free_x') +
  scale_color_viridis( begin = 0.4, end = 0.8, discrete = T) +
  ylab( parse( text = 'list(ME~and~RMSE,mu*"g"~m^{-3})')) +
  xlab( "Month") +
  theme_bw() +
  theme(plot.title = element_text(size = 16, hjust = 0.5),
        axis.text = element_text( size = 14),
        axis.text.x = element_text( angle = 45,
                                    vjust = 1,
                                    hjust = 1),
        axis.title = element_text( size = 20),
        axis.title.x = element_blank(),
        legend.position = 'bottom', #c( .185, .025),	
        legend.title = element_blank(),
        legend.text = element_text(size = 16),
        legend.box.background = element_rect(),
        legend.box.margin = margin(0, 0, 0, 0),
        legend.direction = 'horizontal',
        panel.grid = element_blank(),
        strip.background = element_blank(), #rect(fill='white'),
        strip.text = element_text(size=24))
ggME
ggsave( file = paste0( save.loc, 'ME_RMSE_monthly_2006_2011.png'),
        plot = ggME,
        width =  7 * 1.3,
        height = 7 * 1.7)

# plot mean error and EMSE
ggMB <- ggplot( data = units_month.state.m[metric %in% c( 'MB') & model == 'Model 5'],
                aes( x = as.Date( month),
                     y =  value)) +
  coord_cartesian( ylim = c( -.025,.03)) +
  scale_x_date( date_labels = '%B', breaks = unique( as.Date( units_month.state.m$month))) +
  scale_y_continuous( breaks = seq( -.02, .02, .02)) +
  geom_rect( data = rectangles4,
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
  facet_grid( state.factor ~ year.n, scales = 'free_x') +
  ylab( parse( text = 'list(MB,~mu*"g"~m^{-3})')) +
  theme_bw() +
  theme(plot.title = element_text(size = 16, hjust = 0.5),
        axis.text = element_text( size = 14),
        axis.text.x = element_text( angle = 45,
                                    vjust = 1,
                                    hjust = 1),
        axis.title = element_text( size = 20),
        axis.title.x = element_blank(),
        legend.position = 'bottom', #c( .185, .025),	
        legend.title = element_blank(),
        legend.text = element_text(size = 16),
        legend.box.background = element_rect(),
        legend.box.margin = margin(0, 0, 0, 0),
        legend.direction = 'horizontal',
        panel.grid = element_blank(),
        strip.background = element_blank(), #rect(fill='white'),
        strip.text = element_text(size=24))
ggMB
ggsave( file = paste0( save.loc, 'MB_monthly_2006_2011.png'),
        plot = ggMB,
        width =  7 * 1.3,
        height = 7 * 1.7)


## =================================================== ##
#   investigate
## =================================================== ##
summary( units_month.state.m[metric %in% c( 'R.p') & model == 'Model 5' & state_abbr == 'US'])
summary( units_month.state.m[metric %in% c( 'NME') & model == 'Model 5' & state_abbr == 'US'])
summary( units_month.state.m[metric %in% c( 'NME') & model == 'Model 5' & state_abbr == 'PA'])
summary( units_month.state.m[metric %in% c( 'NME') & model == 'Model 5' & state_abbr == 'KY'])
summary( units_month.state.m[metric %in% c( 'NME') & model == 'Model 5' & state_abbr == 'GA'])
summary( units_month.state.m[metric %in% c( 'R.s') & model == 'Model 5' & state_abbr == 'CA'])

