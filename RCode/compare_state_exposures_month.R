source( '~/Dropbox/Harvard/RFMeval_Local/HyADS_to_pm25/RCode/hyads_to_pm25_functions.R')

## =================================================== ##
#   Read in monthly unit links
## =================================================== ##
loc <- '~/Dropbox/Harvard/RFMeval_Local/HyADS_to_pm25/RData/units_monthly/'
files.idwe  <- paste0( loc, 'popwgt_idwe', rep( c( 2006, 2011), each = 12), '_', 1:12, '.csv')
files.hyads <- paste0( loc, 'popwgt_hyads', rep( c( 2006, 2011), each = 12), '_', 1:12, '.csv')

idwe_month.dt <-  rbindlist( lapply( files.idwe,  fread, drop = c( 'V1', 'ID')))
hyads_month.dt <- rbindlist( lapply( files.hyads, fread, drop = c( 'V1', 'ID')))

# syntax/naming conventions
names.to.change <- c( 'mean_popwgt', 'popwgt', 'mean_pm')
setnames( idwe_month.dt,  names.to.change, paste0( names.to.change, '.idwe'))
setnames( hyads_month.dt, names.to.change, paste0( names.to.change, '.hyads'))

# merge idwe and month
names.to.link <- c( 'state_abbr', 'uID', 'pop_amnt', 'month')
units_month.dt <- merge( idwe_month.dt, hyads_month.dt, by = names.to.link)

## =================================================== ##
#   Do the evaluations
## =================================================== ##
# apply evaluation function
units_month.state <- units_month.dt[, evals.fn( popwgt.idwe, popwgt.hyads), 
                                    by = .( state_abbr, month)]

# melt
units_month.state.m <- melt( units_month.state, id.vars = c( 'state_abbr', 'month'),
                             variable.name = 'metric')

# pick out states, format months
states.use <- c( 'US', 'PA', 'KY', 'GA', 'WI', 'TX', 'CO', 'CA')
units_month.state.m <- units_month.state.m[ state_abbr %in% states.use]
units_month.state.m[, `:=` ( year.n = year( month),
                             month.n =  factor( month.name[ month( month)], levels = month.name),
                             state.factor = factor( state_abbr, levels = states.use))]


## =================================================== ##
#   plot
## =================================================== ##
save.loc <- '~/Dropbox/Harvard/Manuscripts/eval_local_exposure/Review_JESEE/figures/'
# ggplot( data = units_month.state.m[ state_abbr %in% states.use & 
#                                       year( month) == 2011]) + 
#   geom_line( aes( x = as.Date( month),
#                   y = value)) + 
#   facet_grid( metric ~ state_abbr, scales = 'free')

## split up into three plots
# 1. correlations (same plot for both)
# 2. NMB/NME
# 3. MB, RMSE


rect_bot <- seq( -2.75, 2.75, .5)
rect_bot2 <- seq( -3, 3, 1)
rect_bot3 <- seq( -3, 3, .04)
rect_bot4 <- seq( -3, 3, .02)
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
ggcors <- ggplot( data = units_month.state.m[metric %in% c( 'R.p', 'R.s')],
                  aes( x = month.n,
                       y =  value, #substring( value, 2)
                       color = metric,
                       group = metric)) +
  coord_cartesian( ylim = c( -.75,1)) +
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
  facet_grid( state.factor ~ year.n) +
  scale_color_viridis( begin = 0.4, 
                       end = 0.8,
                       discrete = T,
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
        legend.position = c( .185, .025),	
        legend.title = element_blank(),
        legend.text = element_text(size = 16),
        legend.box.background = element_rect(),
        legend.box.margin = margin(0, 0, 0, 0),
        legend.direction = 'horizontal',
        panel.grid = element_blank(),
        strip.background = element_blank(), #rect(fill='white'),
        strip.text = element_text(size=24))

ggsave( file = paste0( save.loc, 'correlations_monthly_2006_2011.png'),
        plot = ggcors,
        width =  7 * 1.3,
        height = 7 * 1.7)

# plot normalized mean error/bias
ggNMB <- ggplot( data = units_month.state.m[metric %in% c( 'NMB')],
                 aes( x = month.n,
                      y =  value, #substring( value, 2)
                      color = metric,
                      group = metric)) +
  coord_cartesian( ylim = c( -1,1)) +
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
  facet_grid( state.factor ~ year.n) +
  ylab( "NMB") +
  theme_bw() +
  theme(plot.title = element_text(size = 16, hjust = 0.5),
        axis.text = element_text( size = 14),
        axis.text.x = element_text( angle = 45,
                                    vjust = 1,
                                    hjust = 1),
        axis.title = element_text( size = 20),
        axis.title.x = element_blank(),
        legend.position = 'none',	
        panel.grid = element_blank(),
        strip.background = element_blank(), #rect(fill='white'),
        strip.text = element_text(size=24))

ggsave( file = paste0( save.loc, 'NMB_monthly_2006_2011.png'),
        plot = ggNMB,
        width =  7 * 1.3,
        height = 7 * 1.7)

# plot normalized mean error/bias
ggNME <- ggplot( data = units_month.state.m[metric %in% c( 'NME')],
                 aes( x = month.n,
                      y =  value, #substring( value, 2)
                      color = metric,
                      group = metric)) +
  coord_cartesian( ylim = c( 0,2)) +
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
  facet_grid( state.factor ~ year.n) +
  ylab( "NMB") +
  theme_bw() +
  theme(plot.title = element_text(size = 16, hjust = 0.5),
        axis.text = element_text( size = 14),
        axis.text.x = element_text( angle = 45,
                                    vjust = 1,
                                    hjust = 1),
        axis.title = element_text( size = 20),
        axis.title.x = element_blank(),
        legend.position = 'none',	
        panel.grid = element_blank(),
        strip.background = element_blank(), #rect(fill='white'),
        strip.text = element_text(size=24))

ggsave( file = paste0( save.loc, 'NME_monthly_2006_2011.png'),
        plot = ggNME,
        width =  7 * 1.3,
        height = 7 * 1.7)

# plot mean error and EMSE
ggME <- ggplot( data = units_month.state.m[metric %in% c( 'ME', 'RMSE')],
                aes( x = month.n,
                     y =  value, #substring( value, 2)
                     color = metric,
                     group = metric)) +
  coord_cartesian( ylim = c( 0,.064)) +
  geom_rect( data = rectangles3,
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
  scale_color_viridis( begin = 0.4, 
                       end = 0.8,
                       discrete = T) +
  facet_grid( state.factor ~ year.n) +
  ylab( parse( text = 'list(ME~and~RMSE,~mu*"g"~m^{-3})')) +
  theme_bw() +
  theme(plot.title = element_text(size = 16, hjust = 0.5),
        axis.text = element_text( size = 14),
        axis.text.x = element_text( angle = 45,
                                    vjust = 1,
                                    hjust = 1),
        axis.title = element_text( size = 20),
        axis.title.x = element_blank(),
        legend.position = c( .15, .09),	
        legend.title = element_blank(),
        legend.text = element_text(size = 16),
        legend.box.background = element_rect(),
        legend.box.margin = margin(0, 0, 0, 0),
        legend.direction = 'horizontal',
        panel.grid = element_blank(),
        strip.background = element_blank(), #rect(fill='white'),
        strip.text = element_text(size=24))

ggsave( file = paste0( save.loc, 'ME_RMSE_monthly_2006_2011.png'),
        plot = ggME,
        width =  7 * 1.3,
        height = 7 * 1.7)

# plot mean error and EMSE
ggMB <- ggplot( data = units_month.state.m[metric %in% c( 'MB')],
                aes( x = month.n,
                     y =  value, #substring( value, 2)
                     color = metric,
                     group = metric)) +
  coord_cartesian( ylim = c( -.024,.006)) +
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
  facet_grid( state.factor ~ year.n) +
  ylab( parse( text = 'list(MB,~mu*"g"~m^{-3})')) +
  theme_bw() +
  theme(plot.title = element_text(size = 16, hjust = 0.5),
        axis.text = element_text( size = 14),
        axis.text.x = element_text( angle = 45,
                                    vjust = 1,
                                    hjust = 1),
        axis.title = element_text( size = 20),
        axis.title.x = element_blank(),
        legend.position = 'none',	
         panel.grid = element_blank(),
        strip.background = element_blank(), #rect(fill='white'),
        strip.text = element_text(size=24))

ggsave( file = paste0( save.loc, 'MB_monthly_2006_2011.png'),
        plot = ggMB,
        width =  7 * 1.3,
        height = 7 * 1.7)


