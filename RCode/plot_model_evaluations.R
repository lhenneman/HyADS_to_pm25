source( '~/Dropbox/Harvard/RFMeval_Local/HyADS_to_pm25/RCode/hyads_to_pm25_functions.R')

## =================================================== ##
#   Load model objects
## =================================================== ##
load( '~/Dropbox/Harvard/RFMeval_Local/HyADS_to_pm25/RData/hyads_to_cmaq_models3.RData')

#coordinate reference system projection string for spatial data
p4s <- "+proj=lcc +lat_1=33 +lat_2=45 +lat_0=40 +lon_0=-97 +a=6370000 +b=6370000"

# get usa mask for masking
# download USA polygon from rnaturalearth
us_states.names <- state.abb[!(state.abb %in% c( 'HI', 'AK'))]
us_states <- st_transform( USAboundaries::us_states(), p4s)
mask.usa <- sf::as_Spatial(us_states)[ us_states$state_abbr %in% us_states.names,]

## =================================================== ##
#   Make plots - total fields
## =================================================== ##
# total fields in 2006
ggplot.a.raster( dats2006.a$cmaq.ddm,
                 preds.ann.hyads06w05$Y.ho.hat.raster$y.hat.gam.cv,
                 preds.ann.idwe06w05$Y.ho.hat.raster$y.hat.gam.cv,
                 ncol. = 3, facet.names = c( 'CMAQ-DDM Hybrid', 'HyADS*', 'IDWE*'),
                 legend.name = parse( text = 'list(2006~coal~source~impacts,mu*"g"~m^{-3})'),
                 bounds = c( 0, 4.5), mask.raster = mask.usa, 
                 theme.obj = theme( strip.text = element_text( size = 20)))

## =================================================== ##
#   Make plots - spatial bias fields
## =================================================== ##
# spatial bias
ggplot.a.raster( preds.ann.hyads05w06$Y.ho.hat.bias.raster$y.hat.lm.cv,
                 preds.ann.idwe05w06$Y.ho.hat.bias.raster$y.hat.lm.cv,
                 ncol. = 3, facet.names = c( 'HyADS*', 'IDWE*'),
                 legend.name = parse( text = 'list(Bias,mu*"g"~m^{-3})'),
                 bounds = c( -.5, .5), mask.raster = mask.usa, 
                 theme.obj = theme( strip.text = element_text( size = 20)))

## =================================================== ##
#   Make plots - covariate contributions
## =================================================== ##
# plots of annual covariate contributions
covs.all.hyads <- names( preds.ann.hyads06w05$Y.ho.terms.gam.raster)
covs.all.idwe  <- names( preds.ann.idwe06w05$Y.ho.terms.gam.raster)
covs.hyads      <- covs.all.hyads[grep( 'hyads', covs.all.hyads)]
covs.hyads.othr <- covs.all.hyads[ !(covs.all.hyads %in% c( covs.hyads, "s.x.y."))]
covs.idwe       <- covs.all.idwe[grep( 'idwe', covs.all.idwe)]
covs.idwe.othr  <- covs.all.idwe[ !(covs.all.idwe %in% c( covs.idwe, "s.x.y."))]


ggplot.a.raster( unstack( stack( sum( subset( preds.ann.hyads06w05$Y.ho.terms.gam.raster, covs.hyads)),
                                 sum( subset( preds.ann.hyads06w05$Y.ho.terms.gam.raster, covs.hyads.othr)),
                                 subset( preds.ann.hyads06w05$Y.ho.terms.gam.raster, 's.x.y.'))),
                 bounds = c( -4,4), ncol. = 3, mask.raster = mask.usa,
                 facet.names = c( 'HyADS', 'Other meteorology', 'Bivariate spline'))

ggplot.a.raster( unstack( stack( sum( subset( preds.ann.idwe06w05$Y.ho.terms.gam.raster, covs.idwe)),
                                 sum( subset( preds.ann.idwe06w05$Y.ho.terms.gam.raster, covs.idwe.othr)),
                                 subset( preds.ann.idwe06w05$Y.ho.terms.gam.raster, 's.x.y.'))),
                 bounds = c( -4,4), ncol. = 3, mask.raster = mask.usa,
                 facet.names = c( 'IDWE', 'Other meteorology', 'Bivariate spline'))


## =================================================== ##
#   Make plots - NME, ME, MB, etc
## =================================================== ##
metrics.ann <- rbind( preds.ann.hyads06w05$metrics[, mod.inpt := 'HyADS'],
                      preds.ann.idwe06w05$metrics[, mod.inpt := 'IDWE'])
metrics.ann.m <- melt( metrics.ann, id.vars = c( 'mod.name', 'mod.inpt'))
metrics.ann.m[ variable == 'R', variable := 'Pearson R']
metrics.ann.m[ variable == 'R.s', variable := 'Spearman R']

# get units and acronyms
mod.names <- data.table( mod.name = c( "lm.cv", "lm.ncv", "gam.cv", "adj.Z"),
                         model.names = c( "Linear", "Linear\nUncontrolled", 
                                          "GAM", "Z Score\nModel"))

metrics.ann.m <- merge( mod.names, metrics.ann.m, by = 'mod.name')

mod.eval <- ggplot( data = metrics.ann.m[ !(mod.name %in% c( 'adj.mean', 'adj.Z.only'))]) + 
  geom_hline( yintercept = 0) +
  geom_point( aes( x = model.names, y = value, 
                   color = mod.inpt, fill = mod.inpt),
              pch = 3, size = 3, stroke = 2, position = position_dodge( width = .4)) + 
  facet_wrap( . ~ variable, ncol = 1, scales = 'free_y') + 
  scale_y_continuous( expand = expansion( mult = c( .4))) + 
  theme_bw() + 
  theme( axis.title = element_blank(),
         axis.text.x = element_text( size = 12),
         axis.text.y = element_text( size = 12),
         legend.title = element_blank(),
         legend.text = element_text( size = 12),
         legend.position = 'bottom',
         panel.grid.major.x = element_line( size = 16),
         panel.grid.minor = element_blank(),
         strip.background = element_rect( fill = NA),
         strip.text = element_text( size = 14))

mod.eval
ggsave( file = '~/Dropbox/Harvard/Manuscripts/eval_local_exposure/Revision_JESEE/figures/cmaq-ddm_modelevals.png',
        mod.eval,
        width  = 6,
        height = 8 )


