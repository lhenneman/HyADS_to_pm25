## ======================================== ##
##      Codechunk #1
## ======================================== ##

library( hyspdisp)
library( ggplot2)
library( viridis)
library( cowplot)

load( '~/Dropbox/Harvard/RFMeval_Local/HyADS_to_PM25/RData/cmaqddm_to_pm25_crossval.RData')

# cmaq
ggplot.a.raster( dats2005.s$cmaq.ddm, dats2006.s$cmaq.ddm, bounds = c( 0, 6),
                 facet.names = c( 'Hybrid CMAQ-DDM - 2005', 'Hybrid CMAQ-DDM - 2006'))

#hyads
ggplot.a.raster( dats2005.s$hyads, dats2006.s$hyads, bounds = c( 0, NA),
                 facet.names = c( 'HyADS - 2005', 'HyADS - 2006'))

# mets
plot_grid(ggplot.a.raster( dats2005.s$temp, bounds = c( 0, NA),
                           facet.names = '2005 Temperature - K'),
          ggplot.a.raster( dats2005.s$apcp, bounds = c( 0, 6),
                           facet.names = '2005 Precip - kg m^{-2}'),
          ggplot.a.raster( dats2005.s$rhum, bounds = c( 20, 80),
                           facet.names = '2005 Relative humidity - %'),
          ggplot.a.raster( dats2005.s$wspd, bounds = c( 0, 3.5),
                           facet.names = '2005 Windspeed - m/s'),
          ncol = 2, nrow = 2,
          labels = NULL, align = 'v', axis = 'lr')

# mean yhat (yhat bar) minut 2005
yhatbar2005 <- mean.lol( single.crossval, 'Y.ho.hat.bias.raster', 'y.hat.lm.cv')
ggplot.a.raster( yhatbar2005, bounds = c( -3, 3),
                 facet.names = c( 'Yhatbar - Hybrid CMAQ-DDM'))

# yhat - y, predict 2006
ggplot.a.raster( pred2006$Y.ho.hat.raster$y.hat.lm.cv - 
                   dats2006.s$cmaq.ddm, bounds = c( -2,2),
                 facet.names = c( 'Yhat - Hybrid CMAQ-DDM'))

 