
#======================================================================#
## import data
#======================================================================#
ddm2005 <- ddm_to_zip( ddm_coal_file = '~//Dropbox/Harvard/RFMeval_Local/CMAQ_DDM/COAL_impacts_2005_update.csv',
                       Year = 2005)
ddm2006 <- ddm_to_zip( ddm_coal_file = '~//Dropbox/Harvard/RFMeval_Local/CMAQ_DDM/COAL_impacts_2006_update.csv',
                       Year = 2006)
names( ddm2005) <- 'cmaq.ddm'
names( ddm2006) <- 'cmaq.ddm'

p4s <- "+proj=lcc +lat_1=33 +lat_2=45 +lat_0=40 +lon_0=-97 +a=6370000 +b=6370000"
hyads2005.dt <- fread( '~/Dropbox/Harvard/RFMeval_Local/HyADS_to_pm25/HyADS_grid/2005grid/HyADS_grid_2005.csv', drop = 'V1')
hyads2006.dt <- fread( '~/Dropbox/Harvard/RFMeval_Local/HyADS_to_pm25/HyADS_grid/2006grid/HyADS_grid_2006.csv', drop = 'V1')
hyads2005 <- rasterFromXYZ( hyads2005.dt[, .( x, y, hyads)], crs = p4s)
hyads2006 <- rasterFromXYZ( hyads2006.dt[, .( x, y, hyads)], crs = p4s)


#======================================================================#
## download met data
#======================================================================#
#define the layer names, do the actual downloading
Sys.setenv(TZ='UTC')
layer.names <- c( "air.2m.mon.mean.nc",
                  "apcp.mon.mean.nc",
                  "rhum.2m.mon.mean.nc",
                  "vwnd.10m.mon.mean.nc",
                  "uwnd.10m.mon.mean.nc")
names( layer.names) <- c( "temp", "apcp", "rhum", "vwnd", "uwnd")

# do the data downloading
list.met <- lapply( layer.names,
                    downloader.fn,
                    dataset = 'NARR')

# take over US
mets2005 <- usa.functioner( 2005, list.met, dataset = 'NARR', return.usa.sub = F)
mets2006 <- usa.functioner( 2006, list.met, dataset = 'NARR', return.usa.sub = F)

#======================================================================#
## reproject, align on same grid
#======================================================================#
hyads2005.p <- projectRaster( hyads2005, ddm2005)
mets2005.p <- projectRaster( mets2005, ddm2005)
hyads2006.p <- projectRaster( hyads2006, ddm2006)
mets2006.p <- projectRaster( mets2006, ddm2006)

# fill NA's with linear interpolation across days
hyads2005.pa <- hyads2005.p #focal( hyads2005.p, w=matrix(1,1,1), fun = mean, na.rm = T, NAonly = T)
hyads2006.pa <- hyads2006.p #focal( hyads2005.p, w=matrix(1,1,1), fun = mean, na.rm = T, NAonly = T)
names( hyads2005.pa) <- names( hyads2005.p)
names( hyads2006.pa) <- names( hyads2006.p)

#======================================================================#
## stack up!
#======================================================================#
dats2005.s <- stack( ddm2005, hyads2005.pa, mets2005.p)
dats2006.s <- stack( ddm2006, hyads2006.pa, mets2006.p)

# crop to US
usa <- rnaturalearth::ne_countries(scale = 110, type = "countries", country = "United States of America", 
                                   geounit = NULL, sovereignty = NULL,
                                   returnclass = c("sp"))
usa.sub <- disaggregate(usa)[6,]
usa.sub <- spTransform(usa.sub, CRSobj = proj4string( dats2005.s))

dats2005.s <- crop( mask(dats2005.s, usa.sub), usa.sub)
dats2006.s <- crop( mask(dats2006.s, usa.sub), usa.sub)

## ========================================================= ##
##                Read in emissions data
## ========================================================= ##
d_cems_cmaq.f <- "~/Dropbox/Harvard/RFMeval_Local/CMAQ_DDM/COAL IMPACTS/INVENTORY/CEM/2005_cemsum.txt"
d_nonegu.f <- "~/Dropbox/Harvard/RFMeval_Local/CMAQ_DDM/COAL IMPACTS/INVENTORY/NONEGU COAL/ptinv_ptnonipm_xportfrac_cap2005v2_2005cs_orl_06jan2011_v4_orl_COAL.txt"

d_cmaq <- fread( d_cems_cmaq.f)
d_nonegu <- fread( d_nonegu.f, skip = "06029", header = F)[,1:63]
d_nonegu.names <- unlist( fread( d_nonegu.f, skip = 'FIPS,PLANTID,', header = F, nrows = 1))
names( d_nonegu) <- d_nonegu.names

d_nonegu.slim <- d_nonegu[ POLCODE == 'SO2', .( XLOC, YLOC, ANN_EMIS)]


## Convert to spatial object, take over CMAQ raster
d_nonegu.sp <- SpatialPointsDataFrame( d_nonegu.slim[, .( XLOC, YLOC)], 
                                       data.frame( d_nonegu.slim[, ANN_EMIS]),
                                       proj4string = CRS( "+proj=longlat +datum=WGS84 +no_defs"))
d_nonegu.sp <- spTransform( d_nonegu.sp, CRS( p4s))
d_nonegu.r <- rasterize( d_nonegu.sp, dats2005.s)$d_nonegu.slim...ANN_EMIS.
plot( d_nonegu.r)
plot( usa.sub, add = T)
d_nonegu.r[is.na(d_nonegu.r[])] <- 0

dats2005.s <- stack( dats2005.s, d_nonegu.r)
dats2006.s <- stack( dats2006.s, d_nonegu.r)
#======================================================================#
## pick out section for analysis
#======================================================================#
square <- extent( 0e6, 1e6, -5e5, 5e5)
square.sm <- extent( 0.25e6, 0.75e6, -2.5e5, 2.5e5)
square.lg <- extent( -.5e6, 1.5e6, -10e5, 10e5)

#======================================================================#
## crop to test area
#======================================================================#
dats2005.s.square    <- crop( dats2005.s, square)
dats2005.s.small     <- crop( dats2005.s, square.sm)
dats2005.s.large     <- crop( dats2005.s, square.lg)

#======================================================================#
## run linear model holdouts
#======================================================================#
# single.crossval.sm <- lapply( 1:50,
#                               lm.hyads.ddm.holdout,
#                               dat.stack = dats2005.s.small)
# single.crossval.la <- lapply( 1:50,
#                               lm.hyads.ddm.holdout,
#                               dat.stack = dats2005.s.large)
cov.names <- c( "temp", "rhum", "vwnd", "uwnd", "wspd", "d_nonegu.slim...ANN_EMIS.")
single.crossval <- lapply( 1:500,
                           lm.hyads.ddm.holdout,
                           dat.stack = dats2005.s,
                           covars.names = cov.names)

pred2006 <- lm.hyads.ddm.holdout( dat.stack = dats2005.s, dat.stack.pred = dats2006.s,
                                  ho.frac = 0, covars.names = cov.names, return.mods = T)
pred2005 <- lm.hyads.ddm.holdout( dat.stack = dats2006.s,dat.stack.pred = dats2005.s,
                                  ho.frac = 0, covars.names = cov.names, return.mods = T)

# save.image( '~/Dropbox/Harvard/RFMeval_Local/HyADS_to_PM25/RData/cmaqddm_to_pm25_crossval.RData')
ggplot.a.raster( pred2006$Y.ho.hat.raster$y.hat.lm.cv - dats2006.s$cmaq.ddm, 
                 pred2006$Y.ho.hat.raster$y.hat.lm.ncv - dats2006.s$cmaq.ddm, 
                 bounds = c( -2,2),
                 facet.names = c( 'Yhat - Hybrid CMAQ-DDM w/ covars', 'Yhat - Hybrid CMAQ-DDM w/o covars'))

ggplot.a.raster( pred2005$Y.ho.hat.raster$y.hat.lm.cv - dats2005.s$cmaq.ddm, 
                 pred2006$Y.ho.hat.raster$y.hat.lm.cv - dats2006.s$cmaq.ddm, 
                 pred2006$Y.ho.hat.raster$y.hat.lm.cv - pred2005$Y.ho.hat.raster$y.hat.lm.cv, 
                 dats2006.s$cmaq.ddm - dats2005.s$cmaq.ddm, 
                 bounds = c( -2,2), nrow. = 2,
                 facet.names = c( '2005: Yhat - Hybrid CMAQ-DDM', '2006: Yhat - Hybrid CMAQ-DDM', 
                                  'Difference, predicted', 'Difference, actual'))

plot( pred2005$Y.ho.hat.raster$y.hat.lm.cv - dats2005.s$cmaq.ddm,
      pred2006$Y.ho.hat.raster$y.hat.lm.cv - dats2006.s$cmaq.ddm,
      xlab = '2006 bias from 2005 model', ylab = '2005 bias from 2006 model')
summary.preds <- lm( values( pred2006$Y.ho.hat.raster$y.hat.lm.cv - dats2006.s$cmaq.ddm) ~ 
                       values( pred2005$Y.ho.hat.raster$y.hat.lm.cv - dats2005.s$cmaq.ddm))
summary( summary.preds)
abline(summary.preds)

mean.lol( single.crossval, 'Y.ho.hat.bias.raster', 'y.hat.lm.cv')
mean.lol( single.crossval, 'Y.ho.hat.bias.raster', 'y.hat.lm.ncv')


single.crossval2 <- lm.hyads.ddm.holdout( dat.stack = dats2005.s,
                                          ho.frac = 0)
plot( single.crossval2$Y.ho.hat.bias.raster$y.hat.lm.cv)
plot( single.crossval2$Y.ho.hat.bias.raster$y.hat.Z)
plot( single.crossval2$Y.ho.hat.bias.raster$y.hat.mean)

cuts=0:6
par(mfrow=c(1,2))
plot( single.crossval2$Y.ho.hat.raster$y.hat.lm.cv, main = 'Y.hat.2006')
plot( dats2005.s$cmaq.ddm, main = 'Y.2006')
par(mfrow=c(1,1))


# to do - 
# 1. limit model to over land
# 2. Build distributions of metrics for different models

# to save - input data and cross validation


#======================================================================#
## extract distributions of metrics
#======================================================================#

metrics.out <- rbindlist( lapply( single.crossval, '[[', 'metrics'))
metrics.out.m <- melt( metrics.out, id.vars = 'mod.name')
metrics.means <- metrics.out.m[, .( med = median( value)), by = .( mod.name, variable)]

ggplot( data = metrics.out.m[ mod.name %ni% 'adj.mean'],
        aes( color = mod.name)) +
  facet_wrap( . ~ variable, ncol = 1, scales = 'free') +
  geom_vline( data = metrics.means[ mod.name %ni% 'adj.mean'], 
              aes( xintercept = med, color = mod.name),
              lty = 2) +
  geom_density( aes( x = value)) + 
  scale_color_discrete( name = 'Model Name') 
#scale_x_continuous( limits = c( .23, .37))

# compare 2005 and 2006

ggplot.a.raster( pred2006$Y.ho.hat.raster$y.hat.lm.cv,
                 dats2006.s$cmaq.ddm, bounds = c( 0,5),
                 facet.names = c( 'Yhat', 'Hybrid CMAQ-DDM'))


#======================================================================#
#======================================================================#
## old analysis - compare spLM and lm
#======================================================================#
#======================================================================#
us_buonds <- sf::as_Spatial( USAboundaries::us_states( ))
us_buonds <- spTransform( us_buonds, CRS( p4s))
par(mfrow=c(1,2))
plot( ddm2005, main = 'DDM, 2005')
plot( square, add = T)
plot( square.sm, add = T)
plot( square.lg, add = T)
plot( us_buonds, add = T, border = 'grey50')
plot( hyads2005.p, main = 'HyADS, 2005')
plot( square, add = T)
plot( square.sm, add = T)
plot( square.lg, add = T)
plot( us_buonds, add = T, border = 'grey50')
par(mfrow=c(1,1))

#======================================================================#
## crop to test area
#======================================================================#
hyads2005.pc    <- crop( hyads2005.p, square)
hyads2005.pc.sm <- crop( hyads2005.p, square.sm)
hyads2005.pc.lg <- crop( hyads2005.p, square.lg)

ddm2005.pc      <- crop( ddm2005, square)
ddm2005.pc.sm   <- crop( ddm2005, square.sm)
ddm2005.pc.lg   <- crop( ddm2005, square.lg)

mets2005.pc     <- crop( mets2005.p, square)
mets2005.pc.sm  <- crop( mets2005.p, square.sm)
mets2005.pc.lg  <- crop( mets2005.p, square.lg)

#======================================================================#
## extract coordinates, values
#======================================================================#
hyads2005.dt    <- data.table( rasterToPoints( hyads2005.pc))
hyads2005.sm.dt <- data.table( rasterToPoints( hyads2005.pc.sm))
hyads2005.lg.dt <- data.table( rasterToPoints( hyads2005.pc.lg))

ddm2005.dt      <- data.table( rasterToPoints( ddm2005.pc))
ddm2005.sm.dt   <- data.table( rasterToPoints( ddm2005.pc.sm))
ddm2005.lg.dt   <- data.table( rasterToPoints( ddm2005.pc.lg))

mets2005.dt     <- data.table( rasterToPoints( mets2005.pc))
mets2005.sm.dt  <- data.table( rasterToPoints( mets2005.pc.sm))
mets2005.lg.dt  <- data.table( rasterToPoints( mets2005.pc.lg))

#======================================================================#
## set up Y & X covariates
#======================================================================#
set.inputs <- function( hyads.dt, ddm.dt, mets.dt = NULL){
  list( Y = ddm.dt$X2005,
        X = data.table( intercept = 1, 
                        hyads = hyads.dt$X2005, 
                        mets.dt),
        coords = as.matrix( hyads.dt[,.( x, y)]))
}

in.cv <- set.inputs( hyads2005.dt, ddm2005.dt, mets2005.dt)
in.sm.cv <- set.inputs( hyads2005.sm.dt, ddm2005.sm.dt, mets2005.sm.dt)
in.lg.cv <- set.inputs( hyads2005.lg.dt, ddm2005.lg.dt, mets2005.lg.dt)
in.no <- set.inputs( hyads2005.dt, ddm2005.dt)


###################################################
### benchmark time tests for different sizes
###################################################
# time goes about with the cube of ncells
nin <- 1
# normal size square
# with covariates
results_sq.cv <- splM.hyads.ddm( coords = in.cv$coords, Y = in.cv$Y, X = in.cv$X,
                                 seed.n = 1)
# without covariates
results.sq.on <- splM.hyads.ddm( coords = in.no$coords, Y = in.no$Y, X = in.no$X,
                                 seed.n = 1)

# with covariates & knots
results_sq.cv.5k <- splM.hyads.ddm( coords = in.cv$coords, Y = in.cv$Y, X = in.cv$X,
                                    knots = c( 5, 5, 0), seed.n = 1)
results_sq.cv.10k <- splM.hyads.ddm( coords = in.cv$coords, Y = in.cv$Y, X = in.cv$X,
                                     knots = c( 10, 10, 0), seed.n = 1)

# with small, large & knots
results_sm.cv <- splM.hyads.ddm( coords = in.sm.cv$coords, Y = in.sm.cv$Y, X = in.sm.cv$X,
                                 seed.n = 1)
results_sm.cv.5k <- splM.hyads.ddm( coords = in.sm.cv$coords, Y = in.sm.cv$Y, X = in.sm.cv$X,
                                    knots = c( 5, 5, 0), seed.n = 1)
results_lg.cv.5k <- splM.hyads.ddm( coords = in.lg.cv$coords, Y = in.lg.cv$Y, X = in.lg.cv$X,
                                    knots = c( 5, 5, 0), seed.n = 1)
results_lg.cv.10k <- splM.hyads.ddm( coords = in.lg.cv$coords, Y = in.lg.cv$Y, X = in.lg.cv$X,
                                     knots = c( 10, 10, 0), seed.n = 1)
results_lg.cv.20k <- splM.hyads.ddm( coords = in.lg.cv$coords, Y = in.lg.cv$Y, X = in.lg.cv$X,
                                     knots = c( 20, 20, 0), seed.n = 1)

###################################################
### put together benchmarks
###################################################
results.list <- list( results_sq.cv, results.sq.on, results_sq.cv.5k, results_sq.cv.10k,
                      results_sm.cv, results_sm.cv.5k, results_lg.cv.5k, results_lg.cv.10k, results_lg.cv.20k)

results.l.dt <- rbindlist( lapply( seq_along(results.list), 
                                   function( n, l, cv, kn){
                                     lm.n <- copy( l[[n]]$metrics)
                                     lt.n <- copy( l[[n]]$runtimes)
                                     lm.n[, `:=` (ncells.tr = lt.n$ncells.tr,
                                                  cv = cv[n])]
                                     lm.n[ mod.name == 'spLM', `:=` ( kn = kn[n],
                                                                      time = sum( lt.n$train,
                                                                                  lt.n$recover,
                                                                                  lt.n$est.y.hat))]
                                     # lm.n[ mod.name == 'lm', `:=` ( cv = cv[n])]
                                     return( lm.n)                          
                                   }, results.list, 
                                   cv = c( 'yes', 'no', rep( 'yes', 7)),
                                   kn = c( 0, 0, 5, 10, 0, 5, 5, 10, 20)))

ggplot( results.l.dt) + 
  geom_point( aes( color = kn, y = time, x = ncells.tr))
ggplot( results.l.dt) + 
  geom_point( aes( color = ncells.tr, y = RMSE, x = kn))
ggplot( results.l.dt) + 
  geom_point( aes( color = mod.name, y = NME, x = ncells.tr), alpha = 0.5)
ggplot( results.l.dt) + 
  geom_point( aes( color = mod.name, y = R, x = ncells.tr), alpha = 0.5)
ggplot( results.l.dt) + 
  geom_point( aes( color = mod.name, y = MB, x = ncells.tr), alpha = 0.5)
ggplot( results.l.dt) + 
  geom_point( aes( color = cv, y = MB, x = mod.name), alpha = 0.5)

###################################################
### old tests
###################################################
set.seed( 1)
results.sq1a <- splM.hyads.ddm( dummy.n,
                                knots = c( 5, 5),
                                coords = as.matrix( hyads2005.dt[,.( x, y)]),
                                Y = ddm2005.dt$X2005,
                                X = cbind( 1, hyads2005.dt$X2005, mets2005.dt$temp, mets2005.dt$apcp))
set.seed( 1)
results.sq2 <- splM.hyads.ddm( dummy.n,
                               knots = c( 20, 20),
                               coords = as.matrix( hyads2005.dt[,.( x, y)]),
                               Y = ddm2005.dt$X2005,
                               X = cbind( 1, hyads2005.dt$X2005, mets2005.dt$temp, mets2005.dt$apcp))
# benchmark different sizes
results.sm <- splM.hyads.ddm( dummy.n,
                              coords = as.matrix( hyads2005.sm.dt[,.( x, y)]),
                              Y = ddm2005.sm.dt$X2005,
                              X = cbind( 1, hyads2005.sm.dt$X2005))
# knots 10 => 15 went from time 3 => 15
results.lg <- splM.hyads.ddm( dummy.n,
                              knots = c( 20, 20, 0),
                              coords = as.matrix( hyads2005.lg.dt[,.( x, y)]),
                              Y = ddm2005.lg.dt$X2005,
                              X = cbind( 1, hyads2005.lg.dt$X2005, mets2005.lg.dt$temp, mets2005.lg.dt$apcp))

# set the seed
set.seed( 12345)
holdouts.eval <- rbindlist( lapply( 1, splM.hyads.ddm))
# thoughts
# 1. log first?
# 2. what to do with the uncertainties in y.hat?
# 3. what about uncertainties in HyADS and CMAQ-DDM?
# 
# To-dp
# 1. benchmark times for different sizes
# 2. can it be run with multiple X variables?

###################################################
### code chunk number 3: spRecover1
###################################################
m.i.rec <- spRecover(m.i, start=burn.in, thin=5, n.report=100)
plot( m.i.rec$p.beta.recover.samples)
beta.hat <- round(summary(m.i.rec$p.beta.recover.samples)$quantiles[c(3,1,5)],6)

###################################################
### code chunk number 4: wHatSurf1
###################################################
# The w(si)’s provide local adjustment (with structured dependence) to the mean
# and capturing the effect of unmeasured or unobserved regressors with spatial
# pattern.
w.hat <- apply(m.i.rec$p.w.recover.samples, 1, median)

w.hat.raster <- rasterFromXYZ( cbind( coords, w.hat))
plot( w.hat.raster)

theta.raster <- rasterFromXYZ( cbind( coords, m.i.rec$p.theta.recover.samples))

m.i.pred <- spPredict(m.i, start=burn.in, thin=2, pred.covars=X.ho, 
                      pred.coords=coords, verbose=FALSE)
y.hat <- apply(m.i.pred$p.y.predictive.samples, 1, quants)
par(mar=c(5,5,5,5))
plot(Y, y.hat[1,], pch=19, cex=0.5, xlab="Observed y", ylab="Predicted y", 
     ylim=range(y.hat), xlim=range(y.hat), cex.lab=2, cex.axis=2)
arrows(Y, y.hat[1,], Y, y.hat[2,], angle=90, length=0.05)
arrows(Y, y.hat[1,], Y, y.hat[3,], angle=90, length=0.05)
lines(-20:20,-20:20, col="blue")

###################################################
### code chunk number 5: runTime1
###################################################
run.time <- round(m.i$run.time[3]/60,3)

###################################################
### code chunk number 7: ppModI
###################################################
m.i <- spLM( Y ~ X - 1, coords=coords, knots=c(5, 5, 0), starting=starting,
             tuning=tuning, priors=priors, cov.model="exponential",
             modified.pp=FALSE, n.samples=n.samples, n.report=2500)


###################################################
### code chunk number 8: ppModII-IV
###################################################
m.ii <- spLM(Y~X-1, coords=coords, knots=c(5, 5, 0), starting=starting,
             tuning=tuning, priors=priors, cov.model="exponential",
             modified.pp=TRUE, n.samples=n.samples, verbose=FALSE)

m.iii <- spLM(Y~X-1, coords=coords, knots=c(10, 10, 0), starting=starting,
              tuning=tuning, priors=priors, cov.model="exponential",
              modified.pp=FALSE, n.samples=n.samples, verbose=FALSE)

m.iv <- spLM(Y~X-1, coords=coords, knots=c(10, 10, 0), starting=starting,
             tuning=tuning, priors=priors, cov.model="exponential",
             modified.pp=TRUE, n.samples=n.samples, verbose=FALSE)

quants <- function(x){
  quantile(x, prob=c(0.5, 0.025, 0.975))
}

ci.print <- function(x, digits=2){
  apply(round(apply(x, 2, quants), digits), 2, function(x){paste(x[1]," (",x[2],", ",x[3],")",  sep="")})
}

m.i <- spRecover(m.i, start=burn.in, thin=5, verbose=FALSE)
m.ii <- spRecover(m.ii, start=burn.in, thin=5, verbose=FALSE)
m.iii <- spRecover(m.iii, start=burn.in, thin=5, verbose=FALSE)
m.iv <- spRecover(m.iv, start=burn.in, thin=5, verbose=FALSE)

sub <- burn.in:n.samples

B <- as.matrix(c(1,5))
sigma.sq <- 2
tau.sq <- 1
phi <- 6

ests <- cbind(c(B[1], B[2], sigma.sq, tau.sq, phi), 
              c(ci.print(m.i$p.beta.samples[sub,]), ci.print(m.i$p.theta.samples[sub,])),
              c(ci.print(m.ii$p.beta.samples[sub,]), ci.print(m.ii$p.theta.samples[sub,])),
              c(ci.print(m.iii$p.beta.samples[sub,]), ci.print(m.iii$p.theta.samples[sub,])),
              c(ci.print(m.iv$p.beta.samples[sub,]), ci.print(m.iv$p.theta.samples[sub,])))

run.times <- c("", round(c(m.i$run.time[3]/60, m.ii$run.time[3]/60, m.iii$run.time[3]/60,  m.iv$run.time[3]/60), 2))

rel.run.times <- c("", round(c(m.i$run.time[3]/m.iv$run.time[3], m.ii$run.time[3]/m.iv$run.time[3], m.iii$run.time[3]/m.iv$run.time[3],  m.iv$run.time[3]/m.iv$run.time[3]), 2))

tab <- rbind(ests, run.times, rel.run.times)

rownames(tab) <- c("$\\beta_0$","$\\beta_1$","$\\sigma^2$", "$\\tau^2$", "$\\phi$", "Time", "Rel. time")

colnames(tab) <- c("true", "i", "ii", "iii", "iv")


###################################################
### code chunk number 9: tab1
###################################################
print(xtable(tab, 
             caption="Candidate predictive process models' parameter estimates, run-time (wall time) in minutes, and run-time relative to model $iv$. Parameter posterior summary 50 (2.5, 97.5) percentiles.", label="tab1", align="cccccc"),
      table.placement="!ht",  caption.placement="bottom", sanitize.text.function=function(x){x})



crs.usa <- crs( "+proj=lcc +lat_1=33 +lat_2=45 +lat_0=40 +lon_0=-97 +a=6370000 +b=6370000")

# download USA polygon from rnaturalearth
usa <- rnaturalearth::ne_countries(scale = 110, type = "countries", country = "United States of America", 
                                   geounit = NULL, sovereignty = NULL,
                                   returnclass = c("sp"))
usa.sub <- disaggregate(usa)[6,]
usa.sub.t <- spTransform(usa.sub, CRSobj = crs.usa)
ext <- extent(usa.sub)

r <- raster( resolution = 1, ext = extent(usa), crs = crs( usa))
crop( r, usa.sub)
crop( r, usa.sub.t)
ext.proj <- projectExtent( usa.sub, p4s)
ext.proj.t <- projectExtent( usa.sub.t, crs( usa))
crop( r, ext.proj.t)

projectExtent( raster( resolution = 1, ext, crs = crs( usa)), crs = crs( p4s))

r <- raster(xmn=-110, xmx=-90, ymn=40, ymx=60, ncols=40, nrows=40)

dat <- data.table( y = 1:10, x1 = runif( 10), x2 = cos( 1:10))
covars.names <- c( 'x1', 'x2')
y.name = 'y'
# check out linear regression models - define them
form.cv <-  as.formula( paste( y.name, '~ (', paste( c( x.name, covars.names), collapse = '+'), ') ^2'))
form.ncv <- as.formula( paste( y.name, '~', x.name))
lm.cv <-  lm( form.cv,  data = dat)
lm.ncv <- lm( form.ncv, data = dat.stack.tr)


