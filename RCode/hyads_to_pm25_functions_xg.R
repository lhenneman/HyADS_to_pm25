# assuming hyads_to_pm25_functions.R is loaded as well as rnaturalearth, USABoundaries,
# and all relevant datasets
# -----------------------------------
# source('./RCode/hyads_to_pm25_functions.R')
library(xgboost)
# ---
# Define the xgboost model function
#  @params rounds -the number of trees to grow, i.e. max # of iterations used to build model
#          stopping_rounds -if after some iteration, stopping_rounds trees have been built 
#                           and no improvement to a given evaluation metric has resulted, 
#                           stop training
#          approx.contrib -boolean indicating whether to approximate covariate contributions
#                          or not when they are calculated; the exact algorithm takes time
#          hyp.params -list of parameters that govern xgboost's training algorithm, see docs
#          ... -additional args to pass to xgb.train(), see docs
# ---
xg.hyads.ddm.holdout <- function( seed.n = NULL,
                                  dat.stack = NULL,
                                  dat.stack.pred = NULL,
                                  y.name = 'cmaq.ddm',
                                  x.name = 'hyads',
                                  name.idwe = 'tot.sum',
                                  covars.names = NULL, 
                                  ho.frac = 0.1,
                                  rounds = 200,
                                  stopping_rounds = 15,
                                  approx.contrib = F,
                                  hyp.params = list(max_depth = 11,
                                                    eta = 0.01,
                                                    eval_metric = 'rmse'),
                                  return.mods = F,
                                  ...) {
  set.seed( seed.n)
  
  # if no covar names provided, use all covariates that aren't x or y
  if( is.null( covars.names))
    covars.names <- names( dat.stack)[names( dat.stack) %ni% c( x.name, y.name)]
  
  # define holdout parameters
  N <- ncell( dat.stack)
  ho <- sample( 1:N, ceiling( ho.frac * N))
  
  # extract coordinates
  dat.coords <- coordinates( dat.stack)
  dat.coords.ho <- data.table( dat.coords[ho,])
  dat.coords.tr <- data.table( dat.coords[-ho,])
  
  # define inputs
  dat.stack.ho <- data.table( values( dat.stack)[ ho,])
  dat.stack.tr <- data.table( values( dat.stack)[-ho,])
  
  # special case for ho is zero
  if( length( ho) == 0){
    # extract coordinates
    dat.coords.ho <- data.table( dat.coords)
    dat.coords.tr <- data.table( dat.coords)
    
    # define inputs
    dat.stack.ho <- data.table( dat.coords.ho, values( dat.stack))
    dat.stack.tr <- data.table( dat.coords.tr, values( dat.stack))
  } 
  if( !is.null( dat.stack.pred)){
    # extract coordinates
    dat.coords.ho <- data.table( coordinates( dat.stack.pred))
    
    # define inputs
    dat.stack.ho <- data.table( dat.coords.ho, values( dat.stack.pred))
  } 
  
  # remove NAs
  dat.stack.ho <- na.omit(dat.stack.ho)
  dat.stack.tr <- na.omit(dat.stack.tr)
  
  # remove coordinates in coordinates dts if they don't exist in dat.stacks
  # (i.e. do an inner join)
  dat.coords.ho <- merge(dat.coords.ho, dat.stack.ho, by = c('x', 'y'))
  dat.coords.tr <- merge(dat.coords.tr, dat.coords.tr, by = c('x', 'y'))
  
  # pass training data into xgboost, only takes matrices
  dat.mat.tr <- as.matrix(dat.stack.tr[,!y.name, with=F])
  dat.mat.ho <- as.matrix(dat.stack.ho[,!y.name, with=F])
  dat.mat.y.tr <- as.matrix(dat.stack.tr[, y.name, with=F])
  dat.mat.y.ho <- as.matrix(dat.stack.ho[, y.name, with=F])
  
  xtrain1c <- xgb.DMatrix(dat.mat.tr, label=dat.mat.y.tr)
  xhold1c <- xgb.DMatrix(dat.mat.ho, label=dat.mat.y.ho)
  
  # train with the given rounds and hyperparameters, supress eval output, stop training
  # if the error on the training does not improve after stopping_rounds rounds
  xgm.cv <- xgb.train(data = xtrain1c, params = hyp.params, nrounds = rounds, verbose = 0, 
                      early_stopping_rounds = stopping_rounds, watchlist=list(validation=xtrain1c),
                      ...)
  
  # make simpler models first
  mean.y.over.x <- mean( unlist( dat.stack.tr[,..y.name]) / unlist( dat.stack.tr[,..x.name]), 
                         na.rm = T)
  mean.y <- mean( unlist( dat.stack.tr[,..y.name]), na.rm = T)
  mean.x <- mean( unlist( dat.stack.tr[,..x.name]), na.rm = T)
  sd.y <- sd( unlist( dat.stack.tr[,..y.name]), na.rm = T)
  sd.x <- sd( unlist( dat.stack.tr[,..x.name]), na.rm = T)
  
  # get the simple predictions
  y.hat.mean <- unlist( dat.stack.ho[,..x.name] * mean.y.over.x)
  y.hat.Z <- unlist( (dat.stack.ho[,..x.name] - mean.x) / sd.x * sd.y + mean.y)
  x.rscale.Z <- unlist( (dat.stack.ho[,..x.name] - mean.x) / sd.x)
  y.rscale.Z <- unlist( (dat.stack.ho[,..y.name] - mean.y) / sd.y)
  
  # full xgboost model, return predictions
  y.ho <- unlist( dat.stack.ho[,..y.name])
  y.hat.xgm.cv <- predict(xgm.cv, newdata = xhold1c)

  # extract covariate contributions and bias
  Y.ho.terms <- predict(xgm.cv, newdata = xhold1c, predcontrib = T, 
                         approxcontrib = approx.contrib)
  Y.ho.terms.dt <- data.table(dat.coords.ho, Y.ho.terms)
  
  # the calculation of standard error is complicated; WIP and ommiting for now
  
  # set up evaluation data.table
  Y.ho.hat <- data.table( dat.coords.ho, y.ho, y.hat.xgm.cv, y.hat.mean, y.hat.Z)
  Y.ho.hat.bias <- data.table( dat.coords.ho, y.ho, 
                               y.hat.xgm.cv = y.hat.xgm.cv - y.ho, 
                               y.hat.mean = y.hat.mean - y.ho, 
                               y.hat.Z = y.hat.Z - y.ho)
  # Y.ho.hat.se <- data.table( dat.coords.ho, y.hat.lm.cv = y.hat.lm.cv$se.fit, y.hat.gam.cv = y.hat.gam.cv$se.fit,  
  #                           y.hat.lm.ncv = y.hat.lm.ncv$se.fit)
  rscale.Z <- data.table( dat.coords.ho, x.rscale.Z, y.rscale.Z)
  
  # rasterize output for plots
  crs.in <- crs( dat.stack)
  Y.ho.hat.raster <- projectRaster( rasterFromXYZ( Y.ho.hat, crs = crs.in), dat.stack)
  Y.ho.terms.raster <- projectRaster( rasterFromXYZ( Y.ho.terms.dt, crs = crs.in), dat.stack)
  # Y.ho.hat.se.raster <- projectRaster( rasterFromXYZ( Y.ho.hat.se, crs = crs.in), dat.stack)
  Y.ho.hat.bias.raster <- projectRaster( rasterFromXYZ( Y.ho.hat.bias, crs = crs.in), dat.stack)
  rscale.Z.raster <- projectRaster( rasterFromXYZ( rscale.Z, crs = crs.in), dat.stack)
  
  # calculate evaluation metrics
  eval.fn <- function( Yhat, Yact, mod.name){
    num.diff <- sum( Yhat - Yact, na.rm = T)
    abs.diff <- sum( abs( Yhat - Yact), na.rm = T)
    denom <- sum( Yact, na.rm = T)
    metrics <- data.table( mod.name = mod.name,
                           NMB = num.diff / denom,
                           NME = abs.diff / denom,
                           MB   = num.diff / length( Yhat),
                           RMSE = sqrt( sum( ( Yhat - Yact) ^ 2, na.rm = T) / length( Yhat)),
                           R = cor( Yhat, Yact, use = 'complete.obs'),
                           R.s = cor( Yhat, Yact, use = 'complete.obs', method = 'spearman'))
    return( metrics)
  }
  metrics.out <- rbind( eval.fn( Y.ho.hat$y.hat.xgm.cv, y.ho, 'xgm.cv'),
                        eval.fn( Y.ho.hat$y.hat.mean, y.ho, 'adj.mean'),
                        eval.fn( Y.ho.hat$y.hat.Z, y.ho, 'adj.Z'),
                        eval.fn( y.rscale.Z, x.rscale.Z, 'adj.Z.only'))
  
  # calculate correlations along quantiles of idwe
  evals.q <- data.table()
  evals.qq <- data.table()
  if( ho.frac == 0){
    vals.idwe <- unlist( data.table( values( dat.stack))[, ..name.idwe])
    vals.eval <- data.table( values( Y.ho.hat.raster))
    for ( s in seq( 0.01, 1, .01)){
      q <- quantile( vals.idwe, s, na.rm = T)
      vals.use <- vals.eval[vals.idwe < q,]
      evals <- eval.fn( vals.use$y.hat.xgm.cv, vals.use$y.ho, x.name)[, s := s]
      evals.q <- rbind( evals.q, evals)
    }
    for ( s in seq( 0.05, 1, .05)){
      q <- quantile( vals.idwe, s, na.rm = T)
      qq <- quantile( vals.idwe, s - .05, na.rm = T)
      vals.use <- vals.eval[vals.idwe < q & vals.idwe > qq,]
      evals <- eval.fn( vals.use$y.hat.xgm.cv, vals.use$y.ho, x.name)[, s := s]
      evals.qq <- rbind( evals.qq, evals)
    }
  }
  
  # listify the models
  if( return.mods)
    out <- list( metrics = metrics.out, 
                 evals.q = evals.q,
                 evals.qq = evals.qq,
                 model.xg = xgm.cv,
                 Y.ho.hat.raster = Y.ho.hat.raster, 
                 Y.ho.hat.bias.raster = Y.ho.hat.bias.raster,
                 Y.ho.terms.raster = Y.ho.terms.raster)
  if( !return.mods)
    out <- list( metrics = metrics.out, 
                 evals.q = evals.q,
                 evals.qq = evals.qq,
                 Y.ho.hat.raster = Y.ho.hat.raster, 
                 Y.ho.hat.bias.raster = Y.ho.hat.bias.raster,
                 Y.ho.terms.raster = Y.ho.terms.raster)
  
  return( out)
}

# do the above for several months
# @params
#   ... -args to pass to xg.hyads.ddm.holdout
month.trainer.xg <- function( name.m = names( mets2005.m)[1], 
                              name.p = names( mets2006.m)[1], 
                              name.x,
                              y.m, 
                              ddm.m = ddm.m.all, 
                              mets.m = mets.m.all,
                              emiss.m = d_nonegu.r,
                              idwe.m = idwe.m,
                              .mask.use = NULL,
                              cov.names = NULL,
                              ...) { 
  
  # create training dataset
  ddm.use <- ddm.m[[name.m]]
  hyads.use <- y.m[[name.m]]
  mets.use <- mets.m[[name.m]]
  emiss.use <- emiss.m
  idwe.use <- idwe.m[[name.m]]
  
  # create prediction dataset
  ddm.use.p <- ddm.m[[name.p]]
  hyads.use.p <- y.m[[name.p]]
  mets.use.p <- mets.m[[name.p]]
  idwe.use.p <- idwe.m[[name.p]]
  
  # fix names
  names( ddm.use)     <- 'cmaq.ddm'
  names( ddm.use.p)   <- 'cmaq.ddm'
  names( hyads.use)   <- name.x
  names( hyads.use.p) <- name.x
  names( idwe.use)   <- 'tot.sum.idwe'
  names( idwe.use.p) <- 'tot.sum.idwe'
  
  # combine each dataset as stacks
  dat.s <- project_and_stack( ddm.use,   hyads.use,   idwe.use,   
                              mets.use, emiss.use, mask.use = .mask.use)
  dat.p <- project_and_stack( ddm.use.p, hyads.use.p, idwe.use.p, 
                              mets.use.p, emiss.use, mask.use = .mask.use)
  
  # do the modeling
  pred <- xg.hyads.ddm.holdout( dat.stack = dat.s, dat.stack.pred = dat.p, x.name = name.x,
                                ho.frac = 0, covars.names = cov.names, name.idwe = 'tot.sum.idwe', 
                                return.mods = T, ...)
  
  return( pred)
}