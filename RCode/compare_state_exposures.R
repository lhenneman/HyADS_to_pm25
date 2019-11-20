# copied from evaluate_RFMs

library( arepa)
library( hyspdisp)
library( data.table)
library( ggplot2)

rm( list = ls())

## ================================================================= ##
## Function
## ================================================================= ##
rankfacs_by_popwgt_location <- function( link.ampd.files = NULL,
                                         link.dt = NULL,
                                         cba.dt = copy( d_CBA_2007_2012),
                                         metrics = c( 'dist.km', 'sox.inverse.dist'),
                                         zip.value = '*',
                                         state.value = '*',
                                         city.value = '*',
                                         popwgt = T){
  gc()
  
  
  # read in HyADS link file and ampd file
  if( !is.null( link.ampd.files) & is.null( link.dt)){
    ziplinks_hyads <- fread(link.ampd.files[['hyads.file']])[, V1 := NULL]
    ziplinks_hyads[, `:=` ( year = as.integer( gsub( '_.*$', '', yearmonth)))]
    
    
    ziplinks_ampd  <- fread(link.ampd.files[['ampd.file']])[, V1 := NULL][uID %in% ziplinks_hyads$uID]
    
    year.dt.creator <- function( unique_yearmonths){
      out <- data.table( year = as.integer( gsub( '_.*$', '', unique_yearmonths)),
                         year_month = unique( unique_yearmonths))
      return( out)
    }
    
    ampd_years.dt <- year.dt.creator( unique( ziplinks_ampd$year_month))
    ziplinks_ampd <- merge( ziplinks_ampd, ampd_years.dt, by = 'year_month') 
    #[,  `:=` ( year = as.integer( gsub( '_.*$', '', year_month)))]
    
    link.dt <- merge( ziplinks_hyads, ziplinks_ampd, 
                      by = c( 'ZIP', 'uID', 'year'),
                      all = T)
    
    link.dt[, ZIP := formatC(ZIP, width = 5, format = "d", flag = "0")]
    print( paste( Sys.time(), 'HyADS/IDWE linked for', 
                  na.omit( unique( link.dt$year_month))))
    
  } 
  
  ## Retrieve list of all zip codeis to segnment by state
  ZIP <- get_zip_codes()[, .( zip, State.zip, City.zip)]
  setnames(ZIP, "zip", "ZIP")
  
  ## Merge ZIP code and census info with link.dt
  link.dt <- merge( link.dt, ZIP, by = 'ZIP')
  link.dt <- merge( link.dt, cba.dt, by = c( 'ZIP', 'year'))
  
  ## limit data table to subset.value in subset.feature 
  zip.search   <- paste0( zip.value,   collapse = '|')
  state.search <- paste0( state.value, collapse = '|')
  city.search  <- paste0( city.value,  collapse = '|')
  
  link.dt.trim <- link.dt[ Reduce( intersect, list( grep( zip.search, ZIP),
                                                    grep( state.search, State.zip),
                                                    grep( city.search, City.zip)))]
  
  ## Replace NA's with 0
  link.dt.trim[ is.na( link.dt.trim)] <- 0
  
  ## Weight metric by popultion
  names.py <- paste0( metrics, '.py')
  # define different functions depending on whether popwgt == F or not
  if( popwgt){
    link.dt.trim[, (names.py) := lapply(metrics, function(x) { TOTPOP_CY * get(x)})]
  } else
    link.dt.trim[, (names.py) := lapply(metrics, function(x) { get(x)})]  
  
  ## Sum pop-weighted metrics by uID
  names.py.mean <- paste0( metrics, '.py.mean')
  uID.pw <- link.dt.trim[, lapply(names.py, function(x) { mean( get( x), na.rm = T) * .N}),
                         by = c("uID", "year")]
  setnames(uID.pw, paste0( 'V', 1:length( names.py.mean)), names.py.mean)
  
  ## Rank facilities by metric in each year
  names.py.rank <- paste0( metrics, '.rank')
  uID.pw[, (names.py.rank) := lapply(names.py.mean, function(x) { frankv( get(x), order = -1)}),
         by = year]
  
  return( uID.pw)  
}

## ================================================================= ##
# Read in census data (use 2007 as 2005 and 2006 proxy)
## ================================================================= ##
d_CBA_2007_2012 <- fread( "~/shared_space/ci3_nsaph/LucasH/RFM_health/intermediates/CBA_07_12_zip_links.csv")[, .(ZIP, 
                                                                                                                  TOTPOP_CY,
                                                                                                                  year)]
d_CBA_2007_2012$ZIP  <- formatC( d_CBA_2007_2012$ZIP, 
                                 width = 5, 
                                 format = "d", 
                                 flag = "0")

d_CBA_2006 <- d_CBA_2007_2012[year == 2005][, year := 2006]
d_CBA_2007_2012 <- rbind( d_CBA_2007_2012, d_CBA_2006)
rm(d_CBA_2006)
## ================================================================= ##
## Read in HyADS unit-zipcode links and other metrics 
## for exposure, format, and merge
## ================================================================= ##
ziplinks_hyads <- fread('~/shared_space/ci3_nsaph/LucasH/RFM_health/intermediates/ZIPexposures.hyspl.disp_byunit_w2011.csv')[, V1 := NULL]
ziplinks_ampd  <- fread('~/shared_space/ci3_nsaph/LucasH/RFM_health/intermediates_source_receptor/ampd_zips_links_5km.csv')[, V1 := NULL]

ziplinks_hyads[, ZIP := formatC(ZIP, width = 5, format = "d", flag = "0")]
ziplinks_ampd[, ZIP := formatC(ZIP, width = 5, format = "d", flag = "0")]

ziplinks_hyads[, uID := gsub('_|-|\\*', '.', uID)]
ziplinks_ampd[, uID := gsub('_|-|\\*', '.', uID)]

ziplinks.dt <- merge( ziplinks_hyads, ziplinks_ampd, 
                      by = c( 'ZIP', 'uID', 'year'),
                      all = T)

## ================================================================= ##
## Rank using population weighting - states
## ================================================================= ##
rank.metrics <- c( 'sox.inverse.dist', 'nox.inverse.dist', 'hyspdisp')

ranks_KY <- rankfacs_by_popwgt_location( link.dt = ziplinks.dt,
                                         cba.dt = d_CBA_2007_2012,
                                         metrics = rank.metrics,
                                         state.value = 'KY')[, state := 'KY']
ranks_NY <- rankfacs_by_popwgt_location( link.dt = ziplinks.dt,
                                         cba.dt = d_CBA_2007_2012,
                                         metrics = rank.metrics,
                                         state.value = 'NY')[, state := 'NY']
ranks_GA <- rankfacs_by_popwgt_location( link.dt = ziplinks.dt,
                                         cba.dt = d_CBA_2007_2012,
                                         metrics = rank.metrics,
                                         state.value = 'GA')[, state := 'GA']
ranks_PA <- rankfacs_by_popwgt_location( link.dt = ziplinks.dt,
                                         cba.dt = d_CBA_2007_2012,
                                         metrics = rank.metrics,
                                         state.value = 'PA')[, state := 'PA']
ranks_WI <- rankfacs_by_popwgt_location( link.dt = ziplinks.dt,
                                         cba.dt = d_CBA_2007_2012,
                                         metrics = rank.metrics,
                                         state.value = 'WI')[, state := 'WI']
ranks_CA <- rankfacs_by_popwgt_location( link.dt = ziplinks.dt,
                                         cba.dt = d_CBA_2007_2012,
                                         metrics = rank.metrics,
                                         state.value = 'CA')[, state := 'CA']
ranks_CO <- rankfacs_by_popwgt_location( link.dt = ziplinks.dt,
                                         cba.dt = d_CBA_2007_2012,
                                         metrics = rank.metrics,
                                         state.value = 'CO')[, state := 'CO']
ranks_TX <- rankfacs_by_popwgt_location( link.dt = ziplinks.dt,
                                         cba.dt = d_CBA_2007_2012,
                                         metrics = rank.metrics,
                                         state.value = 'TX')[, state := 'TX']
ranks_US <- rankfacs_by_popwgt_location( link.dt = ziplinks.dt,
                                         cba.dt = d_CBA_2007_2012,
                                         metrics = rank.metrics)[, state := 'US']
ranks_all <- rbind( ranks_KY, ranks_NY, ranks_GA, ranks_PA, ranks_WI, ranks_US, 
                    ranks_CA, ranks_CO, ranks_TX)

## ================================================================= ##
## Rank using population weighting - cities
## ================================================================= ##
ranks_NY_syracuse <- rankfacs_by_popwgt_location( link.dt = ziplinks.dt,
                                                  cba.dt = d_CBA_2007_2012,
                                                  metrics = rank.metrics,
                                                  state.value = 'NY',
                                                  city.value = 'Syracuse')[, city := 'syracuse']
ranks_KY_lville <- rankfacs_by_popwgt_location( link.dt = ziplinks.dt,
                                                cba.dt = d_CBA_2007_2012,
                                                metrics = rank.metrics,
                                                state.value = 'KY',
                                                city.value = 'Louisville')[, city := 'louisville']
ranks_GA_atlanta <- rankfacs_by_popwgt_location( link.dt = ziplinks.dt,
                                                 cba.dt = d_CBA_2007_2012,
                                                 metrics = rank.metrics,
                                                 state.value = 'GA',
                                                 city.value = 'Atlanta')[, city := 'atlanta']
ranks_PA_pittsburgh <- rankfacs_by_popwgt_location( link.dt = ziplinks.dt,
                                                    cba.dt = d_CBA_2007_2012,
                                                    metrics = rank.metrics,
                                                    state.value = 'PA',
                                                    city.value = 'Pittsburgh')[, city := 'pittsburgh']
ranks_WI_madison <- rankfacs_by_popwgt_location( link.dt = ziplinks.dt,
                                                 cba.dt = d_CBA_2007_2012,
                                                 metrics = rank.metrics,
                                                 state.value = 'WI',
                                                 city.value = 'Madison')[, city := 'madison']
ranks_city <- rbind( ranks_KY_lville, ranks_NY_syracuse, ranks_GA_atlanta, 
                     ranks_PA_pittsburgh, ranks_WI_madison)


## ================================================================= ##
## rank at states without population-weighting
## ================================================================= ##
ranks_KY_npw <- rankfacs_by_popwgt_location( link.dt = ziplinks.dt,
                                             cba.dt = d_CBA_2007_2012,
                                             metrics = rank.metrics,
                                             state.value = 'KY',
                                             popwgt = F)[, state := 'KY']
ranks_GA_npw <- rankfacs_by_popwgt_location( link.dt = ziplinks.dt,
                                             cba.dt = d_CBA_2007_2012,
                                             metrics = rank.metrics,
                                             state.value = 'GA',
                                             popwgt = F)[, state := 'GA']
ranks_PA_npw <- rankfacs_by_popwgt_location( link.dt = ziplinks.dt,
                                             cba.dt = d_CBA_2007_2012,
                                             metrics = rank.metrics,
                                             state.value = 'PA',
                                             popwgt = F)[, state := 'PA']
ranks_NY_npw <- rankfacs_by_popwgt_location( link.dt = ziplinks.dt,
                                             cba.dt = d_CBA_2007_2012,
                                             metrics = rank.metrics,
                                             state.value = 'NY',
                                             popwgt = F)[, state := 'NY']
ranks_WI_npw <- rankfacs_by_popwgt_location( link.dt = ziplinks.dt,
                                             cba.dt = d_CBA_2007_2012,
                                             metrics = rank.metrics,
                                             state.value = 'WI',
                                             popwgt = F)[, state := 'WI']
ranks_CA_npw <- rankfacs_by_popwgt_location( link.dt = ziplinks.dt,
                                             cba.dt = d_CBA_2007_2012,
                                             metrics = rank.metrics,
                                             state.value = 'CA',
                                             popwgt = F)[, state := 'CA']
ranks_US_npw <- rankfacs_by_popwgt_location( link.dt = ziplinks.dt,
                                             cba.dt = d_CBA_2007_2012,
                                             metrics = rank.metrics,
                                             popwgt = F)[, state := 'US']
ranks_states_npw <- rbind( ranks_KY_npw, ranks_NY_npw, ranks_GA_npw, ranks_PA_npw, ranks_WI_npw, ranks_US_npw)

## ================================================================= ##
## Separate data for just a few plants to investigate change with distance
## ================================================================= ##

ziplinks.dt.sub <- ziplinks.dt[uID %in% c( '703.1BLR', '2850.1', '3497.1', '1363.4', '6193.061B')]

plot( ziplinks.dt.sub$dist.km, ziplinks.dt.sub$hyspdisp)
points( ziplinks.dt.sub$dist.km, ziplinks.dt.sub$sox.inverse.dist, col = 'red')

## ================================================================= ##
## Save data
## ================================================================= ##
outdir <- "/nfs/home/H/henneman/shared_space/ci3_nsaph/LucasH/HyADS/unit_rankings"
outdir2 <- "/mnt/user_shared/henneman/henneman_desktop/Users/luh919/Dropbox/Harvard/Meetings_and_People/Colab_Columb/data/HyADS_ranks"
outdir3 <- "/mnt/user_shared/henneman/henneman_laptop/Users/lhenneman/Dropbox/Harvard/RFMeval_Local/HyADS/unit_rankings"


write.csv( ranks_states_npw,
           file = file.path( outdir, 'ranks_states_npw.csv'))
write.csv( ranks_city ,
           file = file.path( outdir, 'ranks_city.csv'))
write.csv( ranks_all,
           file = file.path( outdir, 'ranks_all_w2011.csv'))

write.csv( ziplinks.dt.sub,
           file = file.path( outdir, 'ziplinks_hyads_ampd_sub.csv'))


write.csv( ranks_NY,
           file = file.path( outdir, 'ranks_NY.csv'))
write.csv( ranks_NY_syracuse,
           file = file.path( outdir, 'ranks_NY_syracuse.csv'))
write.csv( ranks_KY,
           file = file.path( outdir, 'ranks_KY.csv'))
write.csv( ranks_lville,
           file = file.path( outdir, 'ranks_KY_lville.csv'))
write.csv( ranks_PA,
           file = file.path( outdir, 'ranks_PA.csv'))
write.csv( ranks_PA_pittsburgh,
           file = file.path( outdir, 'ranks_PA_pittsburgh.csv'))
write.csv( ranks_GA,
           file = file.path( outdir, 'ranks_GA.csv'))
write.csv( ranks_GA_atlanta,
           file = file.path( outdir, 'ranks_GA_atlanta.csv'))
write.csv( ranks_WI,
           file = file.path( outdir, 'ranks_WI.csv'))
write.csv( ranks_WI_madison,
           file = file.path( outdir, 'ranks_WI_madison.csv'))
write.csv( ranks_US_hyads,
           file = file.path( outdir, 'ranks_US_hyads.csv'))
write.csv( ranks_US_hyads.fac,
           file = file.path( outdir, 'ranks_US_hyads.fac.csv'))


## ================================================================= ##
## run multiple times for monthly eposures
## ================================================================= ##
# define list of files
file.id <- function( yearmonth = '2005_01',
                     yearmonth.follower.hyads = '',
                     location = '~/shared_space/ci3_nsaph/LucasH/RFM_health/intermediates_source_receptor/',
                     hyads.string = 'ZIPexposures_byunit_',
                     ampd.string = 'ampd_zips_links_5km'){
  
  yearmonth_ampd <- gsub("_0", '_', yearmonth)
  
  hyads.file <- list.files( file.path( location),
                            pattern = hyads.string,
                            full.names = TRUE)
  ampd.file <- list.files( file.path( location),
                           pattern = ampd.string,
                           full.names = TRUE)
  out <- c( 'hyads.file' = hyads.file[grep( paste0( yearmonth, yearmonth.follower.hyads, '.csv'), 
                                            hyads.file)],
            'ampd.file'  = ampd.file[ grep( paste0( yearmonth_ampd, '.csv'), 
                                            ampd.file)])
  
  return( out)
  
}

month_year_names06 <- c( paste0( '2006_',
                                 formatC( 1:12, width = 2, flag = '0'),
                                 sep = ''))

month_year_names11 <- c(  paste0( '2011_',
                                  formatC( 1:12, width = 2, flag = '0'),
                                  sep = ''))

names( month_year_names06) <- month_year_names06
names( month_year_names11) <- month_year_names11

hyads.ampd.files06 <- lapply( month_year_names06,
                              file.id)
hyads.ampd.files11 <- lapply( month_year_names11,
                              file.id,
                              yearmonth.follower.hyads = 'SO2..tons.')
hyads.ampd.files <- append( hyads.ampd.files06, hyads.ampd.files11)

rank.metrics <- c( 'sox.inverse.dist', 'nox.inverse.dist', 'hyspdisp')

ranks_KY.mo <- lapply( hyads.ampd.files,
                       rankfacs_by_popwgt_location,
                       cba.dt = d_CBA_2007_2012,
                       metrics = rank.metrics,
                       state.value = 'KY')
ranks_KY.mo.dt <- rbindlist( ranks_KY.mo, idcol = "year_month")[, state := 'KY']

# ranks_NY.mo <- lapply( hyads.ampd.files,
#                        rankfacs_by_popwgt_location,
#                        cba.dt = d_CBA_2007_2012,
#                        metrics = rank.metrics,
#                        state.value = 'NY')
# ranks_NY.mo.dt <- rbindlist( ranks_NY.mo, idcol = "year_month")[, state := 'NY']

ranks_GA.mo <- lapply( hyads.ampd.files,
                       rankfacs_by_popwgt_location,
                       cba.dt = d_CBA_2007_2012,
                       metrics = rank.metrics,
                       state.value = 'GA')
ranks_GA.mo.dt <- rbindlist( ranks_GA.mo, idcol = "year_month")[, state := 'GA']

ranks_PA.mo <- lapply( hyads.ampd.files,
                       rankfacs_by_popwgt_location,
                       cba.dt = d_CBA_2007_2012,
                       metrics = rank.metrics,
                       state.value = 'PA')
ranks_PA.mo.dt <- rbindlist( ranks_PA.mo, idcol = "year_month")[, state := 'PA']
#rm(ranks_PA.mo)
gc()
ranks_WI.mo <- lapply( hyads.ampd.files,
                       rankfacs_by_popwgt_location,
                       cba.dt = d_CBA_2007_2012,
                       metrics = rank.metrics,
                       state.value = 'WI')
ranks_WI.mo.dt <- rbindlist( ranks_WI.mo, idcol = "year_month")[, state := 'WI']
#rm(ranks_WI.mo)
gc()
ranks_CA.mo <- lapply( hyads.ampd.files,
                       rankfacs_by_popwgt_location,
                       cba.dt = d_CBA_2007_2012,
                       metrics = rank.metrics,
                       state.value = 'CA')
ranks_CA.mo.dt <- rbindlist( ranks_CA.mo, idcol = "year_month")[, state := 'CA']
#rm(ranks_CA.mo)
gc()
ranks_CO.mo <- lapply( hyads.ampd.files,
                       rankfacs_by_popwgt_location,
                       cba.dt = d_CBA_2007_2012,
                       metrics = rank.metrics,
                       state.value = 'CO')
ranks_CO.mo.dt <- rbindlist( ranks_CO.mo, idcol = "year_month")[, state := 'CO']
#rm(ranks_CO.mo)
gc()
ranks_TX.mo <- lapply( hyads.ampd.files,
                       rankfacs_by_popwgt_location,
                       cba.dt = d_CBA_2007_2012,
                       metrics = rank.metrics,
                       state.value = 'TX')
ranks_TX.mo.dt <- rbindlist( ranks_TX.mo, idcol = "year_month")[, state := 'TX']
#rm(ranks_TX.mo)
gc()
ranks_US.mo <- lapply( hyads.ampd.files,
                       rankfacs_by_popwgt_location,
                       cba.dt = d_CBA_2007_2012,
                       metrics = rank.metrics)
ranks_US.mo.dt <- rbindlist( ranks_US.mo, idcol = "year_month")[, state := 'US']
#rm(ranks_US.mo)


# ranks_all.mo_KY_NY_GA <- rbind( ranks_KY.mo.dt, ranks_NY.mo.dt, ranks_GA.mo.dt )
# write.csv( ranks_all.mo_KY_NY_GA,
#            file = file.path( outdir, 'ranks_all_mo_KY_NY_GA.csv'))

ranks_all.mo <- rbind(  ranks_KY.mo.dt, #ranks_NY.mo.dt, 
                        ranks_GA.mo.dt, ranks_PA.mo.dt, 
                        ranks_WI.mo.dt, ranks_CA.mo.dt, 
                        ranks_TX.mo.dt, 
                        ranks_CO.mo.dt,
                        ranks_US.mo.dt)

outdir <- "/nfs/home/H/henneman/shared_space/ci3_nsaph/LucasH/HyADS/unit_rankings"
#ranks_all.mo_other <- fread( file.path( outdir, 'ranks_all_mo_moststates_2011.csv'))[, V1 := NULL]
#ranks_all.mo <- rbind( ranks_all.mo_other,
#                       ranks_all.mo)

write.csv( ranks_all.mo,
           file = file.path( outdir, 'ranks_all_mo_2006_2011.csv'))


## ================================================================= ##
## calculate average distance to population-weighted state centroids
## ================================================================= ##
library( hyspdisp)
ziplinks_ampd.06 <- merge( units2006, ziplinks_ampd[year == 2006], by = 'uID')
ziplinks_ampd.06[, dist := SOx / sox.inverse.dist]

## Retrieve list of all zip codeis to segnment by state
ZIP <- get_zip_codes()[, .( zip, State.zip, City.zip)]
setnames(ZIP, "zip", "ZIP")

## Merge ZIP code and census info with link.dt
link.dt <- merge( ziplinks_ampd.06, ZIP, by = 'ZIP')
link.dt <- merge( link.dt, d_CBA_2007_2012, by = c( 'ZIP', 'year'))

## Merge ZIP code and census info with link.dt
dists.us <- link.dt[, .( avedist.pew = sum( dist * SOx * TOTPOP_CY, na.rm = T) / sum( SOx * TOTPOP_CY, na.rm = T),
                         avedist.pw = sum( dist * TOTPOP_CY, na.rm = T) / sum( TOTPOP_CY, na.rm = T),
                         avedist = mean( dist, na.rm = T))][, state := 'US']
states <- c( 'PA', 'KY', 'GA', 'WI', 'CO', 'TX', 'CA')
dists.st <- link.dt[ State.zip %in% states, 
                     .( avedist.pew = sum( dist * SOx * TOTPOP_CY, na.rm = T) / sum( SOx * TOTPOP_CY, na.rm = T),
                        avedist.pw = sum( dist * TOTPOP_CY, na.rm = T) / sum( TOTPOP_CY, na.rm = T),
                        avedist = mean( dist, na.rm = T)),
                     by = State.zip]
setnames( dists.st, 'State.zip', 'state')

dists.out <- rbind( dists.us, dists.st)

## do it for hyads
ziplinks_hyads.06 <- merge( units2006, ziplinks_hyads[year == 2006], by = 'uID')
ziplinks_hyads.06[, hyspdisp.raw := hyspdisp / SOx]


## Merge ZIP code and census info with link.dt
link.dt <- merge( ziplinks_hyads.06, ZIP, by = 'ZIP')
link.dt <- merge( link.dt, d_CBA_2007_2012, by = c( 'ZIP', 'year'))

## Merge ZIP code and census info with link.dt
hyads.us <- link.dt[, .( avehyads.pew = sum( hyspdisp.raw * SOx * TOTPOP_CY, na.rm = T) / sum( SOx * TOTPOP_CY, na.rm = T),
                         avehyads.pw = sum( hyspdisp.raw * TOTPOP_CY, na.rm = T) / sum( TOTPOP_CY, na.rm = T),
                         avehyads = mean( hyspdisp.raw, na.rm = T))][, state := 'US']
states <- c( 'PA', 'KY', 'GA', 'WI', 'CO', 'TX', 'CA')
hyads.st <- link.dt[ State.zip %in% states, 
                     .( avehyads.pew = sum( hyspdisp.raw * SOx * TOTPOP_CY, na.rm = T) / sum( SOx * TOTPOP_CY, na.rm = T),
                        avehyads.pw = sum( hyspdisp.raw * TOTPOP_CY, na.rm = T) / sum( TOTPOP_CY, na.rm = T),
                        avehyads = mean( hyspdisp.raw, na.rm = T)),
                     by = State.zip]
setnames( hyads.st, 'State.zip', 'state')

hyads.out <- rbind( hyads.us, hyads.st)

dists.out <- merge( dists.out, hyads.out, by = 'state')

write.csv( dists.out,
           file = file.path( outdir, 'ave_distances.csv'))

## ================================================================= ##
## Read in data - (working on Mac)
## ================================================================= ##
## read in HyADS and emissions - only data
hyads_links <- list.files( '~/Dropbox/Harvard/RFMeval_Local/Comparisons_Intermodel/evaluate_RFMs_intermediates',
                           pattern = '^ranks_..\\.csv',
                           full.names = T)

hyads_links_dt <- rbindlist( lapply( hyads_links,
                                     function( x){
                                       x.dt <- fread( x)[, V1 := NULL]
                                       state <- gsub( '.*ranks_|.csv', '', x)
                                       x.dt[, state := state]
                                       return (x.dt)
                                     }))

## read in adjoint data
adjoint_links <- list.files( '~/Dropbox/Harvard/RFMeval_Local/Adjoint/other_states',
                             full.names = T)
col.names <- c( 'ID',	'Latitude',	'Longitude',	'SOx',	'CO2',	'NOx',
                'Height',	'Diam',	'Velocity',	'Temp',	
                'NOx exposure',	'SOx exposure',	'Total')
adjoint_links_dt <- rbindlist( lapply( adjoint_links,
                                       function( x){
                                         x.dt <- fread( x)[, V1 := NULL]
                                         names( x.dt) <- col.names
                                         x.dt[, uID := gsub('_|-|\\*', '.', ID)]
                                         x.dt[, ID := NULL]
                                         
                                         state <- gsub( '.*adjoint_|.csv', '', x)
                                         state <- ifelse( state == 'USA', 'US', state)
                                         
                                         year <- as.integer( gsub( '.*_units_|_adjoint_.*', '', x))
                                         
                                         x.dt[, `:=` (state = state,
                                                      year = year)]
                                         return (x.dt)
                                       }))

adjoint_links_dt[, `:=` (uID.rank.noxadj = frankv( `NOx exposure`, order = -1),
                         uID.rank.soxadj = frankv( `SOx exposure`, order = -1),
                         uID.rank.totadj = frankv( `Total`, order = -1)), 
                 by = c('state', 'year')]



adj_hysp_ranks <- merge( hyads_links_dt, adjoint_links_dt, by = c( 'uID', 'state', 'year'))
## ================================================================= ##
## Read in data - old version, will update
## ================================================================= ##
hyads <- fread( '~/Dropbox/Harvard/RFMeval_Local/Comparisons_Intermodel/evaluate_RFMs_intermediates/ZIPexposures.hyspl.disp__state_byunit.csv')
hyadsKY <- hyads[ State.zip == 'KY']
hyadsKY[, `:=` (uID = gsub('_|-|\\*', '.', uID),
                V1 = NULL)]

adjoint05 <- fread( '~/Dropbox/Harvard/RFMeval_Local/Adjoint/final_merge_nei_ampd_all_units_2005_adjoint.csv')
adjoint05[, `:=` (year = 2005)]
adjoint12 <- fread( '~/Dropbox/Harvard/RFMeval_Local/Adjoint/final_merge_nei_ampd_all_units_2012_adjoint.csv')
adjoint12[, `:=` (year = 2012)]

adjoint <- rbind( adjoint05,
                  adjoint12)
adjoint[, uID := gsub( ' |"', '', ID)]
adjoint[, uID := gsub('_|-|\\*', '.', uID)]


datKY <- na.omit( merge( hyadsKY, 
                         adjoint,
                         by = c( 'uID', 'year'),
                         all = T))

## ================================================================= ##
## Rank exposures
## ================================================================= ##
datKY[, `:=` (uID.rank.hysp = frankv( hyspdisp.PW, order = -1),
              uID.rank.adj  = frankv( `SOx exposure`, order = -1))]

plot(datKY[year == 2005, uID.rank.hysp], datKY[year == 2005, uID.rank.adj])
plot(datKY[year == 2012, uID.rank.hysp], datKY[year == 2012, uID.rank.adj])

plot(datKY[year == 2005, hyspdisp.PW], datKY[year == 2005, `SOx exposure`])
plot(datKY[year == 2012, hyspdisp.PW], datKY[year == 2012, `SOx exposure`])

cor(datKY[year == 2005, uID.rank.hysp], datKY[year == 2005, uID.rank.adj])
cor(datKY[year == 2012, uID.rank.hysp], datKY[year == 2012, uID.rank.adj])

cor(datKY$hyspdisp, datKY$`SOx exposure`, use = 'complete.obs')

ggplot( data = datKY,
        aes( x = uID.rank.adj,
             y = uID.rank.hysp)) +
  geom_point() +
  facet_wrap( . ~ year, ncol = 2) +
  ylab( "HyADS unit ranks") +
  xlab( "Adjoint unit ranks") +
  theme_bw()

ggplot( data = datKY[SOx > 0],
        aes( x = hyspdisp.PW,
             y = `SOx exposure`)) +
  geom_point() +
  scale_x_continuous( trans = 'log10') +
  scale_y_continuous( trans = 'log10') +
  facet_wrap( . ~ year, ncol = 2) +
  ylab( "HyADS population-weighted exposure") +
  xlab( "Adjoint SO2 population-weighted exposure") +
  theme_bw()














