library( data.table)

## =========================================================== ##
## Functions
## =========================================================== ##
## Load adjoint runs, perform rankings, merge
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
## Read in GC adjoint
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


library( tidycensus)
census_api_key("YOUR API KEY GOES HERE")
us_pops <- get_estimates(geography = "state", product = "population", year = c( 2006, 2011),
                         state = state.name[state.abb %in% unique( adj.other_states$state)])



