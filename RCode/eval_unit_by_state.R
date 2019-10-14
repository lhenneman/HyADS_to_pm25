library( data.table)

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
