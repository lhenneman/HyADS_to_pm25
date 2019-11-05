rm( list = ls())

source( '~/Dropbox/Harvard/RFMeval_Local/HyADS_to_pm25/RCode/hyads_to_pm25_functions.R')

#coordinate reference system projection string for spatial data
p4s <- "+proj=lcc +lat_1=33 +lat_2=45 +lat_0=40 +lon_0=-97 +a=6370000 +b=6370000"

#======================================================================#
## Load saved object of cmaq-ddm / hyads models from hyads_to_pm25_month.R
#======================================================================#
