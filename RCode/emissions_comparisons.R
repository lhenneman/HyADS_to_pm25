library( data.table)

d_cmaq <- fread( "~/Dropbox/Harvard/RFMeval_Local/CMAQ_DDM/COAL IMPACTS/INVENTORY/CEM/2005_cemsum.txt")[HTINPUT > 0]
d_nonegu <- fread( "~/Dropbox/Harvard/RFMeval_Local/CMAQ_DDM/COAL IMPACTS/INVENTORY/NONEGU COAL/ptinv_ptnonipm_xportfrac_cap2005v2_2005cs_orl_06jan2011_v4_orl_COAL.txt")


