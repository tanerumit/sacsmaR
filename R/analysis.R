

# MODEL SETTINGS ---------------------------------------------------------------

source("./R/sac_sma.R")
source("./R/sac_sma_dsm.R")
source("./R/PET functions.R")
source("./R/route_lohamann.R")

# Main file directory
main_filedir <- "C:/Users/Umit/Google Drive/Research/Projects/_DEV/sac_sma_test/"

# Directory containing HRU files
hru_filedir <- paste0(main_filedir,"HRU_climfiles_SACSMA/")

# GridInfo file location
hru_infodir <- paste0(main_filedir,"HRUfiles/")

#Grid information file, each row = a unique hru, each column  = a feature
gridinfo <- read.table(paste0(hru_infodir, 'HRUinfo_arrhon.txt'))

#Calibrated sac-sma parameters
calib_par <- read.table(paste0(hru_infodir,'hru_optpar_0425_arrhon_KGE75.txt'))

# SIMULATION ANALYSIS ----------------------------------------------------------

str_date = as.Date("1995/10/1")
end_date = as.Date("2014/09/30")
grid_lat = gridinfo[[1]]
grid_lon = gridinfo[[2]]
grid_area = gridinfo[[3]]
grid_elev = gridinfo[[4]]
grid_flowlen = gridinfo[[5]]
grid_par  = calib_par
flag_SNOW17 = 0

#profvis({
  
results <- sac_sma_dsm(str_date = as.Date("1995/10/1"), 
                         end_date = as.Date("2014/09/30"), 
                         grid_lat = gridinfo[[1]], 
                         grid_lon = gridinfo[[2]],
                         grid_area = gridinfo[[3]],
                         grid_elev = gridinfo[[4]],
                         grid_flowlen = gridinfo[[5]],
                         grid_par  = calib_par)
  
#})


write.csv(results, "./Data/r_hru_all_nosnow.csv")
write.csv(results, "./Data/r_routing_hrus.csv")
#write.csv(results, "./Data/r_totflow_nosnow.csv")
