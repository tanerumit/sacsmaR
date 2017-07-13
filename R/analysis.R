

# MODEL SETTINGS ---------------------------------------------------------------

source("./R/sac_sma.R")
source("./R/sac_sma_dsm.R")
source("./R/PET functions.R")
source("./R/route_lohamann.R")

# Main file directory
main_filedir <- "/Users/Umit/Google Drive/Research/Projects/_DEV/sac_sma_test/"

# Directory containing HRU files
hru_filedir <- paste0(main_filedir,"HRU_climfiles_SACSMA/")

#Grid information file, each row = a unique hru, each column  = a feature
grid_info <- read.table(paste0(main_filedir, '/HRUfiles/HRUinfo_arrhon.txt'))

#Climate information data
grid_data <- list()
num_hru <- nrow(grid_info)
for (n in 1:num_hru) {
  lat  <- formatC(grid_info[n,1], format = 'f', flag='0', digits = 5)
  lon  <- formatC(grid_info[n,2], format = 'f', flag='0', digits = 5)
  grid_data[[n]] <- read.table(paste0(hru_filedir,lat, "_", lon)) 
}

#Calibrated sac-sma parameters
calib_par <- read.table(paste0(main_filedir,'/HRUfiles/hru_optpar_0425_arrhon_KGE75.txt'))

# SIMULATION ANALYSIS ----------------------------------------------------------

str_date = as.Date("1995/10/1")
end_date = as.Date("2014/09/30")
grid_lat = grid_info[[1]]
grid_lon = grid_info[[2]]
grid_area = grid_info[[3]]
grid_elev = grid_info[[4]]
grid_flowlen = grid_info[[5]]
grid_par  = calib_par
flag_SNOW17 = 0

Rprof("profile.out")
results <- sac_sma_dsm(str_date = as.Date("1995/10/1"), 
                         end_date = as.Date("2014/09/30"), 
                         grid_lat = grid_info[[1]], 
                         grid_lon = grid_info[[2]],
                         grid_area = grid_info[[3]],
                         grid_elev = grid_info[[4]],
                         grid_flowlen = grid_info[[5]],
                         grid_par  = calib_par)
Rprof(NULL)
summaryRprof("profile.out")

write.csv(results, "./Data/r_hru_all_nosnow.csv")
write.csv(results, "./Data/r_routing_hrus.csv")
#write.csv(results, "./Data/r_totflow_nosnow.csv")
