

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Sacramento Soil Moisture Accounting (SAC-SMA) Model
# 
# The SAC-SMA is a conceptual model using a two-layer soil moisture system to
# continuously account for storage and flow through the soil layers. 
# 
# The R version of the model is developed based on the MATLAB code (Sungwook Wi)

# Latest update: 7/3/2017
# # Contact: M. Umit Taner (tanerumit@gmail.com)  
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


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

results <- sac_sma_dsm(str_date     = as.Date("1995/10/1"), 
                       end_date     = as.Date("2014/09/30"), 
                       grid_lat     = gridinfo[[1]], 
                       grid_lon     = gridinfo[[2]],
                       grid_area    = gridinfo[[3]],
                       grid_elev    = gridinfo[[4]],
                       grid_flowlen = gridinfo[[5]],
                       grid_par     = calib_par)



