

# ******************************************************************************

# Sacramento Soil Moisture Accounting (SAC-SMA) Model
# 
# The SAC-SMA is a conceptual model using a two-layer soil moisture system to
# continuously account for storage and flow through the soil layers. 
# 
# The R version of the model is developed based on the MATLAB code (Sungwook Wi)

# Latest update: 6/30/2017
# # Contact: M. Umit Taner (tanerumit@gmail.com)  

# ******************************************************************************

########## MODEL COMPONENTS

source("./R/sac_sma.R")
source("./R/sac_sma_dsm.R")
source("./R/PET functions.R")


########## SET FILE DIRECTORY LOCATIONS 

# Main file directory
parent_filedir <- "C:/Users/Umit/Google Drive/Research/Projects/_DEV/sac_sma_test/"

# Directory containing HRU files
hru_filedirectory <- paste0(parent_filedir,"HRU_climfiles_SACSMA/")

# GridInfo file location
hru_infodir <- paste0(parent_filedir,"HRUfiles/")

######### READ-IN DATA

#Grid information file
#each row represents a unique hru, columns represent hru features
gridinfo <- read.table(paste0(hru_infodir, 'HRUinfo_arrhon.txt'))

grid_lat     <- gridinfo[[1]]   # latitude 
grid_lon     <- gridinfo[[2]]   # longitude  
grid_area    <- gridinfo[[3]]   # area
grid_elev    <- gridinfo[[4]]   # elevation
grid_flowlen <- gridinfo[[5]]   # flow length 

#Calibrated parameters
calib_par <- read.table(paste0(hru_infodir,'hru_optpar_0425_arrhon_KGE75.txt'))

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#Simulate model
results <- sac_sma_dsm(str_date     = as.Date("1995/10/1"), 
                       end_date     = as.Date("2014/09/30"), 
                       grid_lat     = grid_lat, 
                       grid_lon     = grid_lon,
                       grid_area    = grid_area,
                       grid_elev    = grid_elev,
                       grid_flowlen = grid_flowlen,
                       grid_par     = calib_par)

################################################################################


#Compare model outputs
mflow <- read.table("C:\\Users\\Umit\\Desktop\\sacsmatotflow.txt") %>% unlist() %>% as.numeric()
rflow <- hru_simflow

flow_data <- data_frame(Date = sim_date, matlab = mflow, rflow = rflow)

flow_data_long <- flow_data %>%
  gather(key = variable, value = value, -Date)

ggplot(flow_data_long, aes(x = Date, y = value, color = variable)) +
  geom_line(alpha = 0.5)


View(flow_data)

