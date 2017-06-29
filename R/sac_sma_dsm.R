

################################################################################
# Daily Stramflow Simulation using the SAC-SMA version of modeling Distributed System

################################################################################

library(dplyr)
library(readr)
library(ggplot2)

################################################################################

#Main directory storing data files
parent_filedir <- "C:/Users/Umit/Google Drive/Research/Projects/_DEV/sac_sma_test/"

# HRU data file directory. 
hru_filedirectory <- paste0(parent_filedir,"HRU_climfiles_SACSMA/")

gridinfo_dummy <- c("37.45833","-121.75000",5.73380,570.07081,0.00000)
gridinfo <- gridinfo_dummy

hru_lat <- gridinfo[1]
hru_lon <- gridinfo[2]
hru_area <- gridinfo[3] %>% as.numeric()
hru_elev <- gridinfo[4] %>% as.numeric()
hru_flowlen <- gridinfo[5] %>% as.numeric()

str_date <- as.Date("1995/10/1")
end_date <- as.Date("2014/09/30")
sim_date <- seq.Date(str_date, end_date, by = "day")

  
#Calibrated parameters
calib_par <- read.table(paste0(parent_filedir,'HRUfiles/hru_optpar_0425_arrhon_KGE75.txt'))


# Determine whether snow module is included: 0 (not activate), 1(activate)
flag_SNOW17  <- 0

# Initial Storage States in SAC_SMA
uztwc <- 0  # Upper zone tension water storage
uzfwc <- 0  # Upper zone free water storage
lztwc <- 5  # Lower zone tension water storage
lzfsc <- 5  # Lower zone supplementary free water storage
lzfpc <- 5  # Upper zone primary free water storage
adimc <- 0  # Additional impervious area storage

# Initial states for SNOW17
W_i     <- 0  # Accumulated water equivalent of the ice portion of the snow cover (mm)
ATI     <- 0  # Antecedent Temperature Index, deg C
W_q     <- 0  # Liquid water held by the snow (mm)
Deficit <- 0  # Heat Deficit, also known as NEGHS, Negative Heat Storage

#parameters at initial states
inistates <- c(uztwc, uzfwc, lztwc, lzfsc, lzfpc, adimc, W_i, ATI, W_q, Deficit)

# Base Time of Unit Hydrographs of HRU and River in Routing model. 
# (Necessary to be changed only for a VERY LARGE watershed)
KE      <-  12  # Base time for HRU UH (day)
UH_DAY  <-  96  # Base time for river routing UH (day)

# Run Distributed System Modeling version of SAC_SMA for each HRU
num_hru  <- 1 #length(hru_lat);                    # Total number of HRUs

tot_area <- sum(hru_area)                         # Total watershed area
sim_per  <- interval(str_date, end_date) / days(1) + 1 # Simulation period

#Run through each HRU, simulate streamflow 
#(for initial phase, assume only a single HRU)

for (n in 1:num_hru) {

  n <- 1
  
  # Load data file for an HRU. Specify the FULL file name other than Lat and Lon  
  hru_climfile <- paste0(hru_lat[n], "_", hru_lon[n])
  
  griddata <- read.table(paste0(hru_filedirectory,hru_climfile))
  
  # Data extraction for the specified simulation period by "str_date" and "end_date"
  griddata_date <- griddata %>% as_data_frame() %>%
    mutate(Date = paste(V1,V2,V3, sep ="/") %>% as.Date()) %>% pull(Date)

  #Indices for the beginning and ending dates  
  sind <- which(str_date == griddata_date) 
  eind <- which(end_date == griddata_date) 

  #Precipitation and temperature time-series
  hru_prcp <- griddata[sind:eind,4]
  hru_temp <- griddata[sind:eind,5] 
  
  
  hru_par <- calib_par %>% slice(1) %>% unlist() %>% as.numeric()
  
  #Model parameters
  pars_sac <- hru_par[1:27]  # SAC_SMA parameters for a HRU

  
  #Simulate streamflow from SAC-SMA model
  hru_simflow <- sac_sma(S_Date = str_date,
                         E_Date = end_date,
                         Prcp = hru_prcp,
                         Tavg = hru_temp, 
                         Basin_Lat = hru_lat, 
                         Basin_Elev = hru_elev,
                         Par = hru_par,
                         IniState = inistates,
                         flag_snowmodule = 0)
}

