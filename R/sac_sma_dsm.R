

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Sacramento Soil Moisture Accounting (SAC-SMA) Model
# 
# The SAC-SMA is a conceptual model using a two-layer soil moisture system to
# continuously account for storage and flow through the soil layers. 
# 
# The R version of the model is developed based on the MATLAB code (Sungwook Wi)

# Latest update: 7/5/2017

# Current progress: 
#   main sac_sma module: DONE
#   routing module: DONE
#   snow module: not implemented

# # Contact: M. Umit Taner (tanerumit@gmail.com)  
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

### Inputs:

# str_date      = starting date of the simulation, (YYYY/MM/DD)
# end_date      = ending date of the simulation, (YYYY/MM/DD)
# grid_lat      = latitude info for all HRUs, as character vector
# grid_lon      = longitude info for all HRUs, as character vector
# grid_area     = surface area of all HRUs (km2), as numeric vector  
# grid_elev     = elevation info for all HRUs (m), as numeric vector
# grid_flowlen  = flow length info for all HRUs (m),, as numeric vector
# grid_par      = calibrated parameter set for al HRUs, 31 parameter values
# flog_snow17   = flag to indicate if the snow component is active (0 or 1)

#Simulate streamflow for each HRU
#S_Date = str_date; E_Date = end_date; 
#Prcp = hru_prcp; Tavg = hru_temp; Basin_Lat = hru_lat; 
#Basin_Elev = hru_elev; Par = hru_par; IniState = inistates; 
#flag_snowmodule = 0 


#griddata_list <- lapply(1:num_hru, function(x) read.table(hru_names[x])) %>%
#  bind_rows(.id = "hru") %>% as_data_frame() %>%
#  write_csv(., path = "gridclim.csv")


sac_sma_dsm <- function(str_date, end_date, grid_lat, grid_lon, grid_area, grid_elev, 
                        grid_flowlen, grid_par, flag_SNOW17 = 0) {
  
  #Essential R-packages needed
  require(tidyr)
  require(dplyr)
  require(readr)
  require(lubridate)
  
  # Initial Storage States in SAC_SMA
  uztwc <- 0  # Upper zone tension water storage
  uzfwc <- 0  # Upper zone free water storage
  lztwc <- 5  # Lower zone tension water storage
  lzfsc <- 5  # Lower zone supplementary free water storage
  lzfpc <- 5  # Upper zone primary free water storage
  adimc <- 0  # Additional impervious area storage
  inistates <- c(uztwc, uzfwc, lztwc, lzfsc, lzfpc, adimc)
  
  # Include states for snow17 if snow module is 'on'
  if (flag_SNOW17 == 1) {
    W_i     <- 0  # Accumulated water equivalent of the ice portion of the snow cover (mm)
    ATI     <- 0  # Antecedent Temperature Index, deg C
    W_q     <- 0  # Liquid water held by the snow (mm)
    Deficit <- 0  # Heat Deficit, also known as NEGHS, Negative Heat Storage
    inistates <- c(inistates, W_i, ATI, W_q, Deficit)
  }
  
  # Base Time of Unit Hydrographs of HRU and River in Routing model. 
  # (Necessary to be changed only for a VERY LARGE watershed)
  KE      <-  12  # Base time for HRU UH (day)
  UH_DAY  <-  96  # Base time for river routing UH (day)

  # Run Distributed System Modeling version of SAC_SMA for each HRU
  sim_date <- seq.Date(str_date, end_date, by = "day") 
  sim_doy  <- yday(sim_date)
  sim_per  <- length(sim_date)      # Simulation period
  
  num_hru  <- nrow(gridinfo)        # Total number of HRUs
  tot_area <- sum(grid_area)        # Total watershed area
  
  #Lat-long information
  hru_lat  <- sapply(1:num_hru, function(n)
    formatC(grid_lat[n], format = 'f', flag='0', digits = 4))
  hru_lon  <- sapply(1:num_hru, function(n)
    formatC(grid_lon[n], format = 'f', flag='0', digits = 4))
  hru_names <-  paste0(hru_filedir,hru_lat, "_", hru_lon)
  
  griddata_list <- read_csv("gridclim.csv")

  #Run through each HRU, simulate streamflow 
  FLOW <- vector(mode = "numeric", length = sim_per) 
  for (n in 1:num_hru) {
  
    # Read-in data for current HRU
    griddata  <- dplyr::filter(griddata_list, hru == n) %>% select(-hru)
    
    # Find starting and ending indices based on the specified dates
    if (n == 1) {
      grid_date <- as.Date(paste(griddata[[1]], griddata[[2]], griddata[[3]], sep ="/"))
      grid_ind  <- which(str_date == grid_date):which(end_date == grid_date) 
    }
    
    griddata <- as.matrix(griddata)
    
    # FLOW SIMULATION USING SAC-SMA MODEL  
    hru_prcp    <- griddata[grid_ind,4]
    hru_temp    <- griddata[grid_ind,5] 

    hru_simflow <- sac_sma(Prcp = hru_prcp, Tavg = hru_temp, 
      Basin_Lat = as.numeric(hru_lat[n]), Basin_Elev = grid_elev[n], 
      Par = as.numeric(grid_par[n,]), IniState = inistates, doy = sim_doy)
    
    hru_simflow <- hru_simflow * grid_area[n] / tot_area

    # # CHANNEL ROUTING FROM LOHMANN MODEL  
    pars_rout <- as.numeric(grid_par[n,28:31])     
    UH_river  <- route_lohamann(pars = pars_rout, flowlen = grid_flowlen[n])

    # CONVOLUTION (Vectorized form)  
    for(i in 1:sim_per) {
      j <- 1:(KE+UH_DAY-1); j <- j[i-j+1 >= 1]
      FLOW[i] <- FLOW[i] + sum(UH_river[j] * hru_simflow[i-j+1])
    }

  } #close the loop for each hru
  
  return(FLOW)
  
}


