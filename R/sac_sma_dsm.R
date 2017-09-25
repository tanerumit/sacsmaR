

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


sac_sma_dsm <- function(str_date, end_date, grid_lat, grid_lon, grid_area, 
                        grid_elev, grid_par, flag_SNOW17 = 0) {
  
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
  sim_per  <- length(sim_date) # Simulation period
  
  num_hru  <- nrow(grid_info)  # Total number of HRUs
  tot_area <- sum(grid_area)   # Total watershed area
  
  #Simulate streamflow at each hru +++++++++++++++++++++++++++++++++++++++++++++
  hru_simflow <- list()
  for (n in 1:num_hru) {
  
    # Read-in data for current HRU
    hru_data  <- hru_clim[[n]]

    # Find starting and ending indices based on the specified dates
    if (n == 1) {
      grid_date <- as.Date(paste(hru_data[[1]], hru_data[[2]], hru_data[[3]], sep ="/"))
      grid_ind  <- which(str_date == grid_date):which(end_date == grid_date) 
    }
    
    hru_data <- as.matrix(hru_data)
    
    # FLOW SIMULATION USING SAC-SMA MODEL  
    hru_prcp    <- hru_data[grid_ind,4]
    hru_temp    <- hru_data[grid_ind,5] 

    hru_simflow[[n]] <- sac_sma(Prcp = hru_prcp, Tavg = hru_temp, 
      Basin_Lat = grid_lat[n], Basin_Elev = grid_elev[n], 
      Par = as.numeric(grid_par[n,]), IniState = inistates, doy = sim_doy,
      flag_snowmodule = flag_SNOW17)
    
    hru_simflow[[n]] <- hru_simflow[[n]] * grid_area[n] / tot_area
    
  } #close the loop for each hru
  
  #Rout simulated streamflow at each outlet location +++++++++++++++++++++++++++
  outlet_loc <- unique(grid_info_multigage$outlet_id)
  outlet_num <- length(outlet_loc)
  FLOW_ALL <- rep(list(vector(mode = "numeric", length = sim_per)), outlet_num)
  
  #Repeat for all outlet locations
  for (k in 1:outlet_num) {
  
    FLOW <- vector(mode = "numeric", length = sim_per)
    
    grid_info_cur <- grid_info_multigage %>% filter(outlet_id == k)
    num_hru_cur   <- grid_info_cur %>% pull(HRU_id) 
    grid_par_cur  <- grid_par[num_hru_cur,]
    grid_flowlen_cur <- grid_info_cur$`HRU_FlowLen(m)`
    
    for (n in 1:length(num_hru_cur)) {
      
      # CHANNEL ROUTING FROM LOHMANN MODEL 
      pars_rout <- as.numeric(grid_par_cur[n,28:31])     
      UH_river  <- route_lohamann(pars = pars_rout, flowlen = grid_flowlen_cur[n], UH_DAY, KE)
      
      # CONVOLUTION (Vectorized form)  
      for(i in 1:sim_per) {
        j <- 1:(KE+UH_DAY-1) 
        j <- j[i-j+1 >= 1]
        FLOW[i] <- FLOW[i] + sum(UH_river[j] * hru_simflow[[n]][i-j+1])
      }
    }  

    FLOW_ALL[[k]] <- FLOW  
    
  } #close loop for outlet locations
  
  return(FLOW_ALL)
}




