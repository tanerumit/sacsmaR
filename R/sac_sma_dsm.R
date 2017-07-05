

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



sac_sma_dsm <- function(str_date, end_date, grid_lat, grid_lon, grid_area, grid_elev, 
                        grid_flowlen, grid_par, flag_SNOW17 = 0) {
  
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
  sim_per  <- length(sim_date)      # Simulation period
  num_hru  <- nrow(gridinfo)        # Total number of HRUs
  tot_area <- sum(grid_area)        # Total watershed area
  
  #Run through each HRU, simulate streamflow 
  for (n in 1:num_hru) {
  
    # Read-in data for current HRU
    hru_lat <- formatC(grid_lat[n], format = 'f', flag='0', digits = 5)
    hru_lon <- formatC(grid_lon[n], format = 'f', flag='0', digits = 5)  
    griddata  <- read.table(paste0(hru_filedir,hru_lat, "_", hru_lon))
    
    # Find starting and ending indices based on the specified dates
    grid_date <- griddata %>% as_data_frame() %>%
      mutate(Date = paste(V1,V2,V3, sep ="/") %>% as.Date()) %>% pull(Date)
    sind <- which(str_date == grid_date) 
    eind <- which(end_date == grid_date) 
    
    # FLOW SIMULATION USING SAC-SMA MODEL --------------------------------------
    #Set relevant precipitation and temperature time-series

    hru_prcp    <- griddata[sind:eind,4]
    hru_temp    <- griddata[sind:eind,5] 
    
    hru_area    <- grid_area[n]
    hru_elev    <- grid_elev[n]
    hru_flowlen <- grid_flowlen[n]
    hru_par     <- grid_par[n,] %>% as.numeric()
  
    #Simulate streamflow for each HRU
    hru_simflow <- sac_sma(S_Date = str_date, E_Date = end_date, 
      Prcp = hru_prcp, Tavg = hru_temp, Basin_Lat = hru_lat, 
      Basin_Elev = hru_elev, Par = hru_par, IniState = inistates, 
      flag_snowmodule = 0)
    
    hru_simflow <- hru_simflow * hru_area / tot_area
   
    # CHANNEL ROUTING FROM LOHMANN MODEL ---------------------------------------
    pars_rout <- hru_par[28:length(hru_par)] # Routing model parameters
    UH_river  <- route_lohamann(pars = pars_rout, flowlen = hru_flowlen, KE = KE, UH_DAY = UH_DAY)
   
    # MAKE CONVOLUTION FOR BASIN OUTLET STREAMFLOW -----------------------------
    FLOW <- matrix(0, nrow = sim_per, ncol = num_hru) #Simulated flow for each hru

    #Loop through each period
    for(i in 1:sim_per) {
      #Loop through "base time periods"
      for(j in 1:(KE+UH_DAY-1)) {
        if((i-j+1) >= 1) {
          FLOW[i,n] <- FLOW[i,n] + UH_river[j] * hru_simflow[i-j+1]
        }
      }
    }

  } #close the loop for each hru
  
  #Find flow at the outlet by summing flow over HRUs
  TOTFLOW <- apply(FLOW, 1, sum)
  
  return(TOTFLOW)
}

#DF <- data.frame(x = 1:sim_per, y = FLOW[,1])
#ggplot(DF, aes(x,y)) + geom_line()
  