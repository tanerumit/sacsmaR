

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
  num_hru  <- nrow(gridinfo)          # Total number of HRUs
  tot_area <- sum(grid_area)           # Total watershed area
  sim_date <- seq.Date(str_date, end_date, by = "day") 
  sim_per  <- length(sim_date)        # Simulation period

  #EMPTY ARRAYS FOR THE OUTPUTS 
  FLOW <- rep(0, sim_per) #Simulated flow for each hru
  
  #Run through each HRU, simulate streamflow 
  for (n in 1:num_hru) {
  
    #Get latitude & longitude of current hru
    hru_lat <- formatC(grid_lat[n], format = 'f', flag='0', digits = 5)
    hru_lon <- formatC(grid_lon[n], format = 'f', flag='0', digits = 5)  

    # Get climate datafile for the HRU (YYY MM DD Precip Temp)
    griddata  <- read.table(paste0(hru_filedirectory,hru_lat, "_", hru_lon))
    
    # Data extraction for the specified simulation period by "str_date" and "end_date"
    grid_date <- griddata %>% as_data_frame() %>%
      mutate(Date = paste(V1,V2,V3, sep ="/") %>% as.Date()) %>% pull(Date)
  
    #Indices for the beginning and ending dates  
    sind <- which(str_date == grid_date) 
    eind <- which(end_date == grid_date) 
  
    #Precipitation and temperature time-series
    hru_prcp    <- griddata[sind:eind,4]
    hru_temp    <- griddata[sind:eind,5] 
    
    #HRU specific data
    hru_area    <- grid_area[n]
    hru_elev    <- grid_elev[n]
    hru_flowlen <- grid_flowlen[n]
    hru_par     <- grid_par[n,] %>% as.numeric()
  
    #Simulate streamflow for each HRU
    hru_simflow <- sac_sma(S_Date = str_date, E_Date = end_date, Prcp = hru_prcp, 
      Tavg = hru_temp, Basin_Lat = hru_lat, Basin_Elev = hru_elev, Par = hru_par, 
      IniState = inistates, flag_snowmodule = 0)
    
    hru_simflow <- hru_simflow * hru_area[n] / tot_area

    # CHANNEL ROUTING FROM LOHMANN MODEL ---------------------------------------
    pars_rout <- hru_par[28:length(hru_par)] # Routing model parameters
    
    UH_river  <- route_lohamann(pars = pars_rout, 
      flowlen = hru_flowlen, KE = KE, UH_DAY = UH_DAY)

    # MAKE CONVOLUTION FOR BASIN OUTLET STREAMFLOW -----------------------------
    for (i in 1:sim_per) {
      for(j in 1:KE+UH_DAY-1) {
        if(i-j+1 >= 1) {
          FLOW[i] = FLOW[i] + UH_river[j] * hru_simflow[i-j+1]
        }
      }
    }
    
  } #close the loop for each hru

  return(FLOW)
}