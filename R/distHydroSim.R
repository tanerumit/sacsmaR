
str.date = str_date
end.date = end_date
hru.par =  hru_par
hru.info = hru_info
hru.elevband = hru_elevband
clim.dir = clim_dir
snow.flag  = 0


HydroSimDist <- function(str.date = NULL, end.date = NULL, 
  hru.par = NULL, hru.info = NULL, hru.elevband = NULL, clim.dir, snow.flag = 0,
  progress.bar = TRUE) {
  
  # date vectors
  sim_date <- seq.Date(str.date, end.date, by = "day") 
  jdate <- as.numeric(format(sim_date, "%j"))
  ydate <- as.POSIXlt(sim_date)$year + 1900
  ddate <- as.POSIXlt(as.Date(paste0(ydate,"/12/31")))$yday

  # Simulation parameters
  sim_num  <- length(sim_date)  # number of simulation steps (days)
  hru_num  <- nrow(hru.par)  # number of HRUs 
  
  # climate data: YYYY - MM -DD -Precip (mm) - Temp (Deg C)
  # Output list ordered in the order of hru.info
  hru_clim0 <- list()
  for (n in 1:hru_num) {
    lat  <- formatC(as.numeric(hru.info[n,1]), format = 'f', flag='0', digits = 5)
    lon  <- formatC(as.numeric(hru.info[n,2]), format = 'f', flag='0', digits = 5)
    hru_clim0[[n]] <- as.matrix(read.table(paste0(clim.dir,lat, "_", lon)))
    colnames(hru_clim0[[n]]) <- c("year", "month", "day", "prcp", "tavg")
  }

  # Extract the sim period from the climate data for each HRU
  hru_date <- as.Date(paste(hru_clim0[[n]][,1], hru_clim0[[n]][,2], hru_clim0[[n]][,3], sep ="/"))
  grid_ind <- which(str.date == hru_date):which(end.date == hru_date) 
  hru_clim <- lapply(hru_clim0, function(x) as_tibble(x[grid_ind, ]))
  hru_prcp <- lapply(hru_clim, '[[', "prcp")
  hru_tavg <- lapply(hru_clim, '[[', "tavg")

  # Optimal model calibration parameters
  par_sacsma   <- hru.par[,1:16]   # Sac-sma model
  par_petHamon <- hru.par[,17]     # Hamon equation (PET)
  par_snow17   <- hru.par[,18:27]  # Snow17 model
  par_routLah  <- hru.par[,28:31]  # Lohmann routing model
  
  # HRU parameters
  hru_num   <- nrow(hru.info)  # total number of hrus in the watershed
  hru_lat     <- hru.info[, 1] # Latitude of each HRU (Deg)
  hru_lon     <- hru.info[, 2] # Longitude of each HRU (Deg)
  hru_area    <- hru.info[, 3] # Area of each HRU (as %)
  hru_elev    <- hru.info[, 4] # Elevation of each HRU (m)
  hru_flowlen <- hru.info[, 5] # Flow length for each HRU (m)
  tot_area    <- sum(hru_area) # Total area of the watershed system
  
  #  Elevation band parameters
  if (!is.null(hru.elevband)) {
    elevband_num         <- (dim(hru.elevband)[2]-2)/2
    hru_elevband_frac     <- hru.elevband[,3:(elevband_num +2)]
    hru_elevband_elev_avg <- hru.elevband[,(elevband_num +3):ncol(hru.elevband)]
  }

  #Final surfaceflow and baseflow at the outlet
  FLOW_SURF <- vector(mode = "numeric", length = sim_num)   
  FLOW_BASE <- vector(mode = "numeric", length = sim_num)   

  # Loop through each HRU and calculate adjust temp and pet values  
  if(progress.bar == TRUE) pb <- txtProgressBar(min=0, max=hru_num, style=3)
  
  for (h in 1:hru_num) {
    
    # Vectors for surface and baseflow generated at each hru 
    flow_direct_tot <- vector(mode = "numeric", length = sim_num)
    flow_base_tot   <- vector(mode = "numeric", length = sim_num)

    # If elevation bands are not specified
    if (is.null(hru.elevband)) {
    
      # PET using Hamon equation
      hru_pet <- hamon(par = par_petHamon[h], tavg = hru_tavg[[h]], lat = hru_lat[h], jdate = jdate)

      # If snow module is active, calculate snow
      if (snow.flag == 1) {
        hru_mrain <- snow17(par_snow17, hru_prcp[[h]], hru_tavg[[h]], hru_elev[h], jdate) 
      } else { 
        hru_mrain <- hru_prcp[[h]]
      }
 
      # Calculate hru flow for each elevation band  
      out <- sacSma(par = par_sacsma[h,], prcp = hru_mrain, pet = hru_pet)
      flow_direct_tot <- out$surf 
      flow_base_tot   <- out$base 

    } else {

      #Loop through each elevation band
      for (nn in 1:elevband_num) {
      
        if(hru_elevband_frac[h, nn] > 0) {
          
          # Adjusted temperature 
          hru_tavg_adj <- hru_tavg[[h]] - 6.1 * (hru_elev[[h]] - hru_elevband_elev_avg[h, nn]) / 1000
          
          # Adjusted PET
          hru_pet <- hamon(par = par_petHamon[h], tavg = hru_tavg_adj, lat = hru_lat[h], jdate = jdate)
          
          # If snow module is active, calculate snow
          if (snow.flag == 1) {
            hru_mrain <- snow17(par_snow17, hru_prcp[[h]], hru_tavg_adj, hru_elev[h], jdate) 
          } else { 
            hru_mrain <- hru_prcp[[h]]
          }
          
          # Calculate hru flow for each elevation band  
          out <- sacSma(par = par_sacsma[h,], prcp = hru_mrain, pet = hru_pet)
          flow_direct_tot <- flow_direct_tot + out$surf * hru_elevband_frac[h,nn]
          flow_base_tot   <- flow_base_tot   + out$base * hru_elevband_frac[h,nn]

        }

      } # close loop elev. bands

    }

    #---------- Channel flow routing from Lohmann routing model ---------------#  
    out2 <- lohamann(par = par_routLah[h,], inflow.direct = flow_direct_tot, 
      inflow.base = flow_base_tot, flowlen = hru_flowlen[h])

    FLOW_SURF <- FLOW_SURF + out2$surf * hru_area[h] / tot_area
    FLOW_BASE <- FLOW_BASE + out2$base * hru_area[h] / tot_area
    
    
    if(progress.bar == TRUE) setTxtProgressBar(pb, h)
  }
  
  if(progress.bar == TRUE) close(pb)
  return(list(SURF = FLOW_SURF, BASE = FLOW_BASE))
}  
  







