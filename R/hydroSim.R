#' @title Distributed Hydrology-Routing Model
#' @description To be completed...
#' @param par.hamon   parameters passed to hamon module 
#' @param par.snow17  parameters passed to snow17 module
#' @param par.sacsma  parameters passed to sacsma module
#' @param par.lohmann parameters to be passed to lohmann module
#' @param tavg.grid   a list of temperature series for each grid point
#' @param prcp.grid   a list of precipitation series for each grid point
#' @param lat.grid    a list of latitudes for each grid point
#' @param elev.grid   a list of elevation information for each grid point
#' @param area.grid   a list of area information for each grid point
#' @param flength.grid a list of flowlength information for each grid point
#' @param subcat.grid a list of subcatchment outlet points
#' @param jday a time-series of julian day of the year (1-3655)
#' @return flow series at the specified outlet point(s)
#' @details to be completed...
#' @rdname hydroSim
#' @export 
hydroSim <- function(par.hamon, par.snow17, par.sacsma, par.lohmann, tavg.grid, prcp.grid, lat.grid, 
  elev.grid, area.grid, flength.grid, subcat.grid = NULL, jday) {
  
  # Base Time of Unit Hydrographs of HRU and River in Routing model.
  KE <- 12  # Base time for HRU UH (day)
  UH_DAY <- 96  # Base time for river routing UH (day)
  
  sim_num <- length(prcp.grid[[1]])
  hru_num <- length(tavg.grid)
  
  ### CALCULATE PET FOR ALL HRUs
  hru_pet <- lapply(1:hru_num, function(x) hamon(par = par.hamon[x], tavg = tavg.grid[[x]], 
    lat = lat.grid[[x]], jday))
  
  ### CALCULATE SNOW MODULE FOR ALL HRUs
  hru_meltNrain <- lapply(1:hru_num, function(x) snow17(par = par.snow17[x, ], prcp = prcp.grid[[x]], 
    tavg = tavg.grid[[x]], elev = elev.grid[[x]], jday = jday))
  
  ### CALCULATE STREAMFLOW GENERATION AT ALL HRUs
  hru_simflow <- lapply(1:hru_num, function(x) sacsma(par = par.sacsma[x, ], prcp = hru_meltNrain[[x]], 
    pet = hru_pet[[x]], lat = lat.grid[x], elev = elev.grid[x]))
  
  hru_simflow <- lapply(1:hru_num, function(x) hru_simflow[[x]] * area.grid[x]/sum(area.grid))
  
  ### CHANNEL ROUTING
  UH_river <- lapply(1:hru_num, function(x) lohmann(par = par.lohmann[x, ], flength = flength.grid[x]))
  
  ### CONVOLUTION
  j_i <- 1:(KE + UH_DAY - 1)
  j <- lapply(1:sim_num, function(i) j_i[i - j_i + 1 >= 1])
  
  # convolution for a single outlet
  if (is.null(subcat.grid)) {
    
    # Vector to store the final flow computed
    FLOW <- vector("numeric", length = sim_num)
    
    for (n in 1:hru_num) {
      for (i in 1:sim_num) {
        FLOW[i] <- FLOW[i] + sum(UH_river[[n]][j[[i]]] * hru_simflow[[n]][i - j[[i]] + 1])
      }
    }
    
    # convolution for multiple outlets
  } else {
    
    # number of subcatchments
    subc_num <- length(unique(subcat.grid))
    
    # Vector to store the final flow computed
    FLOW <- rep(list(vector(mode = "numeric", length = sim_num)), subc_num)
    
    # Repeat for all outlet locations
    for (k in 1:subc_num) {
      
      # indices of hrus in the subcat
      hru_ind <- which(subcat.grid == k)
      
      # loop through each subcat
      for (n in 1:length(hru_ind)) {
        
        FLOW <- vector("numeric", length = sim_num)
        for (n in 1:hru_num) {
          for (i in 1:sim_num) {
          FLOW[[k]][i] <- FLOW[[k]][i] + sum(UH_river[[n]][j[[i]]] * hru_simflow[[n]][i - j[[i]] + 1])
          }
        }
      }
    }
  }
  
  return(FLOW)
}