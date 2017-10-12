#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param par PARAM_DESCRIPTION
#' @param tavg PARAM_DESCRIPTION
#' @param prcp PARAM_DESCRIPTION
#' @param lat PARAM_DESCRIPTION
#' @param elev PARAM_DESCRIPTION
#' @param dayOfyear PARAM_DESCRIPTION
#' @param lastDayOfYear PARAM_DESCRIPTION
#' @return OUTPUT_DESCRIPTION
#' @details DETAILS
#' @examples 
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @rdname sacHydroSim
#' @export 
sacHydroSim <- function(par_petHamon, par_snow17, par_sacsma, par_routLah,
                        hru_tavg, hru_prcp, hru_lat, hru_elev, hru_area, 
                        hru_flowLength, hru_subcat = NULL, dayOfYear, 
                        lastDayOfYear) {
  
  # Base Time of Unit Hydrographs of HRU and River in Routing model. 
  KE      <-  12  # Base time for HRU UH (day)
  UH_DAY  <-  96  # Base time for river routing UH (day)
  
  sim_num <- length(hru_prcp[[1]])
  hru_num <- length(hru_tavg)

  ### CALCULATE PET FOR ALL HRUs
  hru_pet <- lapply(1:hru_num, function(x) 
    petHamon(coeff = par_petHamon[x], tavg = hru_tavg[[x]], lat = hru_lat[[x]], 
             dayOfYear = dayOfYear))
  
  ### CALCULATE SNOW MODULE FOR ALL HRUs
  hru_meltNrain <- lapply(1:hru_num, function(x) 
    snow17(par = par_snow17[x,], prcp = hru_prcp[[x]], tavg = hru_tavg[[x]], 
           elev = hru_elev[[x]], dayOfYear = dayOfYear, 
           lastDayOfYear =  lastDayOfYear)) 
  
  ### CALCULATE STREAMFLOW GENERATION AT ALL HRUs
  hru_simflow <- lapply(1:hru_num, function(x) 
    sacSma(par = par_sacsma[x,], prcp = hru_meltNrain[[x]], pet = hru_pet[[x]], 
           lat = hru_lat[x], elev = hru_elev[x]))
  hru_simflow <- lapply(1:hru_num, function(x)
                        hru_simflow[[x]] * area[x] / sum(area))
  ### CHANNEL ROUTING
  UH_river <- lapply(1:hru_num, function(x) 
    routeLohamann(par = par_routLah[x,], flowLength = hru_flowLength[x]))
  
  ### CONVOLUTION
  j_i <- 1:(KE+UH_DAY-1) 
  j <- lapply(1:sim_num, function(i) j_i[i-j_i+1 >= 1])
  
  #convolution for a single outlet 
  if (is.null(hru_subat)) {

    #Vector to store the final flow computed
    FLOW <- vector("numeric", length = sim_num)
    
    for (n in 1:hru_num) {
      for(i in 1:sim_num) {
        FLOW[i] <- FLOW[i] + sum(UH_river[[n]][j[[i]]] * hru_simflow[[n]][i-j[[i]]+1])
      }
    }

  #convolution for multiple outlets  
  } else {
    
    #number of subcatchments
    subc_num <- length(unique(hru_subcat))
    
    #Vector to store the final flow computed
    FLOW <- rep(list(vector(mode = "numeric", length = sim_num)), subc_num)
    
    #Repeat for all outlet locations
    for (k in 1:subc_num) {
      
      #indices of hrus in the subcat
      hru_ind <- which(hru_subcat == k)
      
      #loop through each subcat
      for (n in 1:length(hru_ind)) {
        
        FLOW <- vector("numeric", length = sim_num)
        for (n in 1:hru_num) {
          for(i in 1:sim_num) {
            FLOW[[k]][i] <- FLOW[[k]][i] + sum(UH_river[[n]][j[[i]]] * hru_simflow[[n]][i-j[[i]]+1])
          }
        }
      }
    }
  } 
  
  return(FLOW)
}  
  
  
# Current progress: work on the multi-outlet channel routing part