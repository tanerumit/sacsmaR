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
                        tavg, prcp, lat, elev, area, flowLength, 
                        dayOfYear, lastDayOfYear) {
  
  sim_num <- length(prcp[[1]])
  hru_num <- length(tavg)

  ### CALCULATE PET FOR ALL HRUs
  pet <- lapply(1:hru_num, function(x) 
    petHamon(coeff = par_petHamon[x], tavg = tavg[[x]], lat = lat[[x]], 
             dayOfYear = dayOfYear))
  
  ### CALCULATE SNOW MODULE FOR ALL HRUs
  snowMelt <- lapply(1:hru_num, function(x) 
    snow17(par = par_snow17[x,], 
           prcp = prcp[[x]], tavg = tavg[[x]], elev = elev[[x]],
           dayOfYear = dayOfYear, lastDayOfYear =  lastDayOfYear)) 
  
  ### CALCULATE STREAMFLOW GENERATION AT ALL HRUs
  hru_simflow <- lapply(1:hru_num, function(x) 
    sacSma(par = par_sacsma[x,], 
           prcp = snowMelt[[x]], 
           pet = pet[[x]], 
           lat = lat[x], 
           elev = elev[x]))
  
  ### Fractional flow
  hru_simflow <- lapply(1:hru_num, function(x)
                        hru_simflow[[x]] * area[x] / sum(area))
  
  ### CALCULATE ROUTING
  UH_river <- lapply(1:hru_num, function(x) 
    routeLohamann(par = par_routLah[x,], flowLength = flowLength[x]))
  
  ### CONVOLUTION
  # Base Time of Unit Hydrographs of HRU and River in Routing model. 
  KE      <-  12  # Base time for HRU UH (day)
  UH_DAY  <-  96  # Base time for river routing UH (day)
  
  FLOW <- vector("numeric", length = sim_num)
  j_i <- 1:(KE+UH_DAY-1) 
  j <- lapply(1:sim_num, function(i) j_i[i-j_i+1 >= 1])

  for (n in 1:hru_num) {
    for(i in 1:sim_num) {
      FLOW[i] <- FLOW[i] + sum(UH_river[[n]][j[[i]]] * hru_simflow[[n]][i-j[[i]]+1])
    }
  }
  
  return(FLOW)

}
