
#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param par PARAM_DESCRIPTION
#' @param flowLength PARAM_DESCRIPTION
#' @param UH_DAY PARAM_DESCRIPTION, Default: 96
#' @param KE PARAM_DESCRIPTION, Default: 12
#' @return OUTPUT_DESCRIPTION
#' @details DETAILS
#' @examples 
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @rdname routeLohamann
#' @export 
routeLohamann <- function(par, flowLength, UH_DAY = 96, KE = 12) {
  
  hruh_fun <- function(x) dgamma(x, shape = NumRes, scale = 1/K)
  
  # Routing model parameters
  NumRes <- par[1]       # Grid Unit Hydrograph parameter (number of linear reservoirs)
  K      <- par[2]       # Grid Unit Hydrograph parameter (reservoir storage constant)
  VELO   <- par[3]       # wave velocity in the linearized Saint-Venant equation(m/s)
  DIFF   <- par[4]       # diffusivity in the linearized Saint-Venant equation(m2/s)
  
  DT    <- 3600          # Time step in second for solving Saint-Venant equation. This will affect TMAX
  TMAX  <- UH_DAY * 24   # Base time of river routing UH in hour because DT is for an hour
  LE		<- 48*50
  
  # Green's function calculation to solve the linearized Saint-Venant
  UH_DAILY <- vector(mode = "numeric", length = UH_DAY)
  
  if (flowLength == 0) {
    UH_DAILY[1] <- 1
    
  } else {
    
    t <- 0 
    uhm_grid <- vector(mode = "numeric",length = LE)
    
    for(k in 1:LE) {
      
      t <- t + DT 
      pot <- ((VELO*t-flowLength)^2)/(4*DIFF*t)
      
      if(pot <= 69) {
        H <- flowLength/(2*t*sqrt(pi*t*DIFF))*exp(-pot) } 
      else H <- 0
      uhm_grid[k] <- H 
    }
    
    if(sum(uhm_grid) == 0) uhm_grid[1] <- 1.0 
    else uhm_grid <- uhm_grid/sum(uhm_grid)  
    
    UHM <- uhm_grid 
    
    # Derive Daily River Impulse Response Function(Green's function)
    FR <- matrix(data = 0, nrow = TMAX, ncol = 2) 
    FR[1:24,1] <- 1/24   
    # >>>>>> CAN BE VECTORIZED
    for(t in 1:TMAX) {
      L <- 1:(TMAX+24); L <- L[t-L > 0]
      FR[t,2] <- FR[t,2] + sum(FR[t - L,1] * UHM[L])  
    }
    UH_DAILY <- sapply(1:UH_DAY, function(t) sum(FR[(24*t-23):(24*t),2]))
  }
  
  UH_HRU <- vector(mode = "numeric", length = 12)
  for(i in 1:KE) {
    #integration in matlab: integral(fun,xmin,xmax)
    #integration in R: integrate(f, lower, upper, ...)
    UH_HRU[i] <- integrate(hruh_fun, 24*(i-1), 24*i)$value 
  }
  
  # River UH 
  # >>>>>>>>>> CAN BE VECTORIZED
  UH_R <- vector(mode = "numeric", length = KE+UH_DAY-1)
  for(k in 1:KE) {
    for(u in 1:UH_DAY) {
      UH_R[k+u-1] <- UH_R[k+u-1] + UH_HRU[k] * UH_DAILY[u] 
    }
  }
  
  return(UH_R/sum(UH_R)) 
  
}
