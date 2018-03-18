#' @title River channel Routing model
#'
#' @description Mimics the routing model first described by Lohmann et al. (1998), which has been used within the 
#'   Variable Infiltration Capacity model (VIC) 
#'   
#' @param par Vector of length four, describing the following model parameters (see details)
#' @param flength Flow length to the basin outlet (in meters) 
#' @param uh.day Base time for river routing Unit Hydrograph (day) (Default: 96)
#' @param ke  Base time for HRU unit hydrograph (day) (Default: 12)
#' @return Placeholder
#' @rdname lohmann
#' 
#' @details Each grid cell is represented by a node in the channel network.
#' The total runoff and baseflow from each grid cell is first convolved with a unit hydrograph 
#' representing the distribution of travel times of water from its points of origin to the channel network
#' Then, each grid cell's input into the channel network is routed through the channel 
#' using linearized St. Venant's equations
#' 
#' The four parameters of the model are as follows:
#' NumRes: Grid Unit Hydrograph parameter defining the number of linear reservoirs
#' K: Grid Unit Hydrograph parameter defining the reservoir storage constant
#' VELO: wave velocity in the linearized Saint-Venant equation (cubic meters per sec)
#' DIFF: diffusivity in the linearized Saint-Venant equation (sq.meters per sec)
#' 
#' The detailed description of the routing can be found at:
#' Lohmann, D. E. Raschke, B. Nijssen and D.P. Lettenmaier, 1998: Regional scale hydrology: II. Application of the VIC-2L model to the Weser river, Germany, Hydrol. Sci. J., 43(1), 143-158.
#' 
#' @export 
lohmann <- function(par, flength, uh.day = 96, ke = 12) {
  
  hruh_fun <- function(x) dgamma(x, shape = NumRes, scale = 1/K)
  
  # Routing model parameters
  NumRes <- par[1]  # Grid Unit Hydrograph parameter (number of linear reservoirs)
  K <- par[2]  # Grid Unit Hydrograph parameter (reservoir storage constant)
  VELO <- par[3]  # wave velocity in the linearized Saint-Venant equation(m/s)
  DIFF <- par[4]  # diffusivity in the linearized Saint-Venant equation(m2/s)
  
  DT <- 3600  # Time step in second for solving Saint-Venant equation. This will affect TMAX
  TMAX <- uh.day * 24  # Base time of river routing UH in hour because DT is for an hour
  LE <- 48 * 50
  
  # Green's function calculation to solve the linearized Saint-Venant
  UH_DAILY <- vector(mode = "numeric", length = uh.day)
  
  if (flength == 0) {
    UH_DAILY[1] <- 1
  } else {
    
    t <- 0
    uhm_grid <- vector(mode = "numeric", length = LE)
    
    for (k in 1:LE) {
      
      t <- t + DT
      pot <- ((VELO * t - flength)^2)/(4 * DIFF * t)
      
      if (pot <= 69) {
        H <- flength/(2 * t * sqrt(pi * t * DIFF)) * exp(-pot)
      } else H <- 0
      uhm_grid[k] <- H
    }
    
    if (sum(uhm_grid) == 0) 
      uhm_grid[1] <- 1 else uhm_grid <- uhm_grid/sum(uhm_grid)
    
    UHM <- uhm_grid
    
    # Derive Daily River Impulse Response Function(Green's function)
    FR <- matrix(data = 0, nrow = TMAX, ncol = 2)
    FR[1:24, 1] <- 1/24
    # >>>>>> CAN BE VECTORIZED
    for (t in 1:TMAX) {
      L <- 1:(TMAX + 24)
      L <- L[t - L > 0]
      FR[t, 2] <- FR[t, 2] + sum(FR[t - L, 1] * UHM[L])
    }
    UH_DAILY <- sapply(1:uh.day, function(t) sum(FR[(24 * t - 23):(24 * t), 2]))
  }
  
  UH_HRU <- vector(mode = "numeric", length = 12)
  for (i in 1:ke) {
    # integration in matlab: integral(fun,xmin,xmax) integration in R: integrate(f, lower, upper, ...)
    UH_HRU[i] <- integrate(hruh_fun, 24 * (i - 1), 24 * i)$value
  }
  
  # River UH >>>>>>>>>> CAN BE VECTORIZED
  UH_R <- vector(mode = "numeric", length = ke + uh.day - 1)
  for (k in 1:ke) {
    for (u in 1:uh.day) {
      UH_R[k + u - 1] <- UH_R[k + u - 1] + UH_HRU[k] * UH_DAILY[u]
    }
  }
  
  return(UH_R/sum(UH_R))
  
}
