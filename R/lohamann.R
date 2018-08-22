#' Title
#'
#' @param par 
#' @param inflow.direct 
#' @param inflow.base 
#' @param flowlen 
#'
#' @return
#' @export
#'
#' @examples
lohamann <- function(par, inflow.direct, inflow.base, flowlen) {
  
  #browser()
  
  #------------------------------ Lohamann parameters -------------------------#
  N    <- par[1]  # Grid Unit Hydrograph parameter (number of linear reservoirs)
  K    <- par[2]  # Grid Unit Hydrograph parameter (reservoir storage constant)
  VELO <- par[3]  # wave velocity in the linearized Saint-Venant equation(m/s)
  DIFF <- par[4]  # diffusivity in the linearized Saint-Venant equation(m2/s)
  #----------------------------------------------------------------------------#

  #------- Base Time for HRU(watershed subunit) UH and channel routing UH -----#
  KE     <- 12 
  UH_DAY <- 96 
  DT     <- 3600        # Time step in second for solving Saint-Venant equation. This will affect TMAX
  TMAX   <- UH_DAY*24   # Base time of river routing UH in hour because DT is for an hour
  LE		 <- 48*50
  #----------------------------------------------------------------------------#
  
  #----- Derive Daily River Impulse Response Function(Green's function) -------#
  UH_river <- vector(mode = "numeric", length = UH_DAY)

  # if the HUR contains the watershed outlet
  if (flowlen == 0) {
    
    UH_river[1] <- 1
  
  # if not, calculate Green's function to solve the linearized Saint-Venant Eq
  } else {
    
    t <- 0 
    
    uhm_grid <- vector(mode = "numeric",length = LE)
    
    for(k in 1:LE) {
      
      t <- t + DT 
      pot <- ((VELO*t-flowlen)^2)/(4*DIFF*t)
      
      if(pot <= 69) {
        H <- flowlen/(2*t*sqrt(pi*t*DIFF))*exp(-pot) 
      } else {H <- 0}
      
      uhm_grid[k] <- H 
    }
    
    if(sum(uhm_grid) == 0) {
      
      uhm_grid[1] <- 1 
      
    } else {
      
      uhm_grid <- uhm_grid/sum(uhm_grid)
    
    }  
    
    UHM <- uhm_grid 
    
    # Derive Daily River Impulse Response Function(Green's function)
    FR <- matrix(data = 0, nrow = TMAX, ncol = 2) 
    FR[1:24,1] <- 1/24   

    for(t in 1:TMAX) {
      
      ################## CHECK THIS !!!!!!!!!!!!!!!
      L <- 1:(TMAX+24) 
      L <- L[t-L > 0]
      
      FR[t,2] <- FR[t,2] + sum(FR[t - L,1] * UHM[L])  
    }

    UH_river <- sapply(1:UH_DAY, function(t) sum(FR[(24*t-23):(24*t),2]))
  }
  
  #----------------------------------------------------------------------------#
  # HRU's Unit Hydrograph represented by Gamma distribution function
  
  hruh_fun <- function(x) dgamma(x, shape = N, scale = 1/K)
  UH_HRU_direct <- vector(mode = "numeric", length = 12)
  
  for(i in 1:KE) {
    UH_HRU_direct[i] <- integrate(hruh_fun, 24*(i-1), 24*i)$value
  }
  
  UH_HRU_base <-  vector(mode = "numeric", length = KE) 
  UH_HRU_base[1] <- 1
  
  #----------------------------------------------------------------------------#
  
  #-------- Combined UH for HRU's response at the watershed outlet ------------#
  
  UH_direct <- vector(mode = "numeric", length = KE+UH_DAY-1)
  UH_base   <- vector(mode = "numeric", length = KE+UH_DAY-1)
  
  for(k in 1:KE) {
    
    for (u in 1:UH_DAY) {
      
      UH_direct[k+u-1] <- UH_direct[k+u-1] + UH_HRU_direct[k] * UH_river[u] 
      UH_base[k+u-1]   <- UH_base[k+u-1]   + UH_HRU_base[k]   * UH_river[u] 
    }
  }
  
  UH_direct <- UH_direct / sum(UH_direct) 
  UH_base   <- UH_base / sum(UH_base) 
  
  #----------------------------------------------------------------------------#  
  
  #------------ Make Convolution for watershed outlet total flow --------------#
  
  directflow <- vector(mode = "numeric", length = length(inflow.direct))
  baseflow   <- vector(mode = "numeric", length = length(inflow.direct))
  
  for (i in 1:length(inflow.direct)) {
    
      j <- 1:(KE+UH_DAY-1)
      j <- j[i-j+1 >= 1]
      
      directflow[i] <- directflow[i] + sum(UH_direct[j] * inflow.direct[i-j+1])
      baseflow[i] <- baseflow[i] + sum(UH_base[j] * inflow.base[i-j+1])
  }

  return(list(surf = directflow + baseflow, base = baseflow))
}



  ##codegen
  # Routing model for Land Surface model outputs based on Lohmann routing
  # model
  # 
  # Catchment's UH is represented by the Gamma distribution 
  # River routing is based on the linearized Saint-Venant Eq 
  # 
  # Inputs
  #	inflow.direct  =  Direct runoff from catchment (surface runoff, prompt subsurface runoff)
  #	inflow.base    =  Baseflow from catchment (groundwater runoff, delayed subsurface runoff)
  #   flowlen        =  Travel distance of water from catchment outlet to watershed outlet
  #	par      =  UH & Lohmann routing parameters
  #	isOutlet       =  Indicator telling whether the catchment is for watershed outlet 
  #	
  # Parameters
  #	par(1)  =  Catchment's UH shape parameter (N)
  #	par(2)  =  Catchment's UH scale parameter (K)
  #	par(3)  =  Wave velocity in the linearized Saint-Venant equation(m/s) (VELO)
  #	par(4)  =  Diffusivity in the linearized Saint-Venant equation(m2/s) (DIFF)
  # 
  # Outputs
  # runoff   = Runoff at the watershed outlet 
  # baseflow = Baseflow at the watershed outlet 