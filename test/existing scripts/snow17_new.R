

load("../../_DEV/workdir/data/sacTest.Rdata")

# PET calculation
petHamon2 <- function(coeff, lat, tavg, doy) {
  
  var_theta   <- 0.2163108 + 2*atan(0.9671396 * tan(0.0086 * (doy-186))) 
  var_pi      <- asin(0.39795 * cos(var_theta)) 
  daylighthr  = 24-24/pi * acos((sin(0.8333*pi/180) + sin(lat*pi/180) * sin(var_pi))/(cos(lat*pi/180)*cos(var_pi))) 
  
  esat <- 0.611 * exp(17.27 * tavg/(237.3 + tavg))
  pet  <- coeff * 29.8 * daylighthr*(esat/(tavg + 273.2))
  return(pet)
}

# SNOW MODULE 
snow17 <- function(pars, prcp, tavg, elev, statesInput, timeMat) {

  output <- vector(mode = "numeric", length = length(prcp))
  
  #Model parameters
  SCF    <- pars[1]; PXTEMP <- pars[2]; MFMAX  <- pars[3]; MFMIN  <- pars[4] 
  UADJ   <- pars[5]; MBASE  <- pars[6]; TIPM   <- pars[7]; PLWHC  <- pars[8]; 
  NMF    <- pars[9]; DAYGM  <- pars[10]; 
  
  dtt <- 24  # time interval of temperature data
  dtp <- 24  # time interval of prcipitation data

  for (i in 1:length(prcp)) {
    
    Ta <- tavg[i]  # Air temperature at this time step [deg C]
    Pr <- prcp[i]  # Precipitation at this time step [mm]
    
    time <- timeMat[i,]
    
    #State parameters  
    if (i == 1) {
      W_i     <- statesInput[1] 
      ATI     <- statesInput[2] 
      W_q     <- statesInput[3] 
      Deficit <- statesInput[4] 
    }
    
    # FORM OF PRECIPITATION
    if (Ta <= PXTEMP) {
      # Air temperature is cold enough for snow to occur
      SNOW <- Pr; RAIN <- 0; 
    } else {
      # Air temperature is warm enough for rain
      SNOW <- 0; RAIN <- Pr; 
    }              
    
    # ACCUMULATION OF THE SNOW COVER
    Pn  <- SNOW * SCF  # Water equivalent of new snowfall [mm]
    W_i <- W_i + Pn    # Water equivalent of the ice portion of the snow cover [mm]
    E   <- 0           # Excess liquid water in the snow cover
    
    # ENERGY EXCHANGE AT SNOW/AIR SURFACE DURING NON-MELT PERIODS
    
    # Seasonal variation in the non-rain melt factor 
    DAYN <- time$jday 
    end_doy <-  time$end_doy - 365 
    days <- end_doy + 365
    N_Mar21 <- DAYN - (80 + end_doy)
  
    Sv <- (0.5*sin((N_Mar21 * 2 * pi)/days)) + 0.5        # Seasonal variation
    Av <- 1.0                                             # Seasonal variation adjustment, Av<-1.0 when lat < 54N
    Mf <- dtt/6 * ((Sv * Av * (MFMAX - MFMIN)) + MFMIN)   # Seasonally varying non-rain melt factor
    
    # New snow temperature and heat deficit from new snow
    if(Ta < 0) T_snow_new <- Ta else T_snow_new <- 0 
  
    # Change in the heat deficit due to new snowfall [mm], 80 cal/g: 
    #latent heat of fusion, 0.5 cal/g/C: specific heat of ice
    delta_HD_snow <- -(T_snow_new * Pn)/(80/0.5)    
    
    # Heat Exchange due to a temperature gradient
    # change in heat deficit due to a temperature gradient [mm]
    delta_HD_T <- NMF * dtp/6 * Mf/MFMAX * (ATI - T_snow_new)   
    
    # Update ATI[Antecedent Temperature Index]
    if(Pn > 1.5*dtp) {
      ATI <- T_snow_new       #Antecedent temperature index  
    } else {
      TIPM_dtt <- 1.0 - ((1.0 - TIPM)^(dtt/6)) 
      ATI <- ATI + TIPM_dtt * (Ta - ATI)
    }
    
    ATI <- min(ATI,0) 
  
    # SNOW MELT
    T_rain  <- max(Ta,0)    # Temperature of rain (deg C), Ta or 0C, whichever greater
    if(RAIN > 0.25 * dtp) {
      # Rain-on-Snow Melt
      stefan <- 6.12*(10^(-10))                                               # Stefan-Boltzman constant (mm/K/hr)
      e_sat  <- 2.7489*(10^8)*exp((-4278.63/(Ta+242.792)))                    # Saturated vapor pressure at Ta (mb)
      P_atm  <- 33.86*(29.9-(0.335*(elev/100))+(0.00022*((elev/100)^2.4)))    # Atmospheric pressure (mb) where elevation is in HUNDREDS of meters (this is incorrectly stated in the manual)
      term1  <- stefan * dtp * (((Ta+273)^4)-(273^4)) 
      term2  <- 0.0125 * RAIN * T_rain 
      term3  <- 8.5 * UADJ * (dtp/6) * ((0.9*e_sat - 6.11) + (0.00057*P_atm*Ta)) 
      Melt   <- term1 + term2 + term3 
      Melt   <- max(Melt,0)  
      
    } else if ((RAIN <= 0.25 * dtp) && (Ta > MBASE)) {
      # Non-Rain Melt
      Melt <- (Mf * (Ta - MBASE) * (dtp/dtt)) + (0.0125 * RAIN * T_rain) 
      Melt <- max(Melt,0) 
    
    } else {
        Melt <- 0
    }
    
    # Ripeness of the snow cover
    # W_i : water equivalent of the ice portion of the snow cover
    # W_q : liquide water held by the snow
    # W_qx: liquid water storage capacity
    # Qw  : Amount of available water due to melt and rain
    
    Deficit <- max(Deficit + delta_HD_snow + delta_HD_T, 0)    # Deficit <- heat deficit [mm]
    if(Deficit > (0.33*W_i)) {
      # limits of heat deficit
      Deficit <- 0.33 * W_i 
    } 
    
    if (Melt < W_i) {
      W_i  <- W_i-Melt  
      Qw   <- Melt + RAIN 
      W_qx <- PLWHC * W_i 
    
      if((Qw + W_q) > (Deficit + Deficit*PLWHC + W_qx)) { # THEN the snow is RIPE
        
        E   <- Qw + W_q - W_qx - Deficit - (Deficit * PLWHC)       # Excess liquid water [mm]
        W_i <- W_i + Deficit                             # W_i increases because water refreezes as heat deficit is decreased
        W_q <- W_qx + PLWHC * Deficit                    # fills liquid water capacity
        Deficit <- 0 
      
      } else if((Qw + W_q) >= Deficit) {
      
        #& [[Qw + W_q] <= [[Deficit*[1+PLWHC]] + W_qx]]            
        # THEN the snow is NOT yet ripe, but ice is being melted
      
        E <- 0 
        W_i <- W_i + Deficit    # W_i increases because water refreezes as heat deficit is decreased
        W_q <- W_q + Qw - Deficit 
        Deficit <- 0 
  
      } else if((Qw + W_q) < Deficit) {      #elseif [[Qw + W_q] < Deficit]    
      
        # THEN the snow is NOT yet ripe
        E <- 0 
        W_i <- W_i + Qw + W_q   # W_i increases because water refreezes as heat deficit is decreased
        Deficit <- Deficit - Qw - W_q 
      }    
  
    } else  {
      
      Melt <- W_i+W_q # Melt >= W_i
      W_i <- 0 
      W_q <- 0 
      Qw <- Melt + RAIN 
      E  <- Qw 
      # SWE = 0
      
    }              
  
    if(Deficit == 0) {ATI = 0}
  
    # Constant daily amount of melt which takes place at the snow-soil
    # interface whenever there is a snow cover
    if (W_i > DAYGM) {
      
      gmwlos <- (DAYGM/W_i)*W_q 
      gmslos <- DAYGM 
      gmro <- gmwlos + gmslos 
      W_i <- W_i - gmslos 
      W_q <- W_q - gmwlos 
      
      E <- E + gmro 
      SWE <- W_i + W_q 
      
    } else {
      
      gmro <- W_i + W_q 
      W_i <- 0 
      W_q <- 0  
      E <- E + gmro 
      SWE <- 0 
      
    }
    
    output[i] <- E
  }
    
  return(output)
}

# SAC-SMA SIMULATION
sacSim <- function(pars, prcp, tavg, lat, elev, iniState) {
  
  # Capacity Thresholds
  uztwm  <-  pars[1]    # Upper zone tension water capacity [mm]
  uzfwm  <-  pars[2]    # Upper zone free water capacity [mm]
  lztwm  <-  pars[3]    # Lower zone tension water capacity [mm]
  lzfpm  <-  pars[4]    # Lower zone primary free water capacity [mm]
  lzfsm  <-  pars[5]    # Lower zone supplementary free water capacity [mm]
  
  # Recession parsameters
  uzk    <-  pars[6]    # Upper zone free water lateral depletion rate [1/day]
  lzpk   <-  pars[7]    # Lower zone primary free water depletion rate [1/day]
  lzsk   <-  pars[8]    # Lower zone supplementary free water depletion rate [1/day]
  
  # Percolation
  zperc  <-  pars[9]    # Percolation demand scale parsameter [-]
  rexp   <-  pars[10]   # Percolation demand shape parsameter [-]: exponent of the percolation equation
  pfree  <-  pars[11]   # Percolating water split parsameter (decimal fraction): Fraction of water percolating from upper zone directly to lower zone free water storage
  
  # Impervious area
  pctim  <-  pars[12]   # Impervious fraction of the watershed area [decimal fraction)
  adimp  <-  pars[13]   # Additional impervious areas (decimal fraction)
  
  # Others
  riva   <-  pars[14]   # Riparsian vegetation area (decimal fraction)
  side   <-  pars[15]   # The ratio of deep recharge to channel base flow [-]
  rserv  <-  pars[16]   # Fraction of lower zone free water not transferrable to lower zone tension water (decimal fraction)
  
  # Potential Evapotranspiration Computation
  coeff  <-  pars[17]   # Proportionality Coefficient of Hamon Potential Evapotranspiration

  # Initial Storage States (SAC-SMA)
  uztwc <- iniState[1]     # Upper zone tension water storage
  uzfwc <- iniState[2]     # Upper zone free water storage
  lztwc <- iniState[3]     # Lower zone tension water storage
  lzfsc <- iniState[4]     # Lower zone supplementary free water storage
  lzfpc <- iniState[5]     # Upper zone primary free water storage
  adimc <- iniState[6]     # Additional impervious area storage

  # RESERVOIR STATE ARRAY INITIALIZATION
  simflow <- vector(mode = "numeric", length = length(prcp))

  #CALCULATE PET (HAMON EQUATION)  
  pet <- petHamon2(coeff = coeff, lat, tavg, doy)    
  
  #PERFORM HYDROLOGY SIMULATION USING THE SAC-SMA   
  thres_zero  <- 0.00001 # Threshold to be considered as zero
  parea       <- 1 - adimp - pctim

  pr <- prcp 

  for (i in 1:length(pet)) {

    ### COMPUTE ET LOSS FOR A TIME INTERVAL ++++++++++++++++++++++++++++++++++++
    
    # The amount of maximum ET given by input potential evapotranspiration
    edmnd = pet[i]  
    
    ## Compute for different compnents...
    # ET(1), ET from Upper zone tension water storage
    et1 <- edmnd * uztwc/uztwm
    red <- edmnd - et1  # residual ET demand
    uztwc <- uztwc - et1
    
    # ET(2), ET from upper zone free water storage
    et2 <- 0
    
    # in case et1 > uztws, no water in the upper tension water storage
    if (uztwc <= 0) {
      et1 <- et1 + uztwc #et1 = uztwc
      uztwc <- 0
      red <- edmnd - et1
      
      # when upper zone free water content is less than residual ET
      if (uzfwc < red) {
        
        # all content at upper zone free water zone will be gone as ET
        et2 <- uzfwc
        uzfwc <- 0
        red <- red - et2
        if (uztwc < thres_zero) {uztwc <- 0}
        if (uzfwc < thres_zero) {uzfwc <- 0}  
        
        # when upper zone free water content is more than residual ET
      } else {
        et2 <- red  # all residual ET will be gone as ET
        uzfwc <- uzfwc - et2
        red <- 0 
      }
      
      # in case et1 <= uztws, all maximum et (et1) are consumed at uztwc, 
      # so no et from uzfwc (et2=0)
    } else {
      
      # There's possibility that upper zone free water ratio exceeds 
      #upper zone tension water ratio. If so, free water is transferred to 
      #tension water storage
      
      if((uztwc / uztwm) < (uzfwc / uzfwm)) {
        uzrat = (uztwc + uzfwc) / (uztwm + uzfwm)
        uztwc = uztwm * uzrat
        uzfwc = uzfwm * uzrat
      }
      
      if(uztwc < thres_zero) {uztwc = 0}
      if(uzfwc < thres_zero) {uzfwc = 0}
      
    }
    
    # ET(3), ET from Lower zone tension water storage when residual ET > 0
    et3 <- red * lztwc / (uztwm + lztwm) #residual ET is always bigger than ET(3)
    lztwc <- lztwc - et3 
    
    # if lztwc is less than zero, et3 cannot exceed lztws
    if(lztwc < 0) {
      et3   <- et3 + lztwc  # et3 = lztwc  
      lztwc <- 0
    }
    
    # Water resupply from Lower free water storages to Lower tension water storage
    saved  <- rserv * (lzfpm + lzfsm)
    ratlzt <- lztwc / lztwm
    ratlz  <- (lztwc + lzfpc + lzfsc - saved) / (lztwm + lzfpm + lzfsm - saved)
    
    # water is first taken from supplementary water storage for resupply
    if (ratlzt < ratlz) {
      
      del <- (ratlz - ratlzt) * lztwm
      lztwc <- lztwc + del  # Transfer water from lzfss to lztws
      lzfsc <- lzfsc - del
      
      # if tranfer exceeds lzfsc then remainder comes from lzfps
      if(lzfsc < 0) {   
        lzfpc <- lzfpc + lzfsc
        lzfsc <- 0
      }
    }
    
    if(lztwc < thres_zero) {lztwc <- 0}
    
    # ET(5), ET from additional impervious (ADIMP) area
    # ????? no idea where this come from, I think there's a possibility that et5 can be negative values
    et5   <- et1 + (red + et2) * (adimc - et1 - uztwc) / (uztwm + lztwm)  
    adimc <- adimc - et5
    if(adimc < 0) { 
      #et5 cannot exceed adims
      et5 <- et5 + adimc # et5 = adimc
      adimc <- 0
    }
    et5 <- et5 * adimp
    
    # Time interval available moisture in excess of uztw requirements
    twx <- pr + uztwc - uztwm
    
    # all moisture held in uztw- no excess
    if(twx < 0) {
      uztwc <- uztwc + pr
      twx <- 0
      # moisture available in excess of uztw storage
    } else {
      uztwc = uztwm 
    }   
    
    
    # for now twx is excess rainfall after filling the uztwc
    adimc <- adimc + pr - twx
    
    # Compute Impervious Area Runoff
    roimp <- pr * pctim 
    
    # Initialize time interval sums
    sbf   <- 0  # Sum of total baseflow(from primary and supplemental storages)
    ssur  <- 0  # Sum of surface runoff
    sif   <- 0  # Sum of interflow
    sperc <- 0  # Time interval summation of percolation
    sdro  <- 0  # Sum of direct runoff from the additional impervious area
    
    # Determine computational time increments for the basic time interval
    ninc <- floor(1.0 + 0.2*(uzfwc+twx))  # Number of time increments that interval is divided into for further soil-moisture accountng
    dinc <- 1.0 / ninc                    # Length of each increment in days
    pinc <- twx / ninc                    # Amount of available moisture for each increment
    
    # Compute free water depletion fractions for the time increment 
    #(basic depletions are for one day)
    duz   <- 1 - (1 - uzk)^dinc
    dlzp  <- 1 - (1 - lzpk)^dinc
    dlzs  <- 1 - (1 - lzsk)^dinc
    
    # Start incremental for-loop for the time interval
    for(n in 1:ninc) {
      
      adsur <- 0 # Amount of surface runoff. This will be updated.
      
      # Compute direct runoff from adimp area
      ratio <- (adimc - uztwc) / lztwm
      if(ratio < 0) {ratio <- 0} 
      
      # Amount of direct runoff from the additional impervious area     
      addro <- pinc*(ratio^2)
      
      # Compute baseflow and keep track of time interval sum
      # Baseflow from free water primary storage
      bf_p <- lzfpc * dlzp 
      lzfpc <- lzfpc - bf_p
      if(lzfpc <= 0.0001) {
        
        bf_p  <- bf_p + lzfpc
        lzfpc <- 0
      }
      
      sbf <- sbf + bf_p
      
      # Baseflow from free water supplemental storage
      bf_s  <- lzfsc * dlzs  
      lzfsc <- lzfsc - bf_s
      if (lzfsc <= 0.0001) {
        bf_s <- bf_s + lzfsc
        lzfsc <- 0
      }
      
      # Total Baseflow from primary and supplemental storages
      sbf <- sbf + bf_s 
      
      # Compute PERCOLATION- if no water available then skip.
      if((pinc + uzfwc) <= 0.01) {
        uzfwc <- uzfwc + pinc
      } else {
        
        # Limiting drainage rate from the combined saturated lower zone storages
        percm <- lzfpm * dlzp + lzfsm * dlzs 
        perc <- percm * uzfwc / uzfwm
        
        # DEFR is the lower zone moisture deficiency ratio
        defr <- 1.0 - (lztwc + lzfpc + lzfsc)/(lztwm + lzfpm + lzfsm) 
        
        if(defr < 0) {defr <- 0} 
        
        perc <- perc * (1.0 + zperc * (defr^rexp))
        
        # Note. . . percolation occurs from uzfws before pav is added
        
        # Percolation rate exceeds uzfws
        if(perc >= uzfwc) {perc <- uzfwc}
        
        uzfwc <- uzfwc - perc    # Percolation rate is less than uzfws.
        
        # Check to see if percolation exceeds lower zone deficiency.
        check <- lztwc + lzfpc + lzfsc + perc - lztwm - lzfpm - lzfsm
        if(check > 0) {
          perc <- perc - check
          uzfwc <- uzfwc + check
        }    
        
        # SPERC is the time interval summation of PERC
        sperc <- sperc + perc  
        
        
        # Compute interflow and keep track of time interval sum. Note that PINC has not yet been added.
        del <- uzfwc * duz # The amount of interflow
        sif <- sif + del
        uzfwc <- uzfwc - del
        
        # Distribute percolated water into the lower zones. Tension water
        # must be filled first except for the PFREE area. PERCT is
        # percolation to tension water and PERCF is percolation going to
        # free water.
        
        perct <- perc * (1.0 - pfree)  # Percolation going to the tension water storage
        
        if((perct + lztwc) <= lztwm) {
          
          lztwc <- lztwc + perct
          percf <- 0.0 # Pecolation going to th lower zone free water storages
          
        } else {
          
          percf <- lztwc + perct - lztwm
          lztwc <- lztwm  
          
        }
        
        # Distribute percolation in excess of tension requirements among the free water storages.
        percf <- percf + (perc * pfree)
        if(percf != 0) {
          
          # Relative size of the primary storage as compared with total lower zone free water storages.
          hpl <- lzfpm / (lzfpm + lzfsm) 
          
          # Relative fullness of each storage.
          ratlp <- lzfpc / lzfpm
          ratls <- lzfsc / lzfsm
          
          # The fraction going to primary
          fracp <- hpl * 2.0 * (1.0 - ratlp) / (1.0 - ratlp + 1.0 - ratls) 
          if(fracp > 1.0) {fracp <- 1.0} 
          
          percp <- percf * fracp # Amount of the excess percolation going to primary
          percs <- percf - percp # Amount of the excess percolation going to supplemental
          lzfsc <- lzfsc + percs
          
          
          if(lzfsc > lzfsm) {
            percs <- percs - lzfsc + lzfsm
            lzfsc <- lzfsm
          }
          
          lzfpc <- lzfpc + percf - percs
          
          # Check to make sure lzfps does not exceed lzfpm
          if(lzfpc >= lzfpm) {
            excess <- lzfpc - lzfpm
            lztwc <- lztwc + excess
            lzfpc <- lzfpm
          }   
        }
        
        
        # Distribute PINC between uzfws and surface runoff
        if(pinc != 0) {
          
          # check if pinc exceeds uzfwm
          if((pinc + uzfwc) <= uzfwm) {
            
            uzfwc <- uzfwc + pinc  # no surface runoff
          } else {
            sur <- pinc + uzfwc - uzfwm # Surface runoff
            uzfwc <- uzfwm
            
            ssur = ssur + (sur * parea)
            
            # ADSUR is the amount of surface runoff which comes from
            # that portion of adimp which is not currently generating
            # direct runoff. ADDRO/PINC is the fraction of adimp
            # currently generating direct runoff.
            adsur = sur * (1.0 - addro / pinc) 
            ssur = ssur + adsur * adimp
            
          }
        }  
      }
      
      adimc <- adimc + pinc - addro - adsur
      if(adimc > (uztwm + lztwm)) {
        addro = addro + adimc - (uztwm + lztwm)
        adimc = uztwm + lztwm
      }
      
      # Direct runoff from the additional impervious area
      sdro  = sdro + (addro * adimp) 
      
      if(adimc < thres_zero) {adimc <- 0}
      
    } # END of incremental for loop
    
    
    # Compute sums and adjust runoff amounts by the area over which they are generated.
    
    # EUSED is the ET from PAREA which is 1.0 - adimp - pctim
    eused <- et1 + et2 + et3 
    sif <- sif * parea 
    
    # Separate channel component of baseflow from the non-channel component
    tbf <- sbf * parea   # TBF is the total baseflow
    bfcc <- tbf * (1.0 / (1.0 + side))     # BFCC is baseflow, channel component  
    
    # Ground flow and Surface flow
    base <- bfcc                       # Baseflow and Interflow are considered as Ground inflow to the channel
    surf <- roimp + sdro + ssur + sif  # Surface flow consists of Direct runoff and Surface inflow to the channel
    
    # ET(4)- ET from riparian vegetation.
    et4 <- (edmnd - eused) * riva  # no effect if riva is set to zero
    
    # Compute total evapotransporation - TET
    eused <- eused * parea 
    tet <- eused + et4 + et5 
    
    # Check that adims >= uztws
    if(adimc < uztwc) {adimc <- uztwc} 
    
    # Total inflow to channel for a timestep 
    tot_outflow <- surf + base 
    
    tot_outflow <- tot_outflow - et4 
    if(tot_outflow < 0) {tot_outflow <- 0} 
    
    # Main output: Watershed out flow
    simflow[i] <- tot_outflow 
  }
  
  return(simflow)
}
  
  
  
  











W_i     <- iniState[7]   # Accumulated water equivalent of the ice portion of the snow cover (mm)
ATI     <- iniState[8]   # Antecedent Temperature Index, deg C
W_q     <- iniState[9]   # Liquid water held by the snow (mm)
Deficit <- iniState[10]  # Heat Deficit, also known as NEGHS, Negative Heat Storage
snow_state <- c(W_i, ATI, W_q, Deficit)

## Testing the funciton
meltRain1 <- sapply(1:11688, function(i)
  SNOW17(pars = snowpar, 
         prcp = prec_lst[[n]][i], 
         temp = temp_lst[[n]][i], 
         elev = elev, 
         states_input = snow_state, 
         TIME = sim_doy_mat[i,])$meltNrain)


meltRain2 <- snow17(pars = snowpar, 
                    prcp = prec_lst[[n]], 
                    tavg = temp_lst[[n]], 
                    elev = elev, 
                    statesInput = snow_state, 
                    timeMat = sim_doy_mat)
sum(meltRain1)
sum(meltRain2)

all.equal(meltRain1, meltRain2)
