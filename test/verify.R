

# VERIFY STEP 1) --- FLOW AT EACH HRU (BEFORE ROUTING) (OK!)

#Compare model outputs
mflow <- read.table("./Data/matlab_flow_hru_all.txt") %>% as_data_frame()
rflow <- read.csv("./Data/r_hru_all_nosnow.csv") %>% select(-X) %>% as_data_frame()

for (x in 1:23) {

  print(round(sum(mflow[,x]) - sum(rflow[,x])),5)

}

# VERIFY STEP 2) --- RESULTS OF THE ROUTING MODEL   (OK!)

#Compare model outputs
mrouting <- read.table("./Data/matlab_routing_hrus.txt") %>% as_data_frame()
rrouting <- read.csv("./Data/r_routing_hrus.csv") %>% select(-X) %>% as_data_frame()

for (x in 1:23) {
  
  print(round(sum(mrouting[,x]) - sum(rrouting[,x])),5)
  
}


# VERIFY STEP 3) --- CHECK CONVOLUTION!!!!!


# Convolution is also correct (real-time check with matlab)


# VERIFY STEP 4) SNOW-MODEL CHECK


#Simulation period
sim_date <- seq.Date(as.Date("1995/10/1"), as.Date("2014/09/30"), by = "day")

flow_data <- list()

for (n in 1:23) {
  flow_data[[n]] <- data_frame(Date = sim_date, matlab = mflow[[n]], rflow = rflow2[[n]]) %>%
    mutate(diff = abs(matlab - rflow))
}


mflow[,1] - rflow[,1]

tail(mflow[,1])
tail(rflow[,1])

flow_data[[1]] %>%
  filter(diff > 1)

sum(round(flow_data[[1]]$diff,0))

View(flow_data[[1]])

flow_data_long <- flow_data %>%
  gather(key = variable, value = value, -Date)

ggplot2::ggplot(flow_data_long, aes(x = Date, y = value, color = variable)) +
  geom_line(alpha = 0.5)



################################################################################

#AFTER ROUNTING


#Compare model outputs
mflow <- read.table("./Data/matlab_totflow_nosnow.txt") %>% as_data_frame() %>%
  unlist() %>% as.numeric()
rflow <- read.csv("./Data/r_totflow_nosnow.csv") %>% select(-X) %>% as_data_frame() %>%
  unlist() %>% as.numeric()

flow_dat <- data_frame(Date = sim_date, matlab = round(mflow,3), rflow = round(rflow,3)) %>%
    mutate(diff = abs(matlab - rflow))







