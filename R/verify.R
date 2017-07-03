

#Compare model outputs
mflow <- read.table("C:\\Users\\Umit\\Desktop\\sacsmatotflow.txt") %>% unlist() %>% as.numeric()
rflow <- results[,1]

sim_date <- seq.Date(as.Date("1995/10/1"), as.Date("2014/09/30"), by = "day")


flow_data <- data_frame(Date = sim_date, matlab = mflow, rflow = rflow)

flow_data_long <- flow_data %>%
  gather(key = variable, value = value, -Date)

ggplot2::ggplot(flow_data_long, aes(x = Date, y = value, color = variable)) +
  geom_line(alpha = 0.5)


View(flow_data)

