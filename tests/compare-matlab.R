

watershed <- "arrhon"

streamflow_mlab <- read.table(file = paste0("./output/sacsma_matlab_", watershed,".txt"), col.names = FALSE) %>% unlist(use.names = FALSE)
streamflow_r    <- read.table(file = paste0("./output/sacsma_r_", watershed,".txt"), col.names = FALSE) %>% unlist(use.names = FALSE)

df <- data_frame(id = 1:length(streamflow_mlab), matlab = streamflow_mlab, R = streamflow_r)
df <- df[-(1:(365*2)),]

ggplot(df, aes(x = matlab, y = R)) + geom_point()