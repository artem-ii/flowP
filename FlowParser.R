setwd("/Volumes/BigHDD/Inserm/Exps/RoraFoxp3_lung1/14-09-17")
library(readxl)
file <- "SPLEEN(Lung1)-14-09-17.xls"
all_data <- read_excel(file,
                       col_types = c("text", "numeric", "numeric",
                                     "numeric", "numeric", "numeric",
                                     "numeric", "numeric", "numeric",
                                     "numeric", "numeric", "numeric",
                                     "numeric", "numeric", "numeric",
                                     "numeric", "numeric", "numeric",
                                     "numeric", "numeric", "numeric",
                                     "numeric", "numeric", "numeric",
                                     "numeric", "numeric", "numeric",
                                     "numeric", "numeric", "numeric",
                                     "numeric", "numeric", "numeric",
                                     "numeric", "numeric", "numeric",
                                     "numeric", "numeric", "numeric", "numeric",
                                     "numeric", "numeric", "numeric", "numeric",
                                     "numeric", "numeric", "numeric",
                                     "numeric", "numeric", "numeric",
                                     "numeric", "numeric", "numeric",
                                     "numeric", "numeric", "numeric",
                                     "numeric", "numeric", "numeric",
                                     "numeric", "numeric", "numeric",
                                     "numeric", "numeric", "numeric",
                                     "numeric", "numeric", "numeric",
                                     "numeric", "numeric", "numeric",
                                     "numeric", "numeric", "numeric",
                                     "numeric", "numeric", "numeric", "numeric", "numeric",
                                     "numeric", "numeric", "numeric", "numeric",
                                     "numeric", "numeric", "numeric",
                                     "numeric", "numeric", "numeric",
                                     "numeric", "numeric", "numeric",
                                     "numeric", "numeric", "numeric",
                                     "numeric", "numeric", "numeric",
                                     "numeric", "numeric", "numeric",
                                     "numeric", "numeric", "numeric",
                                     "numeric", "numeric", "numeric",
                                     "numeric", "numeric", "numeric",
                                     "numeric", "numeric", "numeric",
                                     "numeric", "numeric", "numeric"))
names(all_data)[1] <- "filename"
names(all_data) <- gsub("CD4 ", "CD4+", names(all_data))
names(all_data) <- gsub("T-bet", "Tbet", names(all_data))
mfis <- names(all_data)[grep("Mean", x = names(all_data))]


df <- popConstruct("/Users/artemii/FlowParser/pops.txt", all_data, verbose=T)
getwd()
l <- tempdf[1]
names(tempdf[[1]]) <- "ohh"
all_data <- all_data[grep("Mean|SD", all_data$filename, invert = T),]
