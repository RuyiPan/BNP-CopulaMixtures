temp <- read.csv("../../Revision/Data/Pollution/polution3A2017rank.csv")
colnames(temp)

polution_rank_month4 <- temp[,c("PM25.4","NO2.4","SO2.4" )]
polution_rank_month6 <- temp[,c("PM25.6","NO2.6","SO2.6" )]
polution_rank_month12 <- temp[,c("PM25.1","NO2.1","SO2.1" )]
colnames(polution_rank_month4)
colnames(polution_rank_month6)

saveRDS(polution_rank_month4, file="../Data/Pollution/pollution_rank_month4.rds")
saveRDS(polution_rank_month6, file="../Data/Pollution/pollution_rank_month6.rds")

nrow(polution_rank_month6)
cor(polution_rank_month12, method="kendall")
