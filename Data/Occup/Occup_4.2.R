occup <- read.csv("/home/panruyi/scratch/BNP_Copula/Data/Occup/occupDIFrank.csv")
occup<-read.csv("../Data/Occup/occupDIFrank.csv")
occup_single_day <- occup %>% filter(DAY==5) %>% 
  filter(CO2!=1,HR!=1,CO2!=0, HR!=0)
data <- as.matrix(occup_single_day[,c(2,3)])
sample_size <- nrow(data)
