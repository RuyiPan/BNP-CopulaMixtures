library(tidyverse)
library(xtable)
library(tibble)
library(dplyr)
library(knitr)
library(HAC)
source("res_helper.R")

folder_name <- "real_wine/"
files <- list.files(folder_name)
range <- seq(100*50, 170*50, by=5)
type_num <- list("A"=9, "F"=5, "C"=3, "J"=7, "G"=1)
LPML <- list()
est_pars_all1 <- list()
est_pars_all2 <- list()
est_pars_all3 <- list()
est_pars_all4 <- list()
est_tau_all <- list()
est_tau_lower <- list()
est_tau_upper <- list()
# c("A", "F", "C", "J", "G")
for (type in c("A", "F", "C", "J", "G")) {
  print(type)
  for (file in files) {
    if (grepl(paste0("type=", type), file) ) {
      print(type)
      temp_res <- readRDS(paste0(folder_name,file))
      temp_taus <- getTau(temp_res$res_den$theta_null, type)
      est_tau_lower[[type]]<- quantile(temp_taus, 0.025)
      est_tau_upper[[type]]<- quantile(temp_taus, 0.975)
      est_tau_all[[type]]<-mean(temp_taus)
      LPML[[type]] <- temp_res$LPML
      temp_num <- NULL
      for (index in range) {
        temp_num <- c(temp_num,
                      length(unique(temp_res$res_mcmc$all_theta[index,])))
      }
      
      # num_comp[job_num] <- Mode(temp_num)
      # num_comp[job_num] <- length(unique(temp_res$res_theta_post$all_pos_theta[1,]))
      
      hname <-  paste0("Figures/", folder_name,"COMP_Type=", type, ".pdf")
      pdf(hname)
      
      # Create the plot again (it needs to be recreated)
      hist(temp_num, prob=T,
           xlab = "m", 
           ylab = "", 
           cex.lab=2,
           font.lab=2,
           main = paste0("Model:", type),
           xaxt = 'n',
           yaxt = 'n',
           cex.main=2, breaks = seq(0, max(temp_num)))
      axis(1, font.axis = 2, cex.axis = 2)
      axis(2, font.axis = 2, cex.axis = 2)
      # Close the file
      dev.off()
      
      
      par_name <-  paste0("Figures/", folder_name,"par_Type=", type, ".pdf")
      pdf(par_name)
      
      # Create the plot again (it needs to be recreated)
      thetas <- temp_res$res_mcmc$all_theta[range,]
      hist(temp_res$res_mcmc$all_theta[range,], prob=T,
           xlab = expression(theta), 
           ylab = "", 
           cex.lab=3,
           font.lab=2,
           main = paste0("Model:", type),
           xaxt = 'n',
           yaxt = 'n',
           cex.main=2, xlim=c(min(thetas), max(thetas)))
      axis(1, font.axis = 2, cex.axis = 2)
      axis(2, font.axis = 2, cex.axis = 2)
      # Close the file
      dev.off()
      
      p <- expand.grid(u = seq(0,1,length.out=100), v = seq(0,1,length.out=100)) %>%
        mutate(Z = temp_res$res_den$den) %>%
        ggplot(aes(u, v, z = Z)) +
        geom_contour() +
        geom_contour_filled(breaks = seq(0,5,by=0.1) ) +
        scale_fill_viridis_d(drop = FALSE) +
        ggtitle(paste0("Model:", type))+
        theme_bw() +
        theme(plot.title = element_text(hjust = 0.5, size=40,face="bold"),aspect.ratio = 1) +
        theme(axis.text=element_text(size=22,face="bold"),
              axis.title=element_text(size=28,face="bold")) +
        theme(legend.position="none")
      
      fname <- paste0("Figures/",folder_name,"dense_Type=", type, ".pdf")
      ggsave(fname, plot = p, width = 7, height = 8.5, units = "in")
      est_pars_all1[[type]] <-   temp_res$est_pars_weight1
      est_pars_all2[[type]] <- temp_res$est_pars_weight2
      est_pars_all3[[type]] <- temp_res$est_pars_weight3
      est_pars_all4[[type]] <- temp_res$est_pars_weight4
    }
  }
}

#LPML
round(unlist(LPML),3)


#estimated parameters
#Frank==2, Gumbel==5
type=5; digit=3
df3 <- as.data.frame(round(est_pars_all3[[type]],digit))
kable(df3,
      format = "latex",
      booktabs = TRUE)


#empprical kendall
cor(temp_res$res_mcmc$data,
    method="kendall")
#estimated kendall
rbind(round(unlist(est_tau_lower),3),
      round(unlist(est_tau_all),3),
      round(unlist(est_tau_upper),3))
