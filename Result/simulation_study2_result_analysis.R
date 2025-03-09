#Result analysis and figure generation
library(tidyverse)
library(xtable)
library(tibble)
library(dplyr)
library(knitr)
library(copula)
source("res_helper.R")

source("../../Code/density.R")
Mode <- function(x) {
  ux <- unique(x)
  ux[which.max(tabulate(match(x, ux)))]
}
folder_name <- "Gaussian_mix/"  #change to  "Gaussian_mix/" for Gaussian copula with two components
files <- list.files(folder_name)
range <- seq(100*50, 170*50, by=5)
type_num <- list("A"=9, "F"=5, "C"=3, "J"=7, "G"=1)
LPML <- list()
est_tau_all <- list()
est_tau_lower <- list()
est_tau_upper <- list()

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
    }
  }
}


#
kable(round(do.call(cbind, LPML),2),
      format = "latex",
      booktabs = TRUE)

#estimated kendall
kable(t(rbind(round(unlist(est_tau_all),3),
              round(unlist(est_tau_lower),3),
              round(unlist(est_tau_upper),3))),
      format = "latex",
      booktabs = TRUE)

df <- round(as.data.frame(temp_res$cc$NID$pars_weight_est),3)
kable(df,
      format = "latex",
      booktabs = TRUE)



#real heatmap
#for Gaussain copula with single component
data <- as.matrix(expand.grid(u = seq(0,1,length.out=100), v = seq(0,1,length.out=100)))
nC1 <- normalCopula(0.7, dim = 2, dispstr = "ex")

dense <- dCopula(data, nC1)
p <- expand.grid(u = seq(0,1,length.out=100), v = seq(0,1,length.out=100)) %>%
  mutate(Z = temp_res$res_den$den) %>%
  ggplot(aes(u, v, z = dense)) +
  geom_contour() +
  geom_contour_filled(breaks = seq(0,5,by=0.1) ) +
  scale_fill_viridis_d(drop = FALSE) +
  ggtitle(paste0("True density"))+
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5, size=40,face="bold"),aspect.ratio = 1) +
  theme(axis.text=element_text(size=22,face="bold"),
        axis.title=element_text(size=28,face="bold")) +
  theme(legend.position="none")
p
fname <- paste0("../Result_new/Figures/Gaussian_single/","dense_True", ".pdf")
ggsave(fname, plot = p, width = 7, height = 8.5, units = "in")

#for Gaussain copula with mixture of two components
nC1 <- normalCopula(0.7, dim = 2, dispstr = "ex")
nC2 <- normalCopula(-0.7, dim = 2, dispstr = "ex")
dense <- (dCopula(data, nC1)+dCopula(data, nC2))/2
p <- expand.grid(u = seq(0,1,length.out=100), v = seq(0,1,length.out=100)) %>%
  mutate(Z = temp_res$res_den$den) %>%
  ggplot(aes(u, v, z = dense)) +
  geom_contour() +
  geom_contour_filled(breaks = seq(0,5,by=0.1) ) +
  scale_fill_viridis_d(drop = FALSE) +
  ggtitle(paste0("True density"))+
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5, size=40,face="bold"),aspect.ratio = 1) +
  theme(axis.text=element_text(size=22,face="bold"),
        axis.title=element_text(size=28,face="bold")) +
  theme(legend.position="none")
p
fname <- paste0("../Result_new/Figures/Gaussian_mix/","dense_True", ".pdf")
ggsave(fname, plot = p, width = 7, height = 8.5, units = "in")
