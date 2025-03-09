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

#True density
for (type in c("A",  "C", "F","J","G")) {
  
  data_points <- expand.grid(seq(0,1,length.out=100),seq(0,1,length.out=100))
  if (type == "A") {
    temp_dense <- lapply(1:10000, function(j){sum(sapply(c(-0.8,0.8), function(x)
      dcopula_new(unlist(data_points[j,]), x,"A")))/2})
  } else if (type == "C") {
    temp_dense <- lapply(1:10000, function(j){sum(sapply(c(-0.5,10), function(x)
      dcopula_new(unlist(data_points[j,]), x,"C")))/2})
  } else if (type == "F") {
    temp_dense <- lapply(1:10000, function(j){sum(sapply(c(-5,5), function(x)
      dcopula_new(unlist(data_points[j,]), x,"F")))/2})
  } else if (type == "J") {
    temp_dense <- lapply(1:10000, function(j){sum(sapply(c(2,10), function(x)
      dcopula_new(unlist(data_points[j,]), x,"J")))/2})
  } else if (type == "G") {
    temp_dense <- lapply(1:10000, function(j){sum(sapply(c(5,10), function(x)
      dcopula_new(unlist(data_points[j,]), x,"G")))/2})
  }
  
  temp_dense <- do.call(rbind, temp_dense)
  
  
  p <- expand.grid(u = seq(0,1,length.out=100), v = seq(0,1,length.out=100)) %>%
    mutate(Z =temp_dense) %>%
    ggplot(aes(u, v, z = Z)) +
    geom_contour() +
    geom_contour_filled(breaks = seq(0,5,by=0.1) ) +
    scale_fill_viridis_d(drop = FALSE) +
    ggtitle(paste0("True:", type))+
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5, size=40,face="bold"),aspect.ratio = 1) +
    theme(axis.text=element_text(size=22,face="bold"),
          axis.title=element_text(size=28,face="bold")) +
    theme(legend.position="none")
  
  fname <- paste0("Figures/realdense_Type=", type, ".pdf")
  ggsave(fname, plot = p,  width = 7, height = 8.5, units = "in")
  
}




#simulation results 
settings <- expand.grid(c("A","F","C","J","G"), c("A","F","C","J","G"))
LPML <- vector()
num_comp <- vector()
folder_name <- "Simulation500/" # change to "Simulation200/" for sample size 200.
files <- list.files(folder_name)
range <- seq(100*50, 170*50, by=5)
for (job_num in c(1:25)) {
  for (file in files) {
    if (grepl(paste0("job_num=", job_num, "sample"), file) ) {
      temp_res <- readRDS(paste0(folder_name,file))
      LPML[job_num] <- temp_res$LPML
      temp_num <- NULL
      for (index in range) {
        temp_num <- c(temp_num,
                      length(unique(temp_res$res_mcmc$all_theta[index,])))
      }
      
      # num_comp[job_num] <- Mode(temp_num)
      # num_comp[job_num] <- length(unique(temp_res$res_theta_post$all_pos_theta[1,]))
      
      stype <- settings[job_num, 1]
      type <- settings[job_num, 2]
      hname <-  paste0("Figures/", folder_name,"COMP_realType=", stype, "Type=", type, ".pdf")
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
      
      
      par_name <-  paste0("Figures/", folder_name,"par_realType=", stype, "Type=", type, ".pdf")
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
      
      fname <- paste0("Figures/",folder_name,"Large_realType=", stype, "Type=", type, ".pdf")
      ggsave(fname, plot = p, width = 7, height = 8.5, units = "in")
      
    }
  }
}

LPML_large <- rbind(LPML[seq(1,25,by=5)],
                    LPML[seq(2,25,by=5)],
                    LPML[seq(3,25,by=5)],
                    LPML[seq(4,25,by=5)],
                    LPML[seq(5,25,by=5)])
colnames(LPML_large) <- c("A",  "F", "C","J", "G")
row.names(LPML_large) <- c("A", "F", "C",  "J", "G")
round(LPML_large,2)



#fetch the components
folder_name <- "Simulation500/" 
files <- list.files(folder_name)
est_pars_all <- list()
num_comp <- NULL
median_comp <- NULL
for (job_num in c(1,7,13,19,25)) {
  for (file in files) {
    if (grepl(paste0("job_num=", job_num, "sample"), file) ) {
      temp_res<- readRDS(paste0(folder_name,file))
      #summarize the number of component and estiamted parameters
      temp_num <- NULL
      for (index in range) {
        temp_num <- c(temp_num,
                      length(unique(temp_res$res_mcmc$all_theta[index,])))
      }
      
      num_comp <- c(num_comp, Mode(temp_num))
      median_comp <- c(median_comp,median(temp_num))
      stype <- settings[job_num, 1]
      type <- settings[job_num, 2]
      
      # Sort the matrix by the specified row
      est_pars_all[[stype]] <- temp_res$est_pars_weight1
    }
  }
}

type=5; digit=3

df <- as.data.frame(round(est_pars_all[[type]],digit))
kable(df,
      format = "latex",
      booktabs = TRUE,
      caption = "Mixture Parameter Estimates with Variable Number of Components")
# Print the LaTeX table code
print(latex_table)





#compute kendall's tau
types <- c("A", "F", "C", "J", "G")
tau_true <- round(c(mean(getTau(c(0.8,-0.8),"A")),
                    mean(getTau(c(5,-5),"F")),
                    mean(getTau(c(10,-0.5),"C")),
                    mean(getTau(c(10,2),"J")),
                    mean(getTau(c(10,5),"G"))),4)


settings <- expand.grid(c("A","F","C","J","G"), c("A","F","C","J","G"))
folder_name <- "Simulation500/"
files <- list.files(folder_name)
# range <- seq(100*50, 170*50, by=5)
range <- seq(100*50, 170*50, by=5)
tau_emp <- NULL
for (job_num in c(1,7,13,19,25)) {
  for (file in files) {
    if (grepl(paste0("job_num=", job_num, "sample"), file) ) {
      temp_res<- readRDS(paste0(folder_name,file))
      #summarize the number of component and estiamted parameters
      data <- temp_res$res_mcmc$data
      tau_emp <- c(tau_emp,cor(data,method="kendall")[1,2])
      plot(data[,1],data[,2])
    }
  }
}

settings <- expand.grid(c("A","F","C","J","G"), c("A","F","C","J","G"))
tau_est <- vector()
tau_est_lower <- vector()
tau_est_upper <- vector()
folder_name <- "Simulation500/"
files <- list.files(folder_name)
# range <- seq(100*50, 170*50, by=5)
range <- seq(100*50, 170*50, by=5)
for (job_num in c(1:25)) {
  for (file in files) {
    if (grepl(paste0("job_num=", job_num, "sample"), file) ) {
      temp_res <- readRDS(paste0(folder_name,file))
      type <- settings[job_num, 2]
      temp_taus <- getTau(temp_res$res_den$theta_null, type)
      tau_est_lower[job_num] <- quantile(temp_taus, 0.025)
      tau_est_upper[job_num] <- quantile(temp_taus, 0.975)
      tau_est[job_num] <- mean(temp_taus)
      # num_comp[job_num] <- Mode(temp_num)
      # num_comp[job_num] <- length(unique(temp_res$res_theta_post$all_pos_theta[1,])
    }
  }
}

tau_est_all<- rbind(tau_est[seq(1,25,by=5)],
                    tau_est[seq(2,25,by=5)],
                    tau_est[seq(3,25,by=5)],
                    tau_est[seq(4,25,by=5)],
                    tau_est[seq(5,25,by=5)])

tau_est_lower_all<- rbind(tau_est_lower[seq(1,25,by=5)],
                          tau_est_lower[seq(2,25,by=5)],
                          tau_est_lower[seq(3,25,by=5)],
                          tau_est_lower[seq(4,25,by=5)],
                          tau_est_lower[seq(5,25,by=5)])

tau_est_upper_all<- rbind(tau_est_upper[seq(1,25,by=5)],
                          tau_est_upper[seq(2,25,by=5)],
                          tau_est_upper[seq(3,25,by=5)],
                          tau_est_upper[seq(4,25,by=5)],
                          tau_est_upper[seq(5,25,by=5)])

tau_all <- cbind(tau_true, tau_emp, tau_est_all)
colnames(tau_all) <- c("True","Empirical", "A", "F", "C", "J", "G")
row.names(tau_all) <- c("A", "F", "C", "J", "G")

df <- as.data.frame(round(tau_all,2))

kable(df, format = "latex",
      booktabs = TRUE)







