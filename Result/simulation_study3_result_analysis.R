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

#high-dimension, 4 variables


folder_name <- "Clayton_Hdim/"
files <- list.files(folder_name)
# range <- seq(100*50, 170*50, by=5)
range <- seq(100*50, 200*50, by=5)
est_pars_all1 <- list()
est_pars_all2 <- list()
est_pars_all3 <- list()
est_pars_all4 <- list()
est_tau_all <- list()
est_tau_lower <- list()
est_tau_upper <- list()
for (job_num in c(1:3)) {
  for (file in files) {
    if (grepl(paste0("job_num=", job_num, "sample"), file) ) {
      temp_res <- readRDS(paste0(folder_name,file))
      LPML[job_num] <- temp_res$LPML
      temp_num <- NULL
      for (index in range) {
        temp_num <- c(temp_num,
                      length(unique(temp_res$res_mcmc$all_theta[index,])))
      }
      
      num_comp[job_num] <- Mode(temp_num)
      num_comp[job_num] <- length(unique(temp_res$res_theta_post$all_pos_theta[1,]))
      
      type="C"
      hname <-  paste0("Figures/", folder_name, "COMP_Type=", type,"Setting",job_num, ".pdf")
      pdf(hname)
      
      # Create the plot again (it needs to be recreated)
      hist(temp_num, prob=T,
           xlab = "m", 
           ylab = "", 
           cex.lab=2,
           font.lab=2,
           main = paste0("Setting ",job_num),
           xaxt = 'n',
           yaxt = 'n',
           cex.main=2, breaks = seq(0, max(temp_num)))
      axis(1, font.axis = 2, cex.axis = 2)
      axis(2, font.axis = 2, cex.axis = 2)
      # Close the file
      dev.off()
      
      
      par_name <-  paste0("Figures/", folder_name, "Pars_Type=", type, "Setting",job_num, ".pdf")
      pdf(par_name)
      
      # Create the plot again (it needs to be recreated)
      thetas <- temp_res$res_mcmc$all_theta[range,]
      hist(temp_res$res_mcmc$all_theta[range,], prob=T,
           xlab = expression(theta), 
           ylab = "", 
           cex.lab=3,
           font.lab=2,
           main = paste0("Setting ",job_num),
           xaxt = 'n',
           yaxt = 'n',
           cex.main=2, xlim=c(min(thetas), max(thetas)))
      axis(1, font.axis = 2, cex.axis = 2)
      axis(2, font.axis = 2, cex.axis = 2)
      # Close the file
      dev.off()
      est_pars_all1[[job_num]] <- temp_res$est_pars_weight1
      est_pars_all2[[job_num]] <- temp_res$est_pars_weight2
      est_pars_all3[[job_num]] <- temp_res$est_pars_weight3
      est_pars_all4[[job_num]] <- temp_res$est_pars_weight4
      temp_taus <- getTau(temp_res$res_den$theta_null, type)
      est_tau_lower[[job_num]]<- quantile(temp_taus, 0.025)
      est_tau_upper[[job_num]]<- quantile(temp_taus, 0.975)
      est_tau_all[[job_num]]<-mean(temp_taus)
      #est_density
      p <- expand.grid(u = seq(0,1,length.out=100), v = seq(0,1,length.out=100)) %>%
        mutate(Z = temp_res$res_den$den) %>%
        ggplot(aes(u, v, z = Z)) +
        geom_contour() +
        geom_contour_filled(breaks = seq(0,5,by=0.1) ) +
        scale_fill_viridis_d(drop = FALSE) +
        ggtitle(paste0("Setting ",job_num))+
        theme_bw() +
        theme(plot.title = element_text(hjust = 0.5, size=40,face="bold"),aspect.ratio = 1) +
        theme(axis.text=element_text(size=22,face="bold"),
              axis.title=element_text(size=28,face="bold")) +
        theme(legend.position="none")
      fname <- paste0("Figures/",folder_name, "est_dense=", type, "Setting",job_num, ".pdf")
      ggsave(fname, plot = p, width = 7, height = 8.5, units = "in")
    }
  }
}


setting=3; digit=3
df1 <- as.data.frame(round(est_pars_all1[[setting]],digit))
row.names(all_df) <- NULL
kable(df1,
      format = "latex",
      booktabs = TRUE)

#estimated kendall's tau
rbind(round(unlist(est_tau_lower),3),
      round(unlist(est_tau_all),3),
      round(unlist(est_tau_upper),3))

#true kendall's tau
tau_true <- round(c(sum(getTau(c(1, 5, 15),"C")*c(0.2,0.3,0.5)),
                    sum(getTau(c(2, 5,10),"C")*c(0.2,0.3,0.5)),
                    sum(getTau(c(2, 7, 15),"C")*c(0.2,0.3,0.5))),4)
tau_true 

df <- as.data.frame(round(est_pars_all2[[5]],4))

latex_table <- xtable(df, digits = c(0, rep(4,4)),caption = "Sample Matrix", label = "tab:sample_matrix")

# Print the LaTeX table code
print(latex_table)


