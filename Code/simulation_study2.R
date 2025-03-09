source("density.R")
source("sampling_helper.R")
source("post_sampling_helper.R")
library(foreach)
library(doParallel)
library(copula)
args=(commandArgs(TRUE))
job_name=args[1]
job_num=as.numeric(args[2])
path=args[3]
#change file to "Gaussian_sim_data/single.rds" to run for simulated gaussian copula with single component
Gaussian_data <- readRDS("/home/panruyi/scratch/BNP_Copula/Data/Gaussian_sim_data/mix.rds")
data <- as.matrix(Gaussian_data)
sample_size <- nrow(data)
settings <- c("A", "F", "C", "J", "G")
type=settings[job_num]

set.seed(123)  # for reproducibility
#hyper-parameters
if (type=="F") {
  pars_G0 <- c(0, 4)
} else if (type == "A") {
  pars_G0 <- c(-1, 1)
} else {
  pars_G0 <- c(4, 1)
}
pars_a0 <- c(1, 20)
pars_b0 <- c(1, 20)


#tuning parameters for adaptive procedure
kappa_theta <- 1
all_kappa_theta <- 1
rate_theta <- 1.1
kappa_a <- 1
all_kappa_a <- 1
rate_a <- 1.1
kappa_b <- 1
all_kappa_b <- 1
rate_b <- 1.1

#number of mcmc iterations
burn_in <- 100
batch_num <- 170
batch_size <- 50
mcmc_iter <- batch_size*batch_num

#sample size, number of theta
n <- nrow(data)

#number of auxilary thetas
r <- 10

#initialize the value for each paramters, theta, a, b
current_theta <- rG0(n, pars_G0, type)
a <- 0.5
b <- 1

#Record the mcmc samples, theta, a, b
all_theta <- current_theta
all_a <- a
all_b <- b


#accpetance rate 
acc_a <- 0
acc_b <- 0
acc_theta <- 0 

#record acceptance rate
all_acc_a <- 0
all_acc_b <- 0
all_acc_theta <- 0 


for (k in 1:batch_num) {
  print(paste0("iteration: ", k))
  print(paste0("acceptance rate for a: ", acc_a))
  print(paste0("acceptance rate for b: ", acc_b))
  print(paste0("acceptance rate for theta: ", acc_theta))
  print(paste0("number of unique theta: ", length(unique(current_theta))))
  count_a <- 0
  count_b <- 0
  count1_theta <- 0  # record  number of proposal
  count2_theta <- 0  #record number of time accepting proposal
  for (t in 1:batch_size) {
    #step 1 and step 2
    
    # draw new theta_i
    for (i in 1:n) {
      # sample auxiliary values from G0
      theta_auxi <- rG0(r, pars_G0, type)
      #Get the frequency of unique theta_-i
      freq <- table(current_theta[-i])
      #number of distinct values
      mi <- length(freq)
      theta_temp <- c(as.numeric(names(freq)), theta_auxi)
      #compute the likelihood for each unique theta and auxiliary theta
       
      dens <- sapply(theta_temp, function(x) dcopula_new(data[i,], x, type)) 
      #get the unnormalized posterior probability for each theta
      pos_theta <- c(as.vector(freq)-a, rep((b+a*mi)/r,r))*dens
      #normalize the probability
      pos_theta_norm <- pos_theta/sum(pos_theta)
      #draw a new theta_i, and update the current theta
      current_theta[i] <- sample(theta_temp,size=1, prob=pos_theta_norm)
    }
    
    
    #step 3
    #cluster configuration, move a bit within each cluster
    
    #compute the unique values of theta
    old_theta <- current_theta
    theta_unique_old <- unique(current_theta)
    m <- length(theta_unique_old) #number of clusters
    count1_theta <- count1_theta + m # increase number of proposal
    for (j in 1:m) {
      theta_old_j <-theta_unique_old[j]
      cc_index <- which(old_theta==theta_old_j)
      data_cc <- data[cc_index, , drop=FALSE]
      theta_prop_j <- rtheta_cc(theta_old_j, kappa_theta, type)
      acc_prob_theta <- min(1, exp(log_f_cc(data_cc, theta_prop_j, pars_G0, type)-
                                     log_f_cc(data_cc, theta_old_j, pars_G0, type)))
      if (runif(1) <= acc_prob_theta ) {
        current_theta[cc_index] <- theta_prop_j
        count2_theta <- count2_theta + 1 # increase number of accepted proposal
      }
      
    }
    
    
    
    all_theta <- rbind(all_theta, current_theta)
    
    #Step 4
    #sample a
    a_prop <- runif(1, max(0, a-kappa_a), min(1, a+kappa_a)) #propose a new a
    cluster_size <- as.vector(table(current_theta)) #get the size for each cluster
    if (a_prop == 0) {
      acc_prob_a <- 0
    } else {
      acc_prob_a <- min(1, exp(log_f_a(a_prop, b, cluster_size, pars_a0, pars_b0)-
                                 log_f_a(a, b, cluster_size, pars_a0, pars_b0)))
    }
    
    if (runif(1) <= acc_prob_a) {
      a <- a_prop
      all_a <- c(all_a, a)
      count_a <- count_a + 1
    } else {
      all_a <- c(all_a, a)
    }
    
    #step 5
    #sample b
    b_prop <- rgamma(1, kappa_b, kappa_b/b)
    acc_prob_b <- min(1, exp(log_f_b(a, b_prop, n, length(cluster_size), pars_b0)+
                               dgamma(b, kappa_b, kappa_b/b_prop, log=T)-
                               log_f_b(a, b, n, length(cluster_size), pars_b0)-
                               dgamma(b_prop, kappa_b, kappa_b/b, log=T)))
    if (runif(1) <= acc_prob_b) {
      b <- b_prop
      all_b <- c(all_b, b)
      count_b <- count_b + 1
    } else {
      all_b <- c(all_b, b)
    }
    
  } 
  
  #Compute the accepetance rate, and update the kappa values
  
  #acceptance rate for theta
  acc_theta <- count2_theta/ count1_theta
  all_acc_theta <- c(all_acc_theta, acc_theta)
  if (acc_theta > 0.45) {
    kappa_theta <- kappa_theta*rate_theta^(sqrt(k))
  } else if (acc_theta < 0.35) {
    kappa_theta <- kappa_theta*rate_theta^(-sqrt(k))
  }
  
  #acceptance rate for a 
  acc_a <- count_a / batch_size
  all_acc_a <- c(all_acc_a, acc_a)
  if (acc_a > 0.45) {
    kappa_a <- kappa_a*rate_a^(-sqrt(k))
  } else if (acc_a < 0.35) {
    kappa_a <- kappa_a*rate_a^(sqrt(k))
  }
  
  #acceptance rate for b
  acc_b <- count_b / batch_size
  all_acc_b <- c(all_acc_b, acc_b)
  if (acc_b > 0.45) {
    kappa_b <- kappa_b*rate_b^(-sqrt(k))
  } else if (acc_b < 0.35) {
    kappa_b <- kappa_b*rate_b^(sqrt(k))
  }
  
  all_kappa_a <- c(all_kappa_a, kappa_a)
  all_kappa_b <- c(all_kappa_b, kappa_b)
  all_kappa_theta <- c(all_kappa_theta, kappa_theta)
}



#
print("mcmc:finish")



#store mcmc results
res_mcmc <- list(data=data, all_theta=all_theta, all_a=all_a, all_b=all_b,
                 all_acc_theta=all_acc_theta,all_acc_a=all_acc_a, all_acc_b=all_acc_b,
                 all_kappa_theta=all_kappa_theta, all_kappa_a=all_kappa_a, all_kappa_b=all_kappa_b)



##LPML
thin <- 5
pos_theta <- res_mcmc$all_theta[seq(burn_in*batch_size,batch_num*batch_size,by=thin),]
GlobalFunctions = ls(globalenv())
registerDoParallel(25)

CPO <- foreach(k=1:n, .export=GlobalFunctions, 
               .packages = c("copula"), .combine="rbind")%dopar%{
                 L <- nrow(pos_theta)
                 dense_inv <- lapply(1:L, function(l) {
                   1 / dcopula_new(unlist(data[k,]), pos_theta[l,k],type)})
                 
                 mean(do.call(rbind,dense_inv))
               }

LPML <- sum(log(1/CPO))


##post sampling
pos_theta <- res_mcmc$all_theta[seq(burn_in*batch_size,batch_num*batch_size,by=thin),]
pos_a <- res_mcmc$all_a[seq(burn_in*batch_size,batch_num*batch_size,by=thin)]
pos_b <- res_mcmc$all_b[seq(burn_in*batch_size,batch_num*batch_size,by=thin)]
L <- nrow(pos_theta)
theta_null <- vector()
for(l in 1:L) {
  freq <- table(pos_theta[l,])
  #number of distinct values
  # print(freq)
  m <- length(freq)
  temp_a <- pos_a[l]
  temp_b <- pos_b[l]
  theta_auxi <- rG0(1, pars_G0, type)
  theta_temp <- c(as.numeric(names(freq)), theta_auxi)
  temp_prob <- c((freq-temp_a)/(temp_b+n),
                 (temp_b+temp_a*m)/(temp_b+n))
  # print(theta_temp)
  # print(temp_prob)
  theta_null[l] <- sample(theta_temp, size=1, prob=temp_prob)
}

print("finish post-sampling")

data_points <- expand.grid(seq(0,1,length.out=100),seq(0,1,length.out=100))
GlobalFunctions = ls(globalenv())
registerDoParallel(25)
den_pos <- foreach(k=1:100,.export=GlobalFunctions, .packages = c("copula"), .combine="rbind")%dopar%{
  L <- length(theta_null)
  range <- (1+(k-1)*100):(k*100)
  dense_est <- lapply(range, function(j){
    dense <- lapply(1:L, function(l) {
      dcopula_new(unlist(data_points[j,]), theta_null[l],type)})
    mean(do.call(rbind,dense))
  })
  
  do.call(rbind,dense_est)
}

registerDoSEQ()


res_den<- list(points=data_points, den=den_pos, theta_null=theta_null)


print("Approach 1: get estimated joint density")



#Approach one, within candiate 
#select the a cluster to start, which determine the number of cluster
#and initial values for each theta
kappa_theta <- res_mcmc$all_kappa_theta[batch_num]
range <- seq(burn_in*batch_size,batch_num*batch_size,by=thin) 
pos_theta <- round(res_mcmc$all_theta[range,],2)
batch_size_post <- 50
batch_burn_in <- 20
batch_num_post <- 60

cc_results <- list()
methods <- c("binder", "omARI", "VI", "NVI", "ID", "NID")
for (m in methods) {
  # Call findcc_salso() with the current method
  cc_obj <- findcc_salso(pos_theta, method = m)
  
  # Sample from the posterior given the current cc object.
  res_theta_post_tmp <- post_sampling(data, cc_obj$cc, cc_obj$pars_unique, 
                                      cc_obj$num_pars, kappa_theta,
                                      batch_num_post, batch_size_post, batch_burn_in,
                                      pars_G0, type, rate_theta)
  
  # Estimate parameter weights.
  est_pars_weight_tmp <- est_pars_weight(res_theta_post_tmp, res_mcmc, range)
  
  # Store all outputs for the current method into the list.
  cc_results[[m]] <- list(
    cc = cc_obj,
    res_theta_post = res_theta_post_tmp,
    pars_weight_est = est_pars_weight_tmp
  )
}


res <- list(res_mcmc=res_mcmc,
            LPML=LPML,
            res_den=res_den,
            cc=cc_results)

filename <- paste0("job_name=", job_name,"job_num=", job_num,
                   "sample_size=", sample_size, 
                   "res_type=",type, ".rds")
saveRDS(res, paste0(path,"/",filename))
