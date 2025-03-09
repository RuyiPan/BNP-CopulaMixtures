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
sample_size=as.numeric(args[4]) #sample_size=500

job_num=3 #change the setting, there are 3 settings in the paper
dim=4
stype="C"
type="C"


set.seed(123)  # for reproducibility
data <- NULL
if (job_num == 1) {
  clayton_copula1 <- claytonCopula(dim = dim, param = 1)
  clayton_copula2 <- claytonCopula(dim = dim, param = 5)
  clayton_copula3 <- claytonCopula(dim = dim, param = 15)
} else if (job_num == 2) {
  clayton_copula1 <- claytonCopula(dim = dim, param = 2)
  clayton_copula2 <- claytonCopula(dim = dim, param = 5)
  clayton_copula3 <- claytonCopula(dim = dim, param = 10)
} else if (job_num == 3) {
  clayton_copula1 <- claytonCopula(dim = dim, param = 2)
  clayton_copula2 <- claytonCopula(dim = dim, param = 7)
  clayton_copula3 <- claytonCopula(dim = dim, param = 15)
} else if (job_num == 4) {
  clayton_copula1 <- claytonCopula(dim = dim, param = 2)
  clayton_copula2 <- claytonCopula(dim = dim, param = 5)
  clayton_copula3 <- claytonCopula(dim = dim, param = 10)
} else if (job_num == 5) {
  clayton_copula1 <- claytonCopula(dim = dim, param = 2)
  clayton_copula2 <- claytonCopula(dim = dim, param = 7)
  clayton_copula3 <- claytonCopula(dim = dim, param = 15)
}


id <- vector()
for (i in 1:sample_size) {
  if (job_num %in% c(1,2,3)) {
    id[i] <- sample(c(1,2,3), size=1, replace = T, prob=c(0.2,0.3,0.5))
  } else if (job_num == 4) {
    id[i] <- sample(c(1,2,3), size=1, replace = T, prob=c(0.5,0.3,0.2))
  } else if (job_num == 5) {
    id[i] <- sample(c(1,2,3), size=1, replace = T, prob=c(0.2,0.5,0.3))
  }
  
  if (id[i] == 1) {
    temp_data <- rCopula(1, clayton_copula1)
  } else if (id[i] == 2) {
    temp_data <- rCopula(1, clayton_copula2)
  } else if (id[i] == 3) {
    temp_data <- rCopula(1, clayton_copula3)
  }
  data <- rbind(data, temp_data)
}


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
batch_num <- 200
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
      # dens <- sapply(theta_temp, function(x) dcopula(data[i,], x, type)) 
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
    # a_prop <- rgamma(1, kappa_a, kappa_a/a) #propose a new a
    # cluster_size <- as.vector(table(current_theta)) #get the size for each cluster
    # if (a_prop == 0) {
    #   acc_prob_a <- 0
    # } else {
    #   acc_prob_a <- min(1, exp(log_f_a(a_prop, b, cluster_size, pars_a0, pars_b0)+
    #                              dgamma(a, kappa_a, kappa_a/a_prop,log=T)-
    #                              log_f_a(a, b, cluster_size, pars_a0, pars_b0)-
    #                              dgamma(a_prop, kappa_a, kappa_a/a, log=T)))
    # }
    
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
    # acc_prob_b <- min(1, f_b(a, b_prop, n, length(cluster_size), pars_b0)*
    #                     dgamma(b, kappa_b, kappa_b/b_prop)/
    #                     f_b(a, b, n, length(cluster_size), pars_b0)/
    #                     dgamma(b_prop, kappa_b, kappa_b/b))
    
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

CPO <- foreach(k=1:(sample_size), .export=GlobalFunctions, 
               .packages = c("copula"), .combine="rbind")%dopar%{
                 L <- nrow(pos_theta)
                 dense_inv <- lapply(1:L, function(l) {
                   1 / dcopula_new(unlist(data[k,]), pos_theta[l,k],type)})
                 
                 mean(do.call(rbind,dense_inv))
               }

LPML <- sum(log(1/CPO))


###post sampling
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




#select the a cluster to start, which determine the number of cluster
#and initial values for each theta
#Approach one, within candiate 
#select the a cluster to start, which determine the number of cluster
#and initial values for each theta
batch_size_post <- 50
batch_num_post <- 40
kappa_theta <- res_mcmc$all_kappa_theta[batch_num]
range <- seq(burn_in*batch_size,batch_num*batch_size,by=thin) 
pos_theta <- res_mcmc$all_theta[range,]



len_theta <- nrow(pos_theta)
temp_num <- NULL
for (id in 1:len_theta) {
  temp_num <- c(temp_num,
                length(unique(pos_theta[id,])))
}

comp_mode <- Mode(temp_num)
if (comp_mode == 1) {
  num_candidate <- c(1,2,3,4,5)
} else if (comp_mode == 2) {
  num_candidate <- c(1,2,3,4,5)
} else {
  num_candidate <- c(-2,-1,0,1,2)+comp_mode
}
candidate_index <- which(temp_num %in% num_candidate)
pos_theta_candidate <- pos_theta[candidate_index,]


cc_old <- findcc_old(pos_theta)

res_theta_post1 <- post_sampling(data, cc_old$cc, cc_old$pars_unique, 
                                 cc_old$num_pars, kappa_theta,
                                 batch_num_post, batch_size_post,
                                 pars_G0, type, rate_theta)

est_pars_weight1 <- est_pars_weight(res_theta_post1, res_mcmc, range)


cc_old_5mod <- findcc_old(pos_theta_candidate)
res_theta_post2 <- post_sampling(data, cc_old_5mod$cc, cc_old_5mod$pars_unique, 
                                 cc_old_5mod$num_pars, kappa_theta,
                                 batch_num_post, batch_size_post,
                                 pars_G0, type, rate_theta)

est_pars_weight2 <- est_pars_weight(res_theta_post2, res_mcmc, range)


salso_cc <- findcc_salso(pos_theta)  

res_theta_post3 <- post_sampling(data, salso_cc$cc, salso_cc$pars_unique, 
                                 salso_cc$num_pars, kappa_theta,
                                 batch_num_post, batch_size_post,
                                 pars_G0, type, rate_theta)

est_pars_weight3 <- est_pars_weight(res_theta_post3, res_mcmc, range)


salso_cc_5mod <- findcc_salso(pos_theta_candidate) 

res_theta_post4 <- post_sampling(data, salso_cc_5mod$cc, salso_cc_5mod$pars_unique, 
                                 salso_cc_5mod$num_pars, kappa_theta,
                                 batch_num_post, batch_size_post,
                                 pars_G0, type, rate_theta)

est_pars_weight4 <- est_pars_weight(res_theta_post4, res_mcmc, range)



res <- list(res_mcmc=res_mcmc,
            LPML=LPML,
            res_den=res_den,
            res_theta_post1=res_theta_post1,
            est_pars_weight1=est_pars_weight1,
            res_theta_post2=res_theta_post2,
            est_pars_weight2=est_pars_weight2,
            res_theta_post3=res_theta_post3,
            est_pars_weight3=est_pars_weight3,
            res_theta_post4=res_theta_post4,
            est_pars_weight4=est_pars_weight4)

filename <- paste0("job_name=", job_name,"job_num=", job_num,
                   "sample_size=", sample_size, 
                   "real_type=",stype,"res_type=",type, ".rds")
saveRDS(res, paste0(path,"/",filename))
