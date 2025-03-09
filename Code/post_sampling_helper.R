#post_sampling_helper
library(salso)
Mode <- function(x) {
  ux <- unique(x)
  ux[which.max(tabulate(match(x, ux)))]
}

compare_all_pairs <- function(vector) {
  n <- length(vector)
  results <- matrix(FALSE, nrow = n, ncol = n)  # Initialize a matrix to store results
  
  for (i in 1:(n - 1)) {
    for (j in (i + 1):n) {
      results[i, j] <- vector[i] == vector[j]
      results[j, i] <- results[i, j]  # Since comparison is symmetric
    }
  }
  
  # Optional: Fill diagonal with TRUE since any element is equal to itself
  diag(results) <- TRUE
  
  return(results)
}


averge_weight_matrix <- function(pos_theta) {
  result_matrix <- list()
  avg_matrix <- matrix(0, nrow=length(pos_theta[1,]), ncol=length(pos_theta[1,]))
  for (i in c(1:nrow(pos_theta))) {
    result_matrix[[i]] <- compare_all_pairs(pos_theta[i,])
    avg_matrix <- avg_matrix + result_matrix[[i]]
  }
  
  avg_weight_matrix <- avg_matrix/nrow(pos_theta)
  
  avg_weight_matrix 
}
#find the index for a good cluster
find_cluster <- function(pos_theta, avg_weight_matrix) {
  
  result_matrix <- list()
  for (i in c(1:nrow(pos_theta))) {
    result_matrix[[i]] <- compare_all_pairs(pos_theta[i,])
  }
  distance <- vector()
  for (i in c(1:nrow(pos_theta))) {
    distance[i] <- sum((result_matrix[[i]] - avg_weight_matrix )^2)
  }
  
  which(distance==min(distance))
}


findcc_old <- function(pos_theta) {
  avg_weight <- averge_weight_matrix(pos_theta)
  index_all <- find_cluster(pos_theta, avg_weight)
  index_small_all <- index_all[which(length(unique(pos_theta[index_all,]))==min(length(unique(pos_theta[index_all,]))))]
  index_final <- index_small_all[length(index_small_all)] #minimal
  cc <- list() #store the pars and corresponding indices
  pars_unique <- sort(unique(pos_theta[index_final,]))
  num_pars <- length(pars_unique)
  for (i in 1:num_pars) {
    cc[[i]]=which(pos_theta[index_final,]==pars_unique[i])
  }
  return(list(cc=cc, pars_unique=pars_unique, num_pars=num_pars))
}

findcc_salso <- function(pos_theta, method="binder") {
  result <- salso(pos_theta, loss = method)
  result_counts <- table(result)
  
  #Sort the frequency table in ascending order (lowest count first)
  sorted_counts <- sort(result_counts, decreasing = T)
  unique_values <- as.numeric(names(sorted_counts))
  cc_num <- length(unique(result))
  cc <- list()
  pars_unique <- c()
  for (i in 1:cc_num) {
    cc[[i]] <- which(result==unique_values[i])
    pars_unique[i] <- mean(pos_theta[,which(result==unique_values[i])])
  }
  num_pars <- length(pars_unique)
  return(list(cc=cc, pars_unique=pars_unique, num_pars=num_pars))
}



post_sampling <- function(data, cc, pars_start, num_pars, kappa_theta,
                          batch_num_post, batch_size_post, batch_burn_in,
                          pars_G0, type, rate_theta) {
  all_pos_theta <- pars_start #record all the sample
  kappa_theta_pos <- rep(kappa_theta, num_pars)
  kappa_theta_pos_all <- NULL
  current_theta <-  pars_start
  for (k in 1:batch_num_post) {
    print(k)
    count2_theta <- rep(0, num_pars)
    for (i in 1:batch_size_post) {
      old_theta <- current_theta
      
      #number of clusters
      # count1_theta <- count1_theta + m # increase number of proposal
      for (j in 1:num_pars) {
        theta_old_j <- current_theta[j]
        cc_index <- cc[[j]]
        data_cc <- data[cc_index, , drop=FALSE]
        theta_prop_j <- rtheta_cc(theta_old_j, kappa_theta_pos[j], type)
        
        acc_prob_theta <- min(1, exp(log_f_cc(data_cc, theta_prop_j, pars_G0, type)-
                                       log_f_cc(data_cc, theta_old_j, pars_G0, type)))
        if (runif(1) <= acc_prob_theta ) {
          current_theta[j] <- theta_prop_j
          count2_theta[j] <- count2_theta[j] + 1 # increase number of accepted proposal
        }
        
      }
      all_pos_theta<- rbind(all_pos_theta, current_theta)
    }
    #acceptance rate for theta
    acc_theta_pos <- count2_theta/batch_size_post
    index_over <- which(acc_theta_pos>0.5)
    index_lower <- which(acc_theta_pos<0.4)
    kappa_theta_pos[index_over] <- kappa_theta_pos[index_over]*rate_theta^(sqrt(k))
    kappa_theta_pos[index_lower] <- kappa_theta_pos[index_lower]*rate_theta^(-sqrt(k))
    kappa_theta_pos_all <- rbind(kappa_theta_pos_all, kappa_theta_pos)
  }
  print("Approach 2: theta_post")
  range <- (batch_burn_in*batch_size_post):(batch_num_post*batch_size_post)
  res_theta_post <- list(all_pos_theta=all_pos_theta[range,,drop=F], cc=cc, kappa_theta_pos_all=kappa_theta_pos_all)
  res_theta_post
}


est_pars_weight <- function(res_theta_post, res_mcmc, range) {
       est_pars <- apply(res_theta_post$all_pos_theta,2,
                                       function(x) {
                                         c(mean(x), quantile(x, c(0.025, 0.975)))})
                     
       as <- res_mcmc$all_a[range]
                     
       bs <- res_mcmc$all_b[range]
       n <- nrow(res_mcmc$data)
       num_pars <- ncol(res_theta_post$all_pos_theta)
       weight <- vector()
       for (j in 1:num_pars) {
         nj <- length(res_theta_post$cc[[j]])
         weight[j] <- mean((nj-as)/(bs+n))
       }
       order_indices <- order(weight,decreasing = TRUE)
       
                     # Sort the matrix by the specified row
        rbind(est_pars, weight)[, order_indices]
}