library(foreach)
library(doParallel)
rho_to_tau <- function(rho) {
  2 * asin(rho) / pi
}
# Example: for rho = 0.7
rho_to_tau(0.7)

loglik <- function(rho, u) {
  z <- qnorm(u)
  z1 <- z[,1]
  z2 <- z[,2]
  # Ensure |rho| < 1
  if (abs(rho) >= 1) return(-Inf)
  n <- length(z1)
  log_det <- -0.5 * n * log(1 - rho^2)
  quad_form <- sum( (rho^2*z1^2 + rho^2*z2^2 - 2 * rho * z1 * z2) ) / (2 * (1 - rho^2))
  # Constant term: -n*log(2π) (can be ignored for MH since it's constant)
  return(log_det - quad_form)
}

loglik2 <- function(rho, u) {
  z <- qnorm(u)
  z1 <- z[1]
  z2 <- z[2]
  # Ensure |rho| < 1
  
  log_det <- -0.5 * log(1 - rho^2)
  quad_form <-  (rho^2*z1^2 + rho^2*z2^2 - 2 * rho * z1 * z2) / (2 * (1 - rho^2))
  # Constant term: -n*log(2π) (can be ignored for MH since it's constant)
  return(log_det - quad_form)
}

############################################################
# 3. METROPOLIS-HASTINGS SAMPLER FOR ρ ~ U(-1,1)
############################################################

mh_sampler <- function(u, n_iter = 10000, init = 0, proposal_sd = 0.05) {
  rho_samples <- numeric(n_iter)
  rho_current <- init
  loglik_current <- loglik(rho_current, u)
  
  for (i in 1:n_iter) {
    # Propose new ρ (random walk)
    rho_proposed <- runif(1,min=-1,max=1)
    loglik_proposed <- loglik(rho_proposed, u)
    # With uniform prior, posterior ∝ likelihood. Compute log acceptance ratio.
    log_ratio <- loglik_proposed - loglik_current
    
    if (log(runif(1)) < log_ratio) {
      # Accept the proposal
      rho_current <- rho_proposed
      loglik_current <- loglik_proposed
    }
    rho_samples[i] <- rho_current
  }
  return(rho_samples)
}

# Run the sampler
set.seed(456)
U <- readRDS("../../Data/Gaussian_sim_data/single.rds")
batch_num <- 170
batch_size <- 50
burn_in_batch <- 100
thin <- 5
n_iter <- batch_num*batch_size 

samples <- mh_sampler(U, n_iter = n_iter, init = 0)

# Remove burn-in 
burn_in <-burn_in_batch*batch_size
range <- seq(burn_in, n_iter, by=thin)
posterior_samples <- samples[range]

round(c(mean(posterior_samples), 
  quantile(posterior_samples, 0.025),
  quantile(posterior_samples, 0.975)),3)

##LPML
CPO <- lapply(1:nrow(U), function(k) {
  mean(1/exp(loglik2(posterior_samples, U[k,])))
})
sum(-log(do.call(rbind,CPO)))

#estimated tau with 95% CI
tau <-rho_to_tau(posterior_samples) 
round(c(mean(tau), 
  quantile(tau, 0.025),
  quantile(tau, 0.975)),3)

set.seed(456)
U <-readRDS("../../Data/Gaussian_sim_data/mix.rds")
batch_num <- 170
batch_size <- 50
burn_in_batch <- 100
thin <- 5
n_iter <- batch_num*batch_size 

samples <- mh_sampler(U, n_iter = n_iter, init = 0)

# Remove burn-in 
burn_in <-burn_in_batch*batch_size
range <- seq(burn_in, n_iter, by=thin)
posterior_samples <- samples[range]

c(mean(posterior_samples), 
  quantile(posterior_samples, 0.025),
  quantile(posterior_samples, 0.975))

##LPML
CPO <- lapply(1:nrow(U), function(k) {
  mean(1/exp(loglik2(posterior_samples, U[k,])))
})
#estimated tau with 95% CI
sum(-log(do.call(rbind,CPO)))
tau <-rho_to_tau(posterior_samples) 
round(c(mean(tau), 
  quantile(tau, 0.025),
  quantile(tau, 0.975)),3)
