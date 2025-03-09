library(MASS)
set.seed(101)  
n1  <- 500     # sample size for Scenario 1
rho <- 0.6     # correlation in the bivariate normal

Sigma1 <- matrix(c(1, rho,
                   rho, 1),
                 nrow = 2, ncol = 2)

Z1 <- mvrnorm(n = n1, mu = c(0, 0), Sigma = Sigma1)

# Convert to [0,1] scale using the standard normal CDF (pnorm)
U1 <- pnorm(Z1)
#or
nC1 <- normalCopula(0.7, dim = 2, dispstr = "ex")
U1 <- rCopula(n1, nC1)
plot(U1[,1], U1[,2])
saveRDS(U1,file="../../Data/Gaussian_sim_data/single.rds")

set.seed(202)  # different seed

# Parameters
n2   <- 500
rho1 <- 0.7    # correlation in component 1
rho2 <- -0.7   # correlation in component 2
pi1  <- 0.5    # mixture weight for component 1
pi2  <- 1 - pi1

# Cov matrices for each component
Sigma2A <- matrix(c(1,  rho1,
                    rho1, 1),
                  nrow = 2, ncol = 2)
Sigma2B <- matrix(c(1,  rho2,
                    rho2, 1),
                  nrow = 2, ncol = 2)

# Allocate space
Z2 <- matrix(NA, nrow = n2, ncol = 2)

# Generate component labels
# labels = 1 => component 1, labels = 2 => component 2
labels <- sample(1:2, size = n2, replace = TRUE, prob = c(pi1, pi2))

# Draw each data point from its assigned bivariate normal
for (i in seq_len(n2)) {
  if (labels[i] == 1) {
    Z2[i, ] <- mvrnorm(n = 1, mu = c(0,0), Sigma = Sigma2A)
  } else {
    Z2[i, ] <- mvrnorm(n = 1, mu = c(0,0), Sigma = Sigma2B)
  }
}

# Convert to copula scale
U2 <- pnorm(Z2)

#or 
labels <- sample(1:2, size = n2, replace = TRUE, prob = c(pi1, pi2))
nC1 <- normalCopula(0.7, dim = 2, dispstr = "ex")
nC2 <- normalCopula(-0.7, dim = 2, dispstr = "ex")

for (i in seq_len(n2)) {
  if (labels[i] == 1) {
    U2[i, ] <- rCopula(1, nC1)
  } else {
    U2[i, ] <- rCopula(1, nC2)
  }
}
# Check a few rows
head(U2)

plot(U2[,1], U2[,2])

saveRDS(U2,file="../../Data/Gaussian_sim_data/mix.rds")


library(copula)


set.seed(123)  # for reproducibility
samples <- rCopula(500, nC1)
plot(samples[,1], samples[,2])


samples1 <- rCopula(250, nC1)
samples2 <- rCopula(250, nC2)
sample_mix <- rbind(samples1,samples2)
plot(sample_mix[,1], sample_mix[,2])
u=c(0.1,0.9) 
rho=-0.7
dCopula(u, nC2, log = T)
loglik(rho, u)
