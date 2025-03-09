#function for the mcmc sampling procedure
rG0 <- function(N, pars, type="C") {
  if (type == "A") {
    runif(N, min=-1, max=1)
  } else if (type == "C") {
    rgamma(N, shape=pars[1], rate=pars[2])-1
  } else if (type == "F") {
    rnorm(N,mean=pars[1],sd=pars[2])
  } else if (type == "G") {
    rgamma(N, shape=pars[1], rate=pars[2])+1
  } else if (type == "J") {
    rgamma(N, shape=pars[1], rate=pars[2])+0.3
  } else {
    warning(paste0("implementation does not include this copula family",type))
  }
}


dG0 <- function(theta, pars, type="C") {
  if (type == "A") {
    dunif(theta, min=-1, max=1)
  } else if (type == "C") {
    dgamma(theta+1, shape=pars[1], rate=pars[2])
  } else if (type == "F") {
    dnorm(theta,mean=pars[1],sd=pars[2])
  } else if (type == "G") {
    dgamma(theta-1, shape=pars[1], rate=pars[2])
  } else if (type == "J") {
    dgamma(theta-0.3, shape=pars[1], rate=pars[2])
  } else {
    warning(paste0("implementation does not include this copula family",type))
  }
}

#test rGcenter
# rGcenter(c(0,0.1),"F")
log_dG0 <- function(theta, pars, type="C") {
  if (type == "A") {
    dunif(theta, min=-1, max=1, log = T)
  } else if (type == "C") {
    dgamma(theta+1, shape=pars[1], rate=pars[2], log = T)
  } else if (type == "F") {
    dnorm(theta,mean=pars[1],sd=pars[2], log = T)
  } else if (type == "G") {
    dgamma(theta-1, shape=pars[1], rate=pars[2], log = T)
  } else if (type == "J") {
    dgamma(theta-0.3, shape=pars[1], rate=pars[2], log = T)
  } else {
    warning(paste0("implementation does not include this copula family",type))
  }
}


#helper function for clustering configuration
f_cc <- function(data, theta, pars_G0, type) {
  # print(sapply(1:nrow(data), function(i) dcopula(data[i,], theta, type)))
  dG0(theta, pars_G0, type)*
    prod(sapply(1:nrow(data), function(x) dcopula_new(data[x,], theta, type)))
  
}

# f_cc(matrix(runif(2*10), 10, 2), 3, c(1,5), "C")
# f_cc(matrix(runif(2*10), 10, 2), 3.1, c(1,5), "C")
log_f_cc <- function(data, theta, pars_G0, type) {
  # print(sapply(1:nrow(data), function(i) dcopula(data[i,], theta, type)))
  log_dG0(theta, pars_G0, type)+
    sum(sapply(1:nrow(data), function(x) log(dcopula_new(data[x,], theta, type))))
  
}

#propose new theta
rtheta_cc <- function(theta, kappa_theta, type) {
  if (type == "A") {
    runif(1, min=max(-1, theta-kappa_theta), max=min(1, theta+kappa_theta))
  } else if (type == "C") {
    runif(1, min=max(-1, theta-kappa_theta), max=theta+kappa_theta)
  } else if (type == "F") {
    runif(1, min=theta-kappa_theta, max=theta+kappa_theta)
  } else if (type == "G") {
    runif(1, min=max(1, theta-kappa_theta), max=theta+kappa_theta)
  } else if (type == "J") {
    runif(1, min=max(0.3, theta-kappa_theta), max=theta+kappa_theta)
  } else {
    warning(paste0("implementation does not include this copula family",type))
  }
}

# rtheta_cc(1, 0.1, "C")
# rtheta_cc(1, 0.1, "F")

#hepler function for posterior, a, and b
f_a <- function(a, b, cluster_size, pars_a0, pars_b0) {
  m <- length(cluster_size)
  # print(b+(1:(m-1))*a)
  # print(gamma(cluster_size-a)/gamma(1-a))
  prod(b+(1:(m-1))*a)*prod(gamma(cluster_size-a)/gamma(1-a))*
    dgamma(b+a, pars_b0[1], pars_b0[2])*
    dbeta(a, pars_a0[1], pars_a0[2])
}

log_f_a <- function(a, b, cluster_size, pars_a0, pars_b0) {
  m <- length(cluster_size)
  # print(b+(1:(m-1))*a)
  # print(gamma(cluster_size-a)/gamma(1-a))
  # print(sum(log(b+(1:(m-1))*a)))
  # print(sum(lgamma(cluster_size-a)-lgamma(1-a)))
  # print(dgamma(b+a, pars_b0[1], pars_b0[2], log=T))
  # print(dbeta(a, pars_a0[1], pars_a0[2],log=T))
  sum(log(b+(1:(m-1))*a))+sum(lgamma(cluster_size-a)-lgamma(1-a))+
    dgamma(b+a, pars_b0[1], pars_b0[2], log=T)+
    dbeta(a, pars_a0[1], pars_a0[2],log=T)
}

f_b <- function(a, b, n, m, pars_b0) {
  gamma(b+1)/gamma(b+n)*prod(b+(1:(m-1))*a)*
    dgamma(b+a, pars_b0[1], pars_b0[2])
}

log_f_b  <- function(a, b, n, m, pars_b0) {
  lgamma(b+1)-lgamma((b+n))+sum(log(b+(1:(m-1))*a))+
    dgamma(b+a, pars_b0[1], pars_b0[2],log=T)
}
# f_a(0.5, 1, c(10,50,20,30,40), c(1,5),c(2,3))
# f_b(0.5, 1, 100, 5, c(2,3))

