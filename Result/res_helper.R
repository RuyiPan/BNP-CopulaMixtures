#res_helper

#Get kendall's tau


kendalls_tau_joe <- function(theta, n_terms = 10000) {
  sum_term <- 0
  for (k in 1:n_terms) {
    sum_term <- sum_term + 1 / (k * (theta * k + 2) * (theta * (k - 1) + 2))
  }
  tau <- 1 - 4 * sum_term
  return(tau)
}
type_num <- list("A"=9, "F"=5, "C"=3, "J"=7, "G"=1)

getTau <- function(theta, type) {
  type_num <- list("A"=9, "F"=5, "C"=3, "J"=7, "G"=1)
  if (type == "A") {
    1-2*{theta+(1-theta)^2*log(1-theta)}/(3*theta^2)
  } else if (type == "F") {
    do.call(rbind, lapply(theta, function(the){
      fc <- frankCopula(param = the)
      tau_val <- tau(fc)
      tau_val
    }))
                           
    
    # Kendall's tau:
    
    # do.call(rbind,lapply(theta, function(the){1-4/the*
    #     (1-integrate(function(t) t / (exp(t) - 1), lower = 0, upper = the)$value/the)}))
  } else if (type == "C"){
    theta/(theta+2)
  } else if (type == "G") {
    1- 1/theta
  } else {
    kendalls_tau_joe(theta)
  }
}


# Install if necessary
# install.packages("copula")



