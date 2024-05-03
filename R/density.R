#densities of archimedean copula
#inculding AMH, CLA, FRA, GUM, JOE

generator <- function(t, theta, type="C") {
  if (type == "A") {
    log((1-theta*(1-t))/t)
  } else if (type == "C") {
    t^{-theta}-1
  } else if (type == "F") {
    -log((exp(-theta*t)-1)/(exp(-theta)-1))
  } else if (type == "G") {
    (-log(t))^theta
  } else if (type == "J") {
    -log(1-(1-t)^theta)
  } else {
    warning(paste0("implementation does not include this copula family",type))
  }
}

generator_inv <- function(t, theta, type="C") {
  if (type == "A") {
    (1-theta)/(exp(t)-theta)
  } else if (type == "C") {
    (1+t)^(-1/theta)
  } else if (type == "F") {
    -log(1-(1-exp(-theta))*exp(-t))/theta
  } else if (type == "G") {
    exp(-t^(1/theta))
  } else if (type == "J") {
    1-(1-exp(-t))^(1/theta)
  } else {
    warning(paste0("implementation does not include this copula family",type))
  }
}

generator_d1 <- function(t, theta, type="C") {
  if (type == "A") {
    (theta-1)/(t*(1-theta*(1-t)))
  } else if (type == "C") {
    -theta*t^(-theta-1)
  } else if (type == "F") {
    (theta*exp(-theta*t))/(exp(-theta*t)-1)
  } else if (type == "G") {
    -theta/t*(-log(t))^(theta-1)
  } else if (type == "J") {
    -theta*(1-t)^(theta-1)/(1-(1-t)^theta)
  } else {
    warning(paste0("implementation does not include this copula family",type))
  }
}

generator_d2 <- function(t, theta, type="C") {
  if (type == "A") {
    ((1-theta)*(1-theta+2*theta*t))/(t*(1-theta*(1-t)))^2
  } else if (type == "C") {
    theta*(theta+1)*t^(-theta-2)
  } else if (type == "F") {
    (theta^2*exp(-theta*t))/(exp(-theta*t)-1)^2
  } else if (type == "G") {
    theta/t^2*(-log(t))^(theta-1)+theta*(theta-1)/t^2*(-log(t))^(theta-2)
  } else if (type == "J") {
    (theta*(theta-1)*(1-t)^(theta-2)+theta*(1-t)^(2*theta-2))/(1-(1-t)^theta)^2
  } else {
    warning(paste0("implementation does not include this copula family",type))
  }
}

#test these functions by numeric calculations
# type<-"J"
# theta <- 2
# generator_d1(0.5, theta, type)
# delta <- 10^(-5)
# (generator(0.5+delta, theta, type)-generator(0.5, theta, type))/delta
# generator_d2(0.5, theta, type)
# (generator_d1(0.5+delta, theta, type)-generator_d1(0.5, theta, type))/delta


#Copula
copula <- function(u=c(0.5,0.5), theta, type="C") {
  generator_inv(generator(u[1], theta, type)+generator(u[2], theta, type),
                theta,type)
}

#density of each copula
dcopula <- function(u=c(0.5,0.5),theta,type="C") {
  C <- copula(u, theta, type)
  -generator_d1(u[1], theta,type)*generator_d1(u[2],theta,type)*
    generator_d2(C,theta,type)/(generator_d1(C,theta,type))^3
}


# check the density function
# library(copula)
# 
# theta = 0.5
# Cop <- amhCopula(theta, dim = 2) #amh, frank, clayton, gumbel, joe
# density_value <- dCopula(c(0.5, 0.5), copula = Cop)
# density_value
# dcopula(c(0.5,0.5),theta,"A")

