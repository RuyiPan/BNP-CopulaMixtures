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
  # print(generator(u[1], theta, type))
  # print(generator(u[2], theta, type))
  generator_inv(generator(u[1], theta, type)+generator(u[2], theta, type),
                theta,type)
}


copula_new <- function(u=c(0.5,0.5), theta, type="C") {
  # print(generator(u[1], theta, type))
  # print(generator(u[2], theta, type))
  if(type=="J") {
    1-((1-u[1])^theta+(1-u[2])^theta-(1-u[1])^theta*(1-u[2])^theta)^{1/theta}
  }else if (type=="A") {
    (1-theta)/((1-theta*(1-u[1]))*(1-theta*(1-u[2]))/u[1]/u[2]-theta)
  } else if (type=="C") {
    (u[1]^(-theta)+u[2]^(-theta)-1)^{-1/theta}
  }
}
#density of each copula
dcopula_old <- function(u=c(0.5,0.5),theta,type="C") {
  C <- copula(u, theta, type)
  # print(C)
 
  dens <- -generator_d1(u[1], theta,type)*generator_d1(u[2],theta,type)*
      generator_d2(C,theta,type)/(generator_d1(C,theta,type))^3
  # if (is.na(dens) | dens==Inf) {
  #   if (type == "A") {
  #     Cop <- amhCopula(theta, dim = 2)
  #   } else if (type == "C") {
  #     Cop <- claytonCopula(theta, dim = 2)
  #   } else if (type == "F") {
  #     Cop <- frankCopula(theta, dim = 2)
  #   } else if (type == "G") {
  #     Cop <- gumbelCopula(theta, dim = 2)
  #   } else if (type == "J") {
  #     Cop <- joeCopula(theta, dim = 2)
  #   } else {
  #     warning(paste0("implementation does not include this copula family",type))
  #   }
  #   dens <- dCopula(u, copula = Cop)
  # }
  dens
  
}
dcopula_new <- function(u=c(0.5,0.5),alpha,type="C") {
  u1 <- u[1]
  u2 <- u[2]
  # print(C)
  if (u[1]==0 | u[2] == 0 | 1-u[1]==0 | 1-u[2]==0){
    return(0)
  }
  
  if (type == "A") {
    # (1-theta+2*theta*C)*C*(1-theta*(1-C))/
    #  (u[1]*(1-theta*(1-u[1])))/(u[2]*(1-theta*(1-u[2])))
   dense<- -((1 + alpha^2 * (-1 + u1) * (-1 + u2) + alpha * 
         (-2 + u1 + u2 + u1 * u2))/(-1 + alpha * (-1 + u1) * (-1 + u2))^3)
   if (dense==Inf) {
     dim <- length(u)
     Cop <- amhCopula(alpha, dim =dim)
     dense <- dCopula(u, copula = Cop)
   }
  } else if (type == "C") {
    
   
    
    if (length(u) > 2) {
        if (alpha <= 0) {
          dense <- 0
        } else {
           dense <- dClayton(u, alpha)
           if (is.na(dense)) {
             dense <- 0
           } else if (dense==Inf) {
             dim <- length(u)
             Cop <- claytonCopula(alpha, dim =dim)
             dense <- dCopula(u, copula = Cop)
           }
        }
       
    } else {
      if ((-1 + u1^(-alpha) + u2^(-alpha)) < 0 & !is.integer(alpha)) {
        dense <- 0 
      } else {
        dense<- (1 + alpha) * u1^(-1 -alpha) * u2^(-1 - alpha) * (-1 + u1^(-alpha) + u2^(-alpha))^(-2 -  1/alpha)
      }
    }
    
  } else if (type == "F") {
    dense<- dFrank(u, alpha)
    if (is.na(dense)) {
      dense <- 0
    } else if (dense==Inf) {
      dim <- length(u)
      Cop <- frankCopula(alpha, dim =dim)
      dense <- dCopula(u, copula = Cop)
    }
    
  } else if (type == "G") {
    dense<- dGumbel(u, alpha)
    if (is.na(dense)) {
      dense <- 0
    } else if (dense==Inf) {
      dim <- length(u)
      Cop <- gumbelCopula(alpha, dim =dim)
      dense <- dCopula(u, copula = Cop)
    }
  } else if (type == "J") {
    # (1-u[1])^(theta-1)*(1-u[2])^(theta-1)*
    #   ((theta-1)*((1-u[1])^theta+(1-u[2])^theta-(1-u[1])^theta*(1-u[2])^theta)^(1/theta-2)+
    #      ((1-u[1])^theta+(1-u[2])^theta-(1-u[1])^theta*(1-u[2])^theta)^(1/theta-1))
    
    dense<-(1-u1)^(alpha-1)*(1-u2)^(alpha-1)*
      ((alpha-1)*((1-u1)^alpha+(1-u2)^alpha-(1-u1)^alpha*(1-u2)^alpha)^(1/alpha-2)+
         ((1-u1)^alpha+(1-u2)^alpha-(1-u1)^alpha*(1-u2)^alpha)^(1/alpha-1))
    if (dense < 0) {
      dense <- 0
    }
  } else {
    warning(paste0("implementation does not include this copula family",type))
  }
  
  # if (is.na(dense)){
  #   if (type == "A") {
  #     Cop <- amhCopula(alpha, dim = 2)
  #   } else if (type == "C") {
  #     Cop <- claytonCopula(alpha, dim = 2)
  #   } else if (type == "F") {
  #     Cop <- frankCopula(alpha, dim = 2)
  #   } else if (type == "G") {
  #     Cop <- gumbelCopula(alpha, dim = 2)
  #   } else if (type == "J") {
  #     Cop <- joeCopula(alpha, dim = 2)
  #   } else {
  #     warning(paste0("implementation does not include this copula family",type))
  #   }
  #   
  #   dense <- dCopula(u, copula = Cop)
  #   print("repalce NA")
  #   print(dense)
  # }
  if (is.na(dense)) {
    dense <- 0
  }
  dense
}

# #dcoup
# type="C"
# u=c(0.5,0.5)
# alpha=-0.3
# dcopula_new2(u,theta,type)
# dcopula(u,theta,type)
dcopula <- function(u=c(0.5,0.5),theta,type="C") {
  if (type == "A") {
    Cop <- amhCopula(theta, dim = 2)
  } else if (type == "C") {
    Cop <- claytonCopula(theta, dim = 2)
  } else if (type == "F") {
    Cop <- frankCopula(theta, dim = 2)
  } else if (type == "G") {
    Cop <- gumbelCopula(theta, dim = 2)
  } else if (type == "J") {
    Cop <- joeCopula(theta, dim = 2)
  } else {
    warning(paste0("implementation does not include this copula family",type))
  }
  dCopula(u, copula = Cop)
}


# check the density function
# library(copula)
# 
# theta = -0.4249283
# Cop <- claytonCopula(theta, dim = 2) #amh, frank, clayton, gumbel, joe
# density_value <- dCopula(c(0.0455565, 0.1585318), copula = Cop)
# density_value
# dcopula(c(0.5,0.5),theta,"A")
# 0.0455565 0.1585318
# theta = -0.4249283
# 0.9558152 0.9607128
# 12.60753
# dcopula_old(c(0.9558152, 0.9607128), 12.60753, "J")
# copula(c(0.9558152, 0.9607128), 12.60753, "J")
# # generator_d2(1,12.60753, "J")
# # generator_d1(1,12.60753, "J")
# generator(0.9558152, 12.60753, "J")
# 
# Cop <- joeCopula(12.60753, dim = 2) #amh, frank, clayton, gumbel, joe
# density_value <- dCopula(c(0.9558152, 0.9607128), copula = Cop)
# density_value
dClayton <- function(u, alpha) {
  dim <- length(u)
  
    prod(1+alpha*c(1:(dim-1)))*prod(u^{-1-alpha})*(-(dim-1)+sum(u^{-alpha}))^(- dim -  1/alpha)
}                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                             

# dim<- 8
# theta <- 22.33639
# Cop <- claytonCopula(theta, dim =dim)
# u <- rep(0.1,dim)
# dCopula(u, copula = Cop)
# dClayton(u, theta )


# Define the general function for dimension n
dGumbel<- function(u, alpha) {

  n <- length(u)
  # Calculate the logarithmic terms
  log_u <- -log(u)
  
  if (n == 1) {
    # Specific case for n = 1
    result <- -(((log_u[1]^alpha)^(1/alpha)) / (exp((log_u[1]^alpha)^(1/alpha)) * u[1] * log(u[1])))
  } else if (n == 2) {
    # Specific case for n = 2
    term1 <- log_u[1]^(-1 + alpha)
    term2 <- -1 + alpha + (log_u[1]^alpha + log_u[2]^alpha)^(1/alpha)
    term3 <- (log_u[1]^alpha + log_u[2]^alpha)^(-2 + 1/alpha)
    term4 <- log_u[2]^(-1 + alpha)
    denominator <- exp((log_u[1]^alpha + log_u[2]^alpha)^(1/alpha)) * u[1] * u[2]
    result <- (term1 * term2 * term3 * term4) / denominator
  } else if (n == 3) {
    # Specific case for n = 3
    term1 <- log_u[1]^(-1 + alpha)
    term2 <- log_u[2]^(-1 + alpha)
    term3 <- 1 + 2 * alpha^2 + 3 * alpha * (-1 + (log_u[1]^alpha + log_u[2]^alpha + log_u[3]^alpha)^(1/alpha)) -
      3 * (log_u[1]^alpha + log_u[2]^alpha + log_u[3]^alpha)^(1/alpha) +
      (log_u[1]^alpha + log_u[2]^alpha + log_u[3]^alpha)^(2/alpha)
    term4 <- (log_u[1]^alpha + log_u[2]^alpha + log_u[3]^alpha)^(-3 + 1/alpha)
    term5 <- log_u[3]^(-1 + alpha)
    denominator <- exp((log_u[1]^alpha + log_u[2]^alpha + log_u[3]^alpha)^(1/alpha)) * u[1] * u[2] * u[3]
    result <- (term1 * term2 * term3 * term4 * term5) / denominator
  } else if (n == 4) {
    # Specific case for n = 4
    term1 <- log_u[1]^(-1 + alpha)
    term2 <- log_u[2]^(-1 + alpha)
    term3 <- log_u[3]^(-1 + alpha)
    term4 <- -1 + 6 * alpha^3 + 11 * alpha^2 * (-1 + (log_u[1]^alpha + log_u[2]^alpha + log_u[3]^alpha + log_u[4]^alpha)^(1/alpha)) +
      6 * alpha * (1 - 3 * (log_u[1]^alpha + log_u[2]^alpha + log_u[3]^alpha + log_u[4]^alpha)^(1/alpha) +
                     (log_u[1]^alpha + log_u[2]^alpha + log_u[3]^alpha + log_u[4]^alpha)^(2/alpha)) +
      7 * (log_u[1]^alpha + log_u[2]^alpha + log_u[3]^alpha + log_u[4]^alpha)^(1/alpha) -
      6 * (log_u[1]^alpha + log_u[2]^alpha + log_u[3]^alpha + log_u[4]^alpha)^(2/alpha) +
      (log_u[1]^alpha + log_u[2]^alpha + log_u[3]^alpha + log_u[4]^alpha)^(3/alpha)
    term5 <- (log_u[1]^alpha + log_u[2]^alpha + log_u[3]^alpha + log_u[4]^alpha)^(-4 + 1/alpha)
    term6 <- log_u[4]^(-1 + alpha)
    denominator <- exp((log_u[1]^alpha + log_u[2]^alpha + log_u[3]^alpha + log_u[4]^alpha)^(1/alpha)) * u[1] * u[2] * u[3] * u[4]
    result <- (term1 * term2 * term3 * term4 * term5 * term6) / denominator
  } else {
    # General case for n > 4 (using a similar pattern)
    message("Haven't implemented, you can refer the HAC package")
  }
  
  return(result)
}



# Define the function
dFrank <- function(u, alpha) {
  n <- length(u)
  E <- exp(alpha)
  
  if (n == 2) {
    u1 <- u[1]
    u2 <- u[2]
    num <- alpha * exp(alpha * (1 + u1 + u2)) * (-1 + E)
    den <- (E - exp(alpha + alpha * u1) + exp(alpha * (u1 + u2)) - exp(alpha + alpha * u2))^2
    return(num / den)
  } else if (n == 3) {
    if (alpha <= 0) {
      return (0)
    } else {
      u1 <- u[1]
      u2 <- u[2]
      u3 <- u[3]
      num <- alpha^2 * exp(alpha * (-2 + u1 + u2 + u3)) * (-1 + E)^2 *
        (-1 + exp(alpha * u1) + exp(alpha * u2) - exp(alpha * (u1 + u2)) + exp(alpha * u3) -
           exp(alpha * (u1 + u3)) - exp(alpha * (u2 + u3)) + exp(alpha * (-2 + u1 + u2 + u3)) -
           2 * exp(alpha * (-1 + u1 + u2 + u3)) + 2 * exp(alpha * (u1 + u2 + u3)))
      den <- (1 - exp(alpha * u1) - exp(alpha * u2) + exp(alpha * (u1 + u2)) - exp(alpha * u3) +
                exp(alpha * (u1 + u3)) + exp(alpha * (u2 + u3)) + exp(alpha * (-2 + u1 + u2 + u3)) -
                2 * exp(alpha * (-1 + u1 + u2 + u3)))^3
      return(num / den)
    }
    
  } else if (n == 4) {
    u1 <- u[1]
    u2 <- u[2]
    u3 <- u[3]
    u4 <- u[4]
    num <- alpha^3 * exp(alpha * (1 + u1 + u2 + u3 + u4)) * (-1 + exp(-alpha))^4 * (-1 + E) *
      (1 - 2 * E + exp(2 * alpha) + (exp(2 * alpha * (-3 + u1 + u2 + u3 + u4)) * (-1 + E)^8) /
         ((-1 + exp(alpha * u1))^2 * (-1 + exp(alpha * u2))^2 * (-1 + exp(alpha * u3))^2 * (-1 + exp(alpha * u4))) -
         (4 * exp(alpha * (-3 + u1 + u2 + u3 + u4)) * (-1 + E)^4) /
         ((-1 + exp(alpha * u1)) * (-1 + exp(alpha * u2)) * (-1 + exp(alpha * u3)) * (-1 + exp(alpha * u4))) +
         (4 * exp(alpha * (-2 + u1 + u2 + u3 + u4)) * (-1 + E)^4) /
         ((-1 + exp(alpha * u1)) * (-1 + exp(alpha * u2)) * (-1 + exp(alpha * u3)) * (-1 + exp(alpha * u4))))
    den <- ((-1 + exp(alpha * u1))^2 * (-1 + exp(alpha * u2))^2 * (-1 + exp(alpha * u3))^2 * (-1 + exp(alpha * u4))^2 *
              (1 - E + (exp(alpha * (-3 + u1 + u2 + u3 + u4)) * (-1 + E)^4) /
                 ((-1 + exp(alpha * u1)) * (-1 + exp(alpha * u2)) * (-1 + exp(alpha * u3)) * (-1 + exp(alpha * u4))))^4)
    return(num / den)
  } else if (n == 5) {
    u1 <- u[1]
    u2 <- u[2]
    u3 <- u[3]
    u4 <- u[4]
    u5 <- u[5]
    num <- alpha^4 * exp(alpha * (1 + u1 + u2 + u3 + u4 + u5)) * (-1 + exp(-alpha))^5 * (-1 + E) *
      (-1 + 3 * E - 3 * exp(2 * alpha) + exp(3 * alpha) +
         (11 * exp(3 * alpha) * (-1 + exp(-alpha))^10) /
         ((-1 + exp(-(alpha * u1)))^2 * (-1 + exp(-(alpha * u2)))^2 * (-1 + exp(-(alpha * u3)))^2 * (-1 + exp(-(alpha * u4)))^2 * (-1 + exp(-(alpha * u5)))^2) +
         (exp(3 * alpha * (-4 + u1 + u2 + u3 + u4 + u5)) * (-1 + E)^15) /
         ((-1 + exp(alpha * u1))^3 * (-1 + exp(alpha * u2))^3 * (-1 + exp(alpha * u3))^3 * (-1 + exp(alpha * u4))^3 * (-1 + exp(alpha * u5))^3) -
         (11 * exp(2 * alpha * (-4 + u1 + u2 + u3 + u4 + u5)) * (-1 + E)^10) /
         ((-1 + exp(alpha * u1))^2 * (-1 + exp(alpha * u2))^2 * (-1 + exp(alpha * u3))^2 * (-1 + exp(alpha * u4))^2 * (-1 + exp(alpha * u5))^2) +
         (11 * exp(alpha * (-4 + u1 + u2 + u3 + u4 + u5)) * (-1 + E)^5) /
         ((-1 + exp(alpha * u1)) * (-1 + exp(alpha * u2)) * (-1 + exp(alpha * u3)) * (-1 + exp(alpha * u4)) * (-1 + exp(alpha * u5))) -
         (22 * exp(alpha * (-3 + u1 + u2 + u3 + u4 + u5)) * (-1 + E)^5) /
         ((-1 + exp(alpha * u1)) * (-1 + exp(alpha * u2)) * (-1 + exp(alpha * u3)) * (-1 + exp(alpha * u4)) * (-1 + exp(alpha * u5))) +
         (11 * exp(alpha * (-2 + u1 + u2 + u3 + u4 + u5)) * (-1 + E)^5) /
         ((-1 + exp(alpha * u1)) * (-1 + exp(alpha * u2)) * (-1 + exp(alpha * u3)) * (-1 + exp(alpha * u4)) * (-1 + exp(alpha * u5))))
    den <- ((-1 + exp(alpha * u1))^2 * (-1 + exp(alpha * u2))^2 * (-1 + exp(alpha * u3))^2 * (-1 + exp(alpha * u4))^2 * (-1 + exp(alpha * u5))^2 *
              (1 - E + (exp(alpha * (-4 + u1 + u2 + u3 + u4 + u5)) * (-1 + E)^5) /
                 ((-1 + exp(alpha * u1)) * (-1 + exp(alpha * u2)) * (-1 + exp(alpha * u3)) * (-1 + exp(alpha * u4)) * (-1 + exp(alpha * u5))))^5)
    return(num / den)
  } else {
    stop("Invalid value for n. Only n = 1, 2, 3, 4 are supported.")
  }
}


