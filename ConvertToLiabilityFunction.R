library(MASS)
library(fMultivar)
library(calculus)
library(numDeriv)

convPearsonToTetrachoric <- function(r_01,K){
  t    <- qnorm(1-K)
  P_11 <- K*K + r_01 * K*(1-K)
  mod  <- optim(par=0,fn=function(rho){
    abs( P_11 - fMultivar::pnorm2d(-t,-t, rho = rho)[1] )
  },method = "Brent",lower=-1,upper=1)
  return(mod$par)
}

ConvertToObserved <- function(csq_l=0.2,hsq_l=0.3,K=0.1,s=0.038){
  t  <- qnorm(1-K)
  f0 <- function(x){
    pnorm2d(-t,-t, rho = csq_l + hsq_l * x)[1]* dnorm(x,mean=0.5,sd=s)
  }
  f1 <- function(x){
    x * pnorm2d(-t,-t, rho = csq_l + hsq_l * x)[1]* dnorm(x,mean=0.5,sd=s)
  }
  I0 <- integral(f=f0,bounds = list(x = c(0,1)),relTol = 1e-4,absTol=1e-16)$value
  I1 <- integral(f=f1,bounds = list(x = c(0,1)),relTol = 1e-4,absTol=1e-16)$value

  # I0 <- integrate(f=f0,lower=0,upper=1)$value
  # I1 <- integrate(f=f1,lower=0,upper=1)$value
  
  hsq_01 <- (I1-0.5*I0)/(K*(1-K)*s*s)
  csq_01 <- ((0.25+s*s)*I0-0.5*I1-s*s*K*K)/(K*(1-K)*s*s)
  return(c(csq_01=csq_01,hsq_01=hsq_01))
}

dx_  <- 0.001
xr_  <- seq(0.3,0.7,by=dx_)
pdx_ <- dnorm(xr_,mean=0.5,sd=0.038)

ConvertToObserved <- function(csq_l=0.2,hsq_l=0.3,K=0.1,s=0.038){
  t  <- qnorm(1-K)
  f0 <- function(x){
    pnorm2d(-t,-t, rho = csq_l + hsq_l * x)[1]
  }
  f1 <- function(x){
    x * pnorm2d(-t,-t, rho = csq_l + hsq_l * x)[1]
  }
  f0_x_ <- sapply(xr_,f0)
  f1_x_ <- sapply(xr_,f1)
  
  I0 <- sum(f0_x_ * pdx_) * dx_
  I1 <- sum(f1_x_ * pdx_) * dx_
  
  hsq_01 <- (I1-0.5*I0)/(K*(1-K)*s*s)
  csq_01 <- ((0.25+s*s)*I0-0.5*I1-s*s*K*K)/(K*(1-K)*s*s) 
  return(c(csq_01=csq_01,hsq_01=hsq_01))
}


ConvertToLiability <- function(csq_01=0.2,hsq_01=0.3,K=0.1,s=0.038){
  t   <- qnorm(1-K)
  Den <- K * (1 - K) * s * s
  pred_vc_01 <- function(theta_l){
    csq_l   <- theta_l[1]
    hsq_l   <- theta_l[2]
    ConvertToObserved(csq_l,hsq_l,K,s)
  }
  model <- optim(par=c(0,0),fn = function(x){
    max( abs( c(csq_01,hsq_01) - pred_vc_01(x)) )
  })
  theta_l <- model$par
  csq_l   <- theta_l[1]
  hsq_l   <- theta_l[2]
  return(c(csq_l=csq_l,hsq_l=hsq_l,Crit=model$value))
}

## numerical derivartives
## ConvFun
## Example with Smoking
K         = 0.355
hsq_01    = 0.224
se_hsq_01 = 0.0189
csq_01    = 0.0999
se_csq_01 = 0.0374
rho       = - 1 / sqrt(1 + 4 * 0.038 * 0.038)
ConvFun <- function(hsq_01,se_hsq_01,
                    csq_01,se_csq_01,
                    K,rho = - 1 / sqrt(1 + 4 * 0.038 * 0.038)){
  # cat(paste0("rho = ",rho,".\n"))
  H_ <- function(theta_01){
    ConvertToLiability(csq_01=theta_01[1],hsq_01=theta_01[2],K=K,s=0.038)
  }
  H <- function(theta_l){
    ConvertToObserved(csq_l=theta_l[1],hsq_l=theta_l[2],K=K,s=0.038)
  }
  CovMat    = matrix(c(se_csq_01,
                       rho*se_csq_01*se_hsq_01,
                       se_csq_01*se_hsq_01,
                       se_hsq_01),2,2)^2
  
  theta_01   <- c(csq_01=csq_01,hsq_01=hsq_01)
  theta_l    <- H_(theta_01)
  #jacobianH_ <- numDeriv::jacobian(func=H_,x=theta_01)
  jacobianH  <- numDeriv::jacobian(func=H,x=theta_l[1:2])
  jacobianH_ <- solve(jacobianH)
  varTheta_l <- jacobianH_%*%CovMat%*%t(jacobianH_)
  return(c(hsq_l=as.numeric(theta_l["hsq_l"]),
           se_hsq_l=sqrt(varTheta_l[2,2]),
           csq_l=as.numeric(theta_l["csq_l"]),
           se_csq_l=sqrt(varTheta_l[1,1]),
           Crit=as.numeric(theta_l["Crit"])
           ))
}

## Smoking
# test <- ConvFun(hsq_01=0.224,se_hsq_01=0.0189,csq_01=0.0999,se_csq_01=0.0374,K=0.355)
# print(test[1:4])
# 
# ## T2D
# test <- ConvFun(hsq_01=0.194,se_hsq_01=0.0401,csq_01=0.0277,se_csq_01=0.0203,K=0.049)
# print(test[1:4])
# 
# ## Asthma
# test <- ConvFun(hsq_01=0.175,se_hsq_01=0.0425,csq_01=0.0326,se_csq_01=0.0215,K=0.181)
# print(test[1:4])

# ## CAD - in EAS
# test <- ConvFun(hsq_01=-2,se_hsq_01=0.273,csq_01=1.1,se_csq_01=0.139,K=0.00825,rho=0.9)
# print(test)

			
