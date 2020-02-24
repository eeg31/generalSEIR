#' get equilibrium proportions
#' @include getBeta.R
#' @include getProps.R
#' @export getEquil


getEquil <- function (mod=NULL, R0=NULL, beta=NULL,
                      omega=0, gamma=0, rho=0, sigma=0, epsilon=0,
                      submodel=NULL,
                      N=1, m=NULL, b){
  if(is.null(R0) & is.null(beta)) {
    print('either R0 or beta must be nonzero')
    return()
  }

  if(is.null(mod) & is.null(submodel)) {
    print('must specify either a model ID or input a submodel object')
    return()
  }

  if(is.null(m)) {
    m <- b
    print('assuming population at equilibrium')
  }

  par <- NA
  if(!is.null(submodel)){
    par <- submodel
  } else {
    par <- getProps(mod, omega, gamma, rho, sigma, epsilon)
  }

  epsilon <- par$epsilon
  omega <- par$omega
  gamma <- par$gamma
  rho <- par$rho
  sigma <- par$sigma
  whichBeta <- par$whichBeta
  whichGamma <- par$whichGamma
  whichSigma <- par$whichSigma

  if(!is.null(R0)) beta <- getBeta(R0, rho, epsilon, sigma, m, gamma, N, whichBeta)

  hasEquil <- TRUE

  beta1 <- 0
  beta2 <- 0
  gamma1 <- 0
  gamma2 <- 0
  sigma1 <- 0
  sigma2 <- 0

  if (whichBeta==1) beta1 <- beta
  if (whichBeta==2) beta2 <- beta
  if (whichGamma==1) gamma1 <- gamma
  if (whichGamma==2) gamma2 <- gamma
  if (whichSigma==1) sigma1 <- sigma
  if (whichSigma==2) sigma2 <- sigma

  Sstar <- (b+sigma1+sigma2+epsilon)*((b+sigma1)*rho+b+gamma1+gamma2) + sigma2*rho
  Sstar <- Sstar/((beta1+beta2)*(b+sigma1+sigma2+epsilon) - (b+sigma1)*beta1- sigma2*beta1)

  if (Sstar>N) {
    Sstar <- N
    hasEquil <- FALSE
  }

  Istar <- (N - Sstar)*(b+sigma2+sigma1+epsilon)*(b+omega)
  Istar <- Istar/((b+sigma1+sigma2+epsilon)*(gamma2+b+omega)+(sigma2+b+omega)*(beta1*Sstar+rho))

  if (Istar>N) {
    #Istar <- N
    hasEquil <- FALSE
  }

  Estar <- beta1*Sstar*Istar+rho*Istar
  Estar <- Estar/(b+sigma1+sigma2+epsilon)

  Rstar <- N - round(Istar) - round(Sstar) - round(Estar)

  return(data.frame(S=round(Sstar),E=round(Estar),I=round(Istar),R=Rstar,hasEquil))
}
