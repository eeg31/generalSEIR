#' stochastic simulation function
#' @include RcppExports.R
#' @include getBeta.R
#' @include getEquil.R
#' @include getProps.R
#' @include models.rda
#' @export simSEIRV

simSEIRV <- function(N, R0=NULL, beta=NULL,
                 S0=NULL, I0=NULL, E0=NULL,
                 model=NULL,
                 gamma=0, sigma=0, rho=0, epsilon=0, omega=0,
                 v=0,
                 b, m=NULL,
                 whichBeta=NULL, whichSigma=NULL, whichGamma=NULL,
                 submodel=NULL,
                 tmax=10, inc=52, equilibrium=FALSE) {
  if(is.null(R0) & is.null(beta)) {
    print('either R0 or beta must be nonzero')
    return()
  }
  if(is.null(model) & is.null(submodel)) {
    print('must specify a model ID, set of specifications, or input a submodel object')
    return()
  }
  if(is.null(m)) {
    m <- b
    print('assuming population at equilibrium (m=b)')
  }
  par <- NA
  if(!is.null(submodel)){
    par <- submodel
  } else {
    par <- getProps(model, omega, gamma, rho, sigma, epsilon)
  }

  nevents <- as.integer(13)
  tmax <- as.integer(tmax)
  inc <- as.integer(inc)

  epsilon <- par$epsilon
  omega <- par$omega
  gamma <- par$gamma
  rho <- par$rho
  sigma <- par$sigma
  whichBeta <- par$whichBeta
  whichGamma <- par$whichGamma
  whichSigma <- par$whichSigma
  if(is.null(beta)) beta <- getBeta(R0, rho, epsilon, sigma, m, gamma, N, par$whichBeta)

  if(equilibrium | is.null(S0)){
    if(is.null(S0)) print('compartment sizes not given; starting at equilibrium proportions')
    state <- getEquil(mod=par$model, R0=R0, beta=beta,
                      omega=omega, gamma=gamma, rho=rho, sigma=sigma, epsilon=epsilon,
                      m=m, b)
    S0 <- max(5, S)
    E0 <- max(0, E)
    I0 <- max(1, I)
  }

  par <- as.numeric(c(beta, gamma, sigma, rho, epsilon, omega, b, m, v))
  spec <- as.integer(c(whichBeta,whichSigma,whichGamma))

  #simulate nsim times and plot together
  init <- as.integer(c(S0, E0, I0, max(N-S0-E0-I0,0), 0, 0))
  sim <- simulateVax(init, par, spec, nevents, tmax, inc)
  S <- sim[1:(tmax*inc)]
  E <- sim[(tmax*inc+1):(2*tmax*inc)]
  I <- sim[(2*tmax*inc+1):(3*tmax*inc)]
  R <- sim[(3*tmax*inc+1):(4*tmax*inc)]
  incidence <- sim[(4*tmax*inc+1):(5*tmax*inc)]
  V <- sim[(5*tmax*inc+1):(6*tmax*inc)]

  N <- S+E+I+R+V
  return(data.frame(N=N, S=S, I=I, E=E, R=R, sero=E+I+R, incidence=incidence, vaccinated=V))
}
