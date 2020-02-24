#' stochastic simulation function
#' @include RcppExports.R
#' @include getBeta.R
#' @include getProps.R
#' @export gSSA

gSSA <- function(N, R0=NULL, beta=NULL,
                 S0, I0, E0,
                 model=NULL,
                 gamma=0, sigma=0, rho=0, epsilon=0, omega=0,
                 b, m=NULL,
                 whichBeta=NULL, whichSigma=NULL, whichGamma=NULL,
                 submodel=NULL,
                 tmax=10, inc=52) {
  if(is.null(R0) & is.null(beta)) {
    print('either R0 or beta must be nonzero')
    return()
  }

  if(is.null(model) & is.null(submodel)) {
    print('must specify either a model ID or input a submodel object')
    return()
  }

  if(is.null(m)) {
    m <- b
    print('assuming population at equilibrium')
  }

  nevents <- 11

  par <- NA
  if(!is.null(submodel)){
    par <- submodel
  } else {
    par <- getProps(model, omega, gamma, rho, sigma, epsilon)
  }

  epsilon <- par$epsilon
  omega <- par$omega
  gamma <- par$gamma
  rho <- par$rho
  sigma <- par$sigma
  whichBeta <- par$whichBeta
  whichGamma <- par$whichGamma
  whichSigma <- par$whichSigma
  if(is.null(beta)) beta <- getBeta(R0, rho, epsilon, sigma, m, gamma, N, par$whichBeta)


  par <- c(beta, gamma, sigma, rho, epsilon, omega, b, m)
  spec <- as.integer(c(whichBeta,whichSigma,whichGamma))

  #simulate nsim times and plot together
  init <- as.integer(c(S0, E0, I0, max(N-S0-E0-I0,0), 0))
  sim <- simulate1(init, par, spec, nevents, tmax, inc)
  S <- sim[1:(tmax*inc)]
  E <- sim[(tmax*inc+1):(2*tmax*inc)]
  I <- sim[(2*tmax*inc+1):(3*tmax*inc)]
  R <- sim[(3*tmax*inc+1):(4*tmax*inc)]
  incidence <- sim[(4*tmax*inc+1):(5*tmax*inc)]

  N <- S+E+I+R
  return(data.frame(N=N,S=S,I=I,E=E,R=R,sero=E+I+R,incidence=incidence))
}
