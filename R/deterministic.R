#' @include getBeta.R
#' @include getProps.R
#' @export deterministic

deterministic <- function(N, R0, beta=NULL,
                          S0, I0, E0,
                          gamma=0, sigma=0, rho=0, epsilon=0, omega=0,
                          b, m=NULL,
                          mod=NULL,
                          submodel=NULL,
                          tmax=10, inc=52) {

  if(is.null(R0) & is.null(beta)) {
    print('either R0 or beta must be nonzero')
    return()
  }

  if(is.null(mod) & is.null(submodel)) {
    print('must specify either model ID or input a submodel object')
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
  if(is.null(beta)) beta <- getBeta(R0, rho, epsilon, sigma, m, gamma, N, par$whichBeta)

  sim <- function(time, state, pars){
    with(as.list(c(state, pars)),{
      dS <- b*N - m*S + omega*R - (whichSigma-2)*sigma*E - (whichGamma-2)*gamma*I - beta*S*I
      dE <- -m*E - (whichBeta-2)*beta*S*I + rho*I - (sigma+epsilon)*E
      dI <- -m*I + (whichBeta-1)*beta*S*I + epsilon*E - gamma*I - rho*I
      dR <- -m*R + (whichSigma-1)*sigma*E + (whichGamma-1)*gamma*I - omega*R
      return(list(c(dS,dE,dI,dR)))
    })
  }

  par <- c(beta=beta, gamma=gamma, sigma=sigma, rho=rho, epsilon=epsilon, omega=omega, b=b, m=m,
           whichBeta=whichBeta, whichSigma=whichSigma, whichGamma=whichGamma)

  init <- c(S = S0,
            E = E0,
            I = I0,
            R = N-S0-E0-I0)

  out <- ode(init, seq(0,tmax-1/inc,1/inc), sim, par)

  return(out)#[,3] is E, [,4] is I, [,5] is R
}


