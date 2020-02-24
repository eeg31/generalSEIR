#' @include getBeta.R
#' @include models.rda
#' @export getProps
#' @export submodel

getProps <- function(model, omega, gamma, rho, sigma, epsilon){
  data("models")

  whichBeta <- models$whichBeta[model]
  whichGamma <- models$whichGamma[model]
  whichSigma <- models$whichSigma[model]
  if(models$omega[model]=="0") omega <- 0
  if(models$rho[model]=="0") rho <- 0
  if(models$rho[model]=="1") gamma <- 0
  if(models$whichSigma[model]==0) {
    sigma <- 0
    whichSigma <- 1
  }
  if(models$whichGamma[model]==0) {
    gamma <- 0
    whichGamma <- 1
  }
  if(rho==0 & whichBeta==2) epsilon <- 0

  return(list(model=model,
              omega=omega, gamma=gamma,
              rho=rho, sigma=sigma,
              epsilon=epsilon,
              whichSigma=whichSigma,
              whichGamma=whichGamma,
              whichBeta=whichBeta))
}

submodel <- function(model,
                     omega, gamma, rho, sigma, epsilon){
  return(getProps(model,
                  omega, gamma, rho, sigma, epsilon))
}
