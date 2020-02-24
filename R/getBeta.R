#'a function to get Beta
#' @export getBeta

getBeta <- function(R0, rho, epsilon, sigma, m, gamma, N, whichBeta){
  beta <- numeric()
  if (whichBeta==1) beta <- R0*((epsilon+sigma+m)*(gamma+m+rho)-epsilon*rho)/(N*epsilon)
  if (whichBeta==2) beta <- R0*((epsilon+sigma+m)*(gamma+m+rho)-epsilon*rho)/(N*(epsilon+sigma+m))
  return(beta)
}
