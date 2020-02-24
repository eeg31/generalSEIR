\name{getBeta}
\alias{getBeta}
\title{calculates beta (transmission rate) from R0 and other infection parameters}
\usage{
getBeta(R0, rho, epsilon, sigma, m, gamma, N, whichBeta)
}
\description{
calculates beta (transmission rate) from R0 and other infection parameters, model specification (where whichBeta=1 indicates S to E transmission and whichBeta=2 indicates S to I transmission), population size (N), and non-infection-related mortality rate (m)
}
