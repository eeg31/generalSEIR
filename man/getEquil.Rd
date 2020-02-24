\name{getEquil}
\alias{getEquil}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
get equilibrium compartment populations
}
\description{
get equilibrium compartment populations from model specification or submodel object and transmission parameter (i.e., R0 or beta) 
}
\usage{
getEquil(mod=NULL, R0=NULL, beta=NULL,
                      omega=0, gamma=0, rho=0, sigma=0, epsilon=0,
                      submodel=NULL,
                      N=1, m=NULL, b)
}