\name{stepModel, tauleap, getRates, simRcpp}
\alias{stepModel}
\title{step model once}
\usage{
stepModel(state, r, spec, nevents)
tauleap(state, r, spec, nevents, tau)
getRates(state, par, nevents, t)
simulate1(state, par, spec, nevents, tmax, inc)
}
\description{
progress model with a single tau leap
increment states according to chosen tau
calculate rates based on current state and time
simulate based on tau leaping for a specified period of time
}
