#include <Rcpp.h>
#include <math.h>
using namespace Rcpp;

//[[Rcpp::export]]
IntegerVector stepModel(IntegerVector state, NumericVector r, IntegerVector spec, int nevents)
{
  IntegerVector var = clone(state);

  int* whichBeta; whichBeta = &spec[0];
  int* whichSigma; whichSigma = &spec[1];
  int* whichGamma; whichGamma = &spec[2];

  int* S; S = &var[0];
  int* E; E = &var[1];
  int* I; I = &var[2];
  int* R; R = &var[3];
  int* incidence; incidence = &var[4];

  double totalRate = sum(r);
  double x = R::runif(0,1);
  x *= totalRate;

  int mindex;
  bool isDone = FALSE;
  for (int i = 0; !isDone; ++i){
    double sum = 0;
    for (int j = 0; j <= i; ++j){
      sum += r[j];
    }
    if (sum > x) {
      mindex = i;
      isDone = TRUE;
    }
  }

  NumericVector nEvents(nevents);
  for (int i = 0; i < nevents; i++) {
    if (i == mindex) {
      nEvents[i] = 1;
    } else {
      nEvents[i] = 0;
    }
  }

  *incidence = *incidence + nEvents[7];
  *S = *S + nEvents[0] - nEvents[1] - nEvents[5] + nEvents[10];
  *E = *E - nEvents[2] - nEvents[6] - nEvents[7] + nEvents[9];
  *I = *I - nEvents[3] + nEvents[7] - nEvents[8] - nEvents[9];
  *R = *R - nEvents[4] - nEvents[10];
  if (*whichBeta == 1) {
    *E += nEvents[5];
  } else {
    *I += nEvents[5];
    *incidence += nEvents[5];
  }
  if (*whichSigma == 1) {
    *S += nEvents[6];
  } else {
    *R += nEvents[6];
  }
  if (*whichGamma == 1){
    *S += nEvents[8];
  } else {
    *R += nEvents[8];
  }

  return var;
}

//[[Rcpp::export]]
IntegerVector tauleap(IntegerVector state, NumericVector r, IntegerVector spec, int nevents, double tau){
  IntegerVector var = clone(state);

  int* whichBeta; whichBeta = &spec[0];
  int* whichSigma; whichSigma = &spec[1];
  int* whichGamma; whichGamma = &spec[2];

  int* S; S = &var[0];
  int* E; E = &var[1];
  int* I; I = &var[2];
  int* R; R = &var[3];
  int* incidence; incidence = &var[4];

  NumericVector nEvents(nevents);
  for (int i = 0; i < nevents; i++) {
    nEvents[i] = R::rpois(r[i] * tau);
  }

  *incidence = *incidence + nEvents[7];
  *S = *S + nEvents[0] - nEvents[1] - nEvents[5] + nEvents[10];
  *E = *E - nEvents[2] - nEvents[6] - nEvents[7] + nEvents[9];
  *I = *I - nEvents[3] + nEvents[7] - nEvents[8] - nEvents[9];
  *R = *R - nEvents[4] - nEvents[10];
  if (*whichBeta == 1) {
    *E += nEvents[5];
  } else {
    *I += nEvents[5];
    *incidence += nEvents[5];
  }
  if (*whichSigma == 1) {
    *S += nEvents[6];
  } else {
    *R += nEvents[6];
  }
  if (*whichGamma == 1){
    *S += nEvents[8];
  } else {
    *R += nEvents[8];
  }

  return var;
}

//[[Rcpp::export]]
NumericVector getRates(IntegerVector var, NumericVector par, int nevents)
{
  double* beta; beta = &par[0];
  double* gamma; gamma = &par[1];
  double* sigma; sigma = &par[2];
  double* rho; rho = &par[3];
  double* epsilon; epsilon = &par[4];
  double* omega; omega = &par[5];
  double* b; b = &par[6];
  double* m; m = &par[7];

  NumericVector rates(nevents);

  int N = sum(var);
  int* S; S = &var[0];
  int* E; E = &var[1];
  int* I; I = &var[2];
  int* R; R = &var[3];

  rates[0] = *b * N;
  rates[1] = *m * *S;
  rates[2] = *m * *E;
  rates[3] = *m * *I;
  rates[4] = *m * *R;
  rates[5] = *beta * *S * *I;
  rates[6] = *sigma * *E;
  rates[7] = *epsilon * *E;
  rates[8] = *gamma * *I;
  rates[9] = *rho * *I;
  rates[10] = *omega * *R;

  return rates;
}

//[[Rcpp::export]]
NumericVector simulate1(IntegerVector state, NumericVector par, IntegerVector spec, int nevents, int tmax, int inc)
{
  NumericVector out(state.size() * tmax * inc);
  out[0] = state[0];
  out[tmax * inc] = state[1];
  out[2 * tmax * inc] = state[2];
  out[3 * tmax * inc] = state[3];
  out[4 * tmax * inc] = state[4];

  int nt = 1;
  int minstate = 10;
  double t = 0.0;
  bool hasEI = true;
  bool allpos;
  bool thingsHappen = true;
  double tau;
  double maxtau = 1.0/365/2;
  double mintau = 1.0/365/1000;
  IntegerVector statenew(state.size());
  IntegerVector stateold(state.size());
  NumericVector rates(nevents);

  while (nt < tmax * inc && hasEI && thingsHappen){
    rates = getRates(state,par,nevents);

    if (sum(rates) == 0) {
      thingsHappen = false;
      //Rcout << "nothing happens";
      break;
    }

    tau = maxtau;

    if (state[0] < minstate){
      tau = 1/max(NumericVector::create(1/tau,rates[1],rates[5]));
    }
    if (state[1] < minstate){
      tau = 1/max(NumericVector::create(1/tau,rates[2],rates[6],rates[7]));
    }
    if (state[2] < minstate){
      tau = 1/max(NumericVector::create(1/tau,rates[3],rates[8],rates[9]));
    }
    if (state[3] < minstate){
      tau = 1/max(NumericVector::create(1/tau,rates[4],rates[10]));
    }

    bool stepped = false;
    while (!stepped && hasEI){
      statenew = tauleap(state,rates,spec,nevents,tau);

      allpos = true;
      for (int i=0; i < 5; i++){
        if (statenew[i] < 0) {
          allpos = false;
          ////Rcout << " stateold = ";
          ////Rcout << stateold;
          //Rcout << " ! negative value in new state ! ";
        }
      }

      if (allpos) {
        for (int i=0; i < 5; i++){
          stateold[i] = state[i];
          state[i] = statenew[i];
          ////Rcout << " tau leaped, t= ";
          ////Rcout << t;

        }
        stepped = true;
      } else {
        double sumR = sum(rates);
        double taunew = R::rexp(1.0/sumR);
        if ((taunew > maxtau) | ((tau/2.0) > mintau)) {
          tau = tau/2.0;
          //Rcout << " tau halved ";

        } else {
          stateold = state;
          ////Rcout << "; old state = ";
          ////Rcout << state;
          tau = taunew;
          state = stepModel(state,rates,spec,nevents);
          //Rcout << "stepped";
          ////Rcout << "; old state = ";
          ////Rcout << stateold;
          ////Rcout << "; rates = ";
          ////Rcout << rates;
          ////Rcout << ", tau = ";
          ////Rcout << tau;//print diagnostic info
          ////Rcout << ", t = ";
          ////Rcout << t;
          stepped = true;
        }
      }
      if (state[1]==0 && state[2]==0) {
        hasEI = false;
      }
    }

    t += tau;

    if(t * inc > nt){
      out[nt] = state[0];
      out[tmax * inc + nt] = state[1];
      out[2 * tmax * inc + nt] = state[2];
      out[3 * tmax * inc + nt] = state[3];
      out[4 * tmax * inc + nt] = state[4];

      state[4] = 0;
      //Rcout << state;
      ++nt;
    }
  }
  return out;
}

//[[Rcpp::export]]
IntegerVector simFixed(IntegerVector var, NumericVector par, IntegerVector spec, int nevents, int tmax)
{
  double t = 0;
  bool hasEI = TRUE;

  while (t < tmax && hasEI){
    NumericVector rates = getRates(var,par,nevents);

    double sumR = sum(rates);
    double tau = R::rexp(1/sumR);

    var = stepModel(var,rates,spec,nevents);

    t += tau;

    if (var[1]==0 && var[2]==0) {
      hasEI = FALSE;
    }
  }
  return var;
}
