#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
NumericMatrix RandomMatrixCpp(IntegerVector dim, int max_value) {
  IntegerVector v = seq_len(max_value+1) - 1;
  int n = dim[0];
  int m = dim[1];
  NumericMatrix X(n,m);
  for(int i=0;i<n;i++){
    for(int j=0;j<m;j++){
      X(i,j) = sample(v,1)[0];
    }
  }
  return(X);
}

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
NumericVector ConditionalProbs(NumericMatrix X, IntegerVector position,
                                  int max_value, NumericMatrix cMat, NumericMatrix vMat,
                                  NumericVector V) {
  IntegerVector v = seq_len(max_value+1) - 1;
  int n = X.nrow();
  int m = X.ncol();
  int cons = cMat.nrow();
  int x = position[0] -1;
  int y = position[1] -1;
  int neix,neiy;
  NumericVector p(max_value + 1);
  IntegerVector vals = seq_len(max_value+1) - 1;
  double U;
  int dif;

  for(int value=0;value<=max_value;value++){
    U = V[value];
    for(int ne=0;ne<cons;ne++){
      neix = x + cMat(ne,0);
      neiy = y + cMat(ne,1);
      if(neix < n && neix >=0 && neiy < m && neiy>=0){
        dif = vals[value] - X(neix,neiy);
        U = U + vMat(ne,dif+max_value);
      }
    }
    p[value] = exp(U);
  }
  p = p/sum(p);
  return(p);
}

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
NumericMatrix rGRF(NumericMatrix cMat, NumericMatrix vMat, NumericVector V,
                   int max_value,int macrosteps,
                   NumericMatrix initial){
  int n,m;
  n = initial.nrow();
  m = initial.ncol();
  int pixnum = n*m;
  int x,y;
  IntegerVector values = seq_len(max_value+1) - 1;
  NumericVector cProbs(max_value + 1);
  IntegerVector coords(2);

  NumericMatrix X(n,m);
  X = Rcpp::clone(initial);

  IntegerVector runpath(pixnum);
  IntegerVector pixset = seq_len(pixnum) - 1;
  for(int ms=0;ms<macrosteps;ms++){
    runpath = sample(pixset,pixnum,false);
    for(int i=0;i<pixnum;i++){
      x = (runpath[i]/m) +1;
      y = (runpath[i]%m) +1;
      coords[0] = x;
      coords[1] = y;
      cProbs = ConditionalProbs(X,coords,max_value,cMat,vMat,V);
      X(x-1,y-1) = sample(values,1,false,cProbs)[0];
    }
  }

  return(X);
}

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
double log_plik(NumericMatrix X, NumericMatrix cMat, NumericMatrix vMat,
                NumericVector V, int max_value){
  int n,m,xval;
  n = X.nrow();
  m = X.ncol();
  IntegerVector coords(2);
  double loglik = 0.0;

  for(int x=0;x<n;x++){
    for(int y=0;y<m;y++){
      coords[0] = x+1;
      coords[1] = y+1;
      xval = X(x,y);
      loglik += log(ConditionalProbs(X,coords,max_value,cMat,vMat,V)[xval]);
    }
  }
  return(loglik);
}

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
NumericMatrix DifHistogramcpp(NumericMatrix X,NumericMatrix cMat,int G){
  int C = cMat.nrow();
  int n = X.nrow();
  int m = X.ncol();
  int neix,neiy,dife;
  NumericMatrix H(C,2*G+1);
  for(int x=0;x<n;x++){
    for(int y=0;y<m;y++){
      for(int c=0;c<C;c++){
        neix = x + cMat(0,c);
        neiy = y + cMat(1,c);
        if(neix>=0 && neix<n && neiy<m && neiy>=0){
          dife = X(x,y) - X(neix,neiy);
          H(c,(dife+G))++;
        }
      }
    }
  }
  return(H);
}
