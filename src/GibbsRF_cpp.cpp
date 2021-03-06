#include <Rcpp.h>
using namespace Rcpp;
using namespace R;

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
NumericVector ConditionalProbs(NumericMatrix X, IntegerVector position,  int C, NumericMatrix cMat, NumericMatrix vMat, NumericVector V) {
  
  int n = X.nrow(); int m = X.ncol();
  int n_neighbors = cMat.nrow();
  int x = position[0] -1; int y = position[1] -1;
  int neix, neiy;
  NumericVector p(C + 1);
  IntegerVector vals = seq_len(C+1) - 1;
  double U;
  int dif;

  for(int value = 0; value <= C; value++){
    U = V[value];
    for(int ne=0; ne < n_neighbors; ne++){
      neix = x + cMat(ne,0); neiy = y + cMat(ne,1);
      if(neix < n && neix >=0 && neiy < m && neiy>=0){
        dif = X(neix,neiy) - vals[value] ;
        U = U + vMat(ne,dif+C);
      }
    }
    p[value] = exp(U);
  }
  p = p/sum(p);
  return(p);
}

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
NumericMatrix rGRF(NumericMatrix cMat, NumericMatrix vMat, NumericVector V, int max_value,int macrosteps, NumericMatrix initial){
  
  int N = initial.nrow(); int M = initial.ncol();
  NumericMatrix X(N,M); X = Rcpp::clone(initial);
  int pixnum = N*M;
  int x,y;
  IntegerVector values = seq_len(max_value+1) - 1;
  NumericVector cProbs(max_value + 1);
  IntegerVector coords(2);

  
  IntegerVector runpath(pixnum);
  IntegerVector pixset = seq_len(pixnum) - 1;
  for(int ms=0;ms<macrosteps;ms++){
    runpath = sample(pixset,pixnum,false);
    for(int i=0;i<pixnum;i++){
      x = (runpath[i]/M) +1;
      y = (runpath[i]%M) +1;
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
        neix = x + cMat(c,0);
        neiy = y + cMat(c,1);
        if(neix>=0 && neix<n && neiy<m && neiy>=0){
          dife = X(x,y) - X(neix,neiy);
          H(c,(dife+G))++;
        }
      }
    }
  }
  return(H);
}

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
NumericMatrix MAPclassICM(NumericMatrix Y,NumericMatrix cMat, NumericVector V,
                       NumericMatrix vMat,int G, NumericVector mus, NumericVector sigmas,
                       NumericMatrix X, int iter){
  int N = Y.nrow();
  int M = Y.ncol();
  NumericMatrix labels = Rcpp::clone(X);
  NumericVector condProb(G+1);
  IntegerVector position(2);
  NumericVector cand(1);
  double max_value;
  int max_index=-1;

  for(int i=0;i<iter;i++){
    for(int x=0;x<N;x++){
      for(int y=0;y<M;y++){
        position[0] = x +1;
        position[1] = y +1;
        cand[0] = Y(x,y);
        condProb = ConditionalProbs(labels,position,G,cMat,vMat,V);
        max_value = 0;
        for(int k=0;k<(G+1);k++){
          condProb[k] = condProb[k]* dnorm(cand,mus[k],sigmas[k])[0];
          if(condProb[k]>max_value){
            max_value = condProb[k];
            max_index=k;
          }
        }
        labels(x,y) = max_index;
      }
    }
  }
  return(labels);
}


// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
NumericMatrix HMEM_CondProb(NumericMatrix Y,NumericMatrix X, NumericMatrix cMat,
                            NumericMatrix vMat, NumericVector V, int G,
                            double mu, double sigma, int candidate_value){
  NumericMatrix Px = Rcpp::clone(X);
  int N = X.nrow();
  int M = X.ncol();
  IntegerVector position(2);
  NumericVector cand(1);
  for(int x=0;x<N;x++){
    for(int y=0;y<M;y++){
      position[0] = x+1;
      position[1] = y+1;
      cand[0] = Y(x,y);
      Px(x,y) = ConditionalProbs(X,position,G,cMat,vMat,V)[candidate_value] *
        dnorm(cand,mu,sigma)[0];
    }
  }
  return(Px);
}


// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
NumericMatrix Hidden_CondSample(NumericMatrix Y,NumericMatrix X, NumericMatrix cMat,
                            NumericMatrix vMat, NumericVector V, int G,
                            NumericVector mu, NumericVector sigma){
  NumericMatrix Nx = Rcpp::clone(X);
  int N = X.nrow();
  int M = X.ncol();
  int x,y;
  IntegerVector values = seq_len(G+1) - 1;
  NumericVector probs(G+1);
  IntegerVector position(2);
  NumericVector cand(1);
  double media,dev;
  int pixnum = N*M;

  IntegerVector runpath(pixnum);
  IntegerVector pixset = seq_len(pixnum) - 1;
  runpath = sample(pixset,pixnum,false);
    for(int i=0;i<pixnum;i++){
      x = (runpath[i]/M) +1;
      y = (runpath[i]%M) +1;
      position[0] = x;
      position[1] = y;
      probs =   ConditionalProbs(Nx, position,G,cMat,vMat,V);
      for(int j=0;j<(G+1);j++){
        media = mu[j];
        dev=sigma[j];
        cand[0] = Y(x-1,y-1);
        probs[j] =  probs[j] * dnorm(cand,media,dev)[0];
      }
      probs = probs/sum(probs);
      Nx(x-1,y-1) = sample(values,1,false,probs)[0];
    }
  return(Nx);
}
