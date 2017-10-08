#' Gets real-valued fourier transform coefficent for a specific pair of frequencies and functions based on the (matrix) fft.
#' @param ft the fourier transform (matrix) object of a image.
#' @param n,m the frequency to get coefficients for.
#' @param f1,f2 names of the functions for the coefficient, can be either "cos" or "sin".
#' @author Victor Freguglia Souza
#' @return the numeric value of the coefficient.
#' @export
getFrequencyCoef2 = function(ft,n,m,f1,f2){
  N = dim(ft)[1] + 1
  M = dim(ft)[2] + 1
  if(f1=="cos" & f2=="cos"){
    if(n>0){
      coef = (Re(ft[n+1,m+1]) + Re(ft[N-n,m+1]))/2
      return(coef)
    }
    if(m>0){
      coef = (Re(ft[n+1,m+1]) + Re(ft[n+1,M-m]))/2
      return(coef)
    }
    coef = Re(ft[1,1])
    return(coef)
  }

  if(f1=="sin" & f2=="sin" &(n*m)>0){
    coef = -(Re(ft[n+1,m+1]) - Re(ft[N-n,m+1]))/2
    return(coef)
  }

  if(f1=="cos" & f2=="sin" & m>0){
    if(n>0){
      coef = -(Im(ft[n+1,m+1])+ Im(ft[N-n,m+1]))/2
      return(coef)
    }
    coef = -Im(ft[n+1,m+1])
    return(coef)
  }

  if(f1=="sin" & f2=="cos" & n>0){
    if(m>0){
      coef = -(Im(ft[n+1,m+1])+ Im(ft[n+1,M-m]))/2
      return(coef)
    }
    coef = -Im(ft[n+1,m+1])
    return(coef)
  }
}

