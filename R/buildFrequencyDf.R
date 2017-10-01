#' Calculates coefficients for real-valued fourier transform of an image.
#' @param dfY a data.frame representation of an image.
#' @value a data.frame with rows corresponding to combinations of frequencies and functions and their respective real-valued fourier transform coefficient.
#' @export
#' @author Victor Freguglia Souza
buildFrequencyDf = function(dfY){
  N = max(dfY$x) + 1
  M = max(dfY$y) + 1
  x = 0:ceiling(max(dfY$x)/2)
  y = 0:ceiling(max(dfY$y)/2)
  func = c("cos","sin")
  coeftab = expand.grid(n=x,m=y,f1=func,f2=func) %>% tbl_df
  coeftab = coeftab %>%
    filter(!(n==0 & f1=="sin")) %>%
    filter(!(m==0 & f2=="sin")) %>%
    filter(!((2*n/N)%%1==0 & f1=="sin")) %>%
    filter(!((2*m/M)%%1==0 & f2=="sin"))
  betas = numeric(dim(coeftab)[1])
  nY = dfY
  nY$x = nY$x + 1
  nY$y = nY$y + 1
  matr = matrix(0, nrow=length(unique(nY$x)), ncol=length(unique(nY$y)))
  matr[cbind(nY$x, nY$y)] = nY$value
  ft = fft(matr) %>% as.matrix
  for(i in 1:(length(betas))){
    betas[i] = getFrequencyCoef2(ft,coeftab$n[i],
                                 coeftab$m[i],as.character(coeftab$f1[i]),
                                 as.character(coeftab$f2[i]))
    if(i%%1000==0) cat('\r',i)
  }
  coeftab = coeftab %>% mutate(coefs = betas)
  return(coeftab)
}
