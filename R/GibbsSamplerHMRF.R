#' Gibbs Sampler for Hidden Markov Random Field model.
#' @param Y The observed values (matrix).
#' @param gModel a GibbsModel object.
#' @param imu a vector indicating initial mean values.
#' @param isigma a vector indicating initial SD values.
#' @param iseg initial segmentation.
#' @param X a matrix with n (number of observations) rows and p (number of fourier basis functions) columns. Use DeisgnMat() with the result of betaPLS() to obtain this matrix.
#' @export

GibbsSamplerHMRF = function(Y,gModel,imu,isigma2,iseg,coeftab,size = 10000){
  if(!is.matrix(Y)) stop("Y must be a matrix of integers.")

  dfY = Y %>% as.cimg %>% cimg2df
  nc = (gModel$G) + 1 #The number of classes
  mus = matrix(NA,nrow=size,ncol=nc) #Store the sampled means
  sigma2s = mus
  current_mu = imu
  current_sigma2 = isigma2
  L = iseg

  for(i in 1:size){
    cat('\r',i)
    #Update betas

    #Update mus and sigmas
    for(k in 1:nc){
      mus[i,k] = rnorm(1,mean = mean(Y[L==(k-1)]),sd = sqrt(current_sigma2[k]/sum(L==(k-1))))

      sigma2s[i,k] = rinvgamma(1, ( sum( L==(k-1)) -2 )/2, sum((Y[L==(k-1)] - mus[i,k])^2)/2 )
    }
    current_mu = mus[i,]
    current_sigma2 = sigma2s[i,]

    #Update segmentation
    L = Hidden_CondSample(Y,iseg,gModel$cMat,gModel$vMat,gModel$V,gModel$G,current_mu,sqrt(current_sigma2))
  }
  return(list(mus,sigma2s,L))
}
