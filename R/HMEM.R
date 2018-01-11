#' Segmentation of noisy image via Hidden Markov Field EM algorithm
#'
#' Provides EM estimation of parameters and segmentation.
#'
#' @author Victor Freguglia Souza
#' @export

HMEM = function(Y,gModel,classes=2,biasField=FALSE,filter_freqs = c(4,4),
                maxiter=150,mindif=10^-3,show_progress=TRUE,
                initial_mu,initial_sigma){
  #Perform initial parameter estimation and segmentation
  cMat = gModel$cMat
  vMat = gModel$vMat
  cMat2 = complete_cMat(cMat)
  vMat2 = complete_vMat(vMat)
  G = gModel$G
  V = gModel$V
  sigmas = initial_sigma
  mus = initial_mu
  iter = 0
  dif = 100

  bias = Y %>% LowPass(freqs = filter_freqs)
  if(!biasField){bias = bias*0}
  X = apply(Y-bias,c(1,2),function(x){
    ps = dnorm(x,mean=mus,sd = sigmas)
    return(which(ps==max(ps)))
  }) %>% as.matrix
  N = dim(X)[1]
  M = dim(X)[2]

  while(iter<maxiter && dif>mindif){
    #Estimate the class labels
    X =  MAPclassICM(Y-bias,cMat2,V,vMat2,G,mus,sigmas,X,5)


    #Calculate posterior distribution
    Px = array(dim=c(N,M,classes))
    for(k in 1:classes){
      Px[,,k] = HMEM_CondProb(Y-bias,X,cMat2,vMat2,V,G,mus[k],sigmas[k],k-1)
    }
    for(i in 1:N){
      for(j in 1:M){
        Px[i,j,] = Px[i,j,]/sum(Px[i,j,])
      }
    }


    #Update parameters
    n_mus = apply(matrix(1:classes),1,function(x){return(sum(Px[,,x]*(Y-bias))/sum(Px[,,x]))})
    n_sigmas = apply(matrix(1:classes),1,function(x){
      return(sum(Px[,,x]*(Y-bias-n_mus[x])^2)/sum(Px[,,x]))}) %>% sqrt

    #Estimate the Bias Field
    if(biasField){
      meanRes = matrix(0,N,M)
      psi = meanRes
      for(x in 1:N){
        for(y in 1:M){
          meanRes[x,y] = sum(Px[x,y,]*(Y[x,y]-n_mus)/sigmas^2)
          psi[x,y] = sum(Px[x,y,]/sigmas^2)
        }
      }
      bias = LowPass(meanRes,filter_freqs)/LowPass(psi,filter_freqs)
    }

    #check stopping conditions
    dif = max(max(abs(n_mus-mus)),max(abs(n_sigmas-sigmas)))
    iter = iter+1
    mus = n_mus
    sigmas = n_sigmas
    if(show_progress){
      cat("\r","Current mean: ",mus, sigmas, ". Iteration: ",iter,"            ")
      X %>% as.cimg %>% plot
    }
  }

  m = rbind(mus,sigmas) %>% data.frame
  colnames(m) = c(1:classes)
  rownames(m) = c("mean","sd")
  return(list(seg = X,theta=m,bias = bias))
}



#cMat = c(1,0,0,1,4,0,0,4) %>% matrix(ncol=2,byrow=T)
#vMat = c(-1,-1,1,-1,-1,-1,-1,1,-1,-1,.2,.2,-.2,.2,.2,.2,.2,-.2,.2,.2) %>% matrix(ncol=5,byrow=T)
#V = c(0,0,0)
#gm = GibbsModel(2,cMat,V,vMat)
#X = rGibbsRF(gm)
#mus = c(5,7.5,10)
#sds = c(1,.8,1)
#b = (X*0 + rnorm(length(X),sd=50)) %>% LowPass(freqs=c(2,2))
#b = b - (mean(b))
#Y = apply(X,c(1,2),function(x){return(rnorm(1,mean=mus[x+1],sd = sds[x+1]))})
#Z = Y+b
