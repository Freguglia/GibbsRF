#' Segmentation of noisy image via Hidden Markov Field EM algorithm
#'
#' Provides EM estimation of parameters and segmentation.
#'
#' @author Victor Freguglia Souza
#' @export

HMEM = function(Y,gModel,classes=2,biasField=FALSE,filter_freqs = c(4,4),
                maxiter=150,mindif=10^-3,show_progress=TRUE){
  #Perform initial parameter estimation and segmentation
  cMat = gModel$cMat
  vMat = gModel$vMat
  V = gModel$V
  sigmas = rep(0.2,classes)
  quants = quantile(Y,probs = seq(0,1,1/classes)) %>% as.numeric
  X = apply(Y,c(1,2),function(x){
    for(i in 1:classes){
      if(x>=quants[i] && x<=quants[i+1]){return(i-1)}
    }
  }) %>% as.matrix
  mus = apply(matrix(c(1:classes)),1,function(x){return(mean(Y[X==(x-1)]))})
  N = dim(X)[1]
  M = dim(X)[2]
  iter = 0
  dif = 100
  while(iter<maxiter || dif>mindif){

    #Estimate the Bias Field
    if(biasField){

    }


    #Calculate the likelihood distribution
    g = dnorm(Y,mean=mus[X+1],sd = sds[X+1])

    #Estimate the class labels
    for(j in 1:5){
      for(x in 1:N){
        for(y in 1:M){
          condprob = log(dnorm(Y[x,y],mean=mus,sd = sigmas)) +
            log(ConditionalProbs(X,c(x,y),(classes-1),cMat,vMat,V))
          X[x,y] = which(condprob==max(condprob)) -1
        }
      }
    }
    #Calculate posterior distribution
    Px = array(dim=c(N,M,classes))
    for(x in 1:N){
      for(y in 1:M){
        Px[x,y,] = dnorm(Y[x,y],mean=mus[1:classes],sd=sigmas[1:classes]) *
          ConditionalProbs(X,c(x,y),(classes-1),cMat,vMat,V)
        Px[x,y,] =Px[x,y, ]/sum(Px[x,y,])
      }
    }

    #Update parameters
    n_mus = apply(matrix(1:classes),1,function(x){return(sum(Px[,,x]*Y)/sum(Px[,,x]))})
    n_sigmas = apply(matrix(1:classes),1,function(x){
      return(sum(Px[,,x]*(Y-mus[x])^2)/sum(Px[,,x]))}) %>% sqrt

    #check stopping conditions
    dif = max(c(abs(n_mus-mus),abs(n_sigmas-sigmas)))
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
  return(list(seg = X,theta=m))
}



#cMat = c(1,0,0,1,4,0,0,4) %>% matrix(ncol=2,byrow=T)
#vMat = c(-1,-1,1,-1,-1,-1,-1,1,-1,-1,.2,.2,-.2,.2,.2,.2,.2,-.2,.2,.2) %>% matrix(ncol=5,byrow=T)
#V = c(0,0,0)
#gm = GibbsModel(2,cMat,V,vMat)
#X = rGibbsRF(gm)
#mus = c(5,8,11)
#sds = c(1,.9,.8)
#Y = apply(X,c(1,2),function(x){return(rnorm(1,mean=mus[x+1],sd = sds[x+1]))})
