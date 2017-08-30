#' Estimation of potentials values via Stochastic Approximation.
#'
#' Estimates potentials vias Stochastic Approximation algorithm.
#' @param X The observed random field (matrix).
#' @param gModel A GibbsModel object (V and vMat can be NULL).
#' @param type Model type. Either "general", "symetric" or "equal". Check GibbsMPLE for details.
#' @param initial Wich method to get initial estimates. Currently only "MPLE" available.
#' @param MC_size Number of fields to sample at each step.
#' @param iter Maximum number of iterations.
#' @param macrosteps Number of macrosteps per random field simulation.
#' @return A GibbsModel object estimated potentials.
#' @author Victor Freguglia Souza
#' @examples
#' StocAp(example.X,example.GibbsModel,"symetric",MC_size = 5)
#' @export


StocAp = function(X,gModel,type,initial="MPLE",MC_size = 1,iter = 10000,macrosteps = 40){
  cMat = gModel$cMat
  G = gModel$G
  Tvec = vecStat(X,cMat,G,type)

  if(initial == "MPLE"){
    cat('\r','Getting MPLE for initialization')
    theta0 = GibbsMPLE(X,gModel,type) %>% vecTheta(type=type)
  }
  grad = 10
  i = 0
  theta = theta0
  StatList = matrix(0,nrow=MC_size,ncol=length(theta))
  while(i<iter&&max(abs(grad))>.0005){
    a = 10/(i/20+3)
    cat('\r','Iteration Number: ',i,' current theta: ',theta)
    gm = gModelTheta(theta,cMat,G,type)
    for(j in 1:MC_size){
      S = rGibbsRF(gm,macrosteps,NULL,dim(X))
      Svec = vecStat(S,cMat,G,type)
      StatList[j,] = Svec
    }
    E = colMeans(StatList)
    grad = Tvec - E
    n_theta = theta + a*grad
    theta = n_theta
    i = i+1
    if(i%%10==0){S %>% as.cimg %>% plot}
  }
  return(gm)
}

vecStat = function(X,cMat,G,type){
  stat = NULL
  stat[1:(G+1)] = table(X)/length(X)
  C = nrow(cMat)
  H = DifHistogram(X,cMat)
  H = H/rowSums(H)
  for(i in 1:C){
    if(type == "general"){
      stat = c(stat,H[i,])
    }
    if(type == "symetric"){
      vec = H[i,(G+1):(2*G+1)]
      vec[2:(G+1)] = vec[2:(G+1)] + H[i,G:1]
      stat = c(stat,vec)
    }
    if(type == "equal"){
      stat = c(stat,H[i,(G+1)],sum(H[i,-(G+1)]))
    }
  }
  return(stat)
}

vecTheta = function(gModel,type){
  theta = gModel$V
  vMat = gModel$vMat
  G = gModel$G
  C = nrow(vMat)
  for(i in 1:C){
    if(type == "general"){
      theta = c(theta,vMat[i,])
    }
    if(type == "symetric"){
      theta = c(theta,vMat[i,(G+1):(2*G+1)])
    }
    if(type == "equal"){
      theta = c(theta,vMat[i,(G+1):(G+2)])
    }
  }
  return(theta)
}

gModelTheta = function(theta,cMat,G,type){
  V = theta[1:(G+1)]
  vMat = matrix(0,nrow=nrow(cMat),ncol = (2*G+1))
  end = (G+1)
  for(i in 1:nrow(cMat)){
    if(type == "general"){
      vMat[i,] = theta[(end+1):(end+(2*G+1))]
      end = end + (2*G+1)
    }
    if(type == "symetric"){
      vMat[i,(G+1)] = theta[end+1]
      vec = theta[(end+2):(end+G+1)]
      vMat[i,(G+2):(2*G+1)] = vec
      vMat[i,1:G] = vec[length(vec):1]
      end = end + (G+1)
    }
    if(type == "equal"){
      vMat[i,(G+1)] = theta[end+1]
      vMat[i,-(G+1)] = theta[end+2]
      end = end+2
    }
  }
  return(GibbsModel(G,cMat,V,vMat))
}

#vMat = c(-2,-.5,2.5,-.5,-2,-2,-.5,2.5,-.5,-2,-.2,1,-.8,1,-.2) %>%
#matrix(byrow=TRUE,ncol=5)
#cMat = matrix(c(1,0,0,1,4,4),ncol=2,byrow=TRUE)
#gm = GibbsModel(2,cMat,c(0,0,0),vMat)
#X = rGibbsRF(gm)

