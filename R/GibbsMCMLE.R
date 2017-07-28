#' Monte-Carlo Maximum Likelihood Estimate of Potentials
#'
#' Provides MC Maximum Likelihood Estimation of Potentials for a given image.
#'
#' @author Victor Freguglia Souza
#' @param X the random field (matrix object).
#' @param gModelX a GibbsModel for X (V and vMat can be NULL).
#' @param type Type of function for Potentials. See GibbsMPLE() for reference.
#' @param gModelRef the reference model to simulate from for Monte-Carlo Likelihood function approximation. Either a GibbsModel object or a character with "zero" for an Independent Random Field, or "MPLE" for MPL estimated model.
#' @param nsamples an integer indicating the number of Monte-Carlo samples.
#' @param showProgress indicates whether the current step being done should be printed or not.
#' @param ... additional parameters to be used in the optimization function (optim()).
#' @return a GibbsModel object with the estimated potentials.
#' @examples
#' GibbsMCMLE(example.X,example.GibbsModel,"symetric","zero") %>% Potentials
#' @export


GibbsMCMLE = function(X,gModelX,type="symetric",gModelRef="zero",nsamples=10000,showProgress=TRUE,
                      macrosteps=70,...){
  if(class(gModelX)!="GibbsModel"){stop("Please provide GibbsModel object in gModelX.")}
  if(class(gModelRef)=="character"){
    if(gModelRef=="MPLE"){cat('\r', 'Calculating MPLE for reference model');
      phiModel = GibbsMPLE(X,gModelX,type,...)}
    if(gModelRef=="zero"){
      warning("Independent Random Field may cause numerical problems for large lattices.")
      phiModel = gModelX
      phiModel$V = numeric(phiModel$G+1)
      phiModel$vMat = matrix(0,nrow=nrow(phiModel$cMat),ncol=2*phiModel$G+1)}
  }
  n = nrow(X) ; m = ncol(X)
  G = gModelX$G ; cMat = gModelX$cMat
  if((phiModel$cMat != cMat) ||(phiModel$G != G)){stop("gModelX and gModelRef cMat and G must match.")}

  dX = DifHistogram(X,cMat)
  phiV = phiModel$V
  phivMat = phiModel$vMat

  dS = array(dim=c(nrow(dX),ncol(dX),nsamples))
  tS = matrix(0,nsamples,G+1)
  if(!(class(gModelRef)=="character" && gModelRef=="zero")){
    for(i in 1:nsamples){
      cat('\r', 'sampling field number',i,'                                ')
      S = rGibbsRF(phiModel,macrosteps,NULL,dim=c(n,m))
      dS[,,i] = DifHistogram(S,cMat)
      tS[i,] = table(S)
    }

    ll = function(params){
      gm = vec2gModel(params,cMat,G,type)
      V = gm$V
      vMat = gm$vMat
      ll = sum((vMat)*dX) + sum((V)*table(X))
      s = 1:nsamples %>% as.matrix
      logh = apply(s,1,function(x) return(sum(dS[,,x]*(vMat-phivMat)))) +
        apply(s,1,function(x) return(sum(tS[x,]*(V-phiV))))
      meio = mean(logh)
      dife = logh - meio
      second_term = meio*log(mean(exp(dife)))
      ll = ll - second_term
      return(ll)
    }
  }
  if(class(gModelRef)=="character" && gModelRef=="zero"){
    probs = apply(matrix(c(-G:G),nrow=1),1,function(x){(1+G-abs(x))/(1+G)^2})
    sizes = rowSums(dX)
    for(i in 1:nsamples){

    }
  }
  cat('\r','Starting optimization process. May take a while','                ')
  opt = optim(gModel2vec(phiModel,type),ll,control = list(fnscale=-1))
  optpar = opt$par
  cat('\r','Done','                                                     ')

  result = vec2gModel(optpar,cMat,G,type)

  return(result)
}
