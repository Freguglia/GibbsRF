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
#' GibbsMCMLE(example.X,example.GibbsModel,"symetric","MPLE") %>% Potentials
#' @export


GibbsMCMLE = function(X,gModelX,type="symetric",gModelRef="MPLE",nsamples=100,showProgress=TRUE,...){
  S = array(dim=c(nrow(X),ncol(X),nsamples))
  if(class(gModelX)!="GibbsModel"){stop("Please provide GibbsModel object in gModelX.")}
  if(class(gModelRef)=="character"){
    if(gModelRef=="MPLE"){phiModel = GibbsMPLE(X,gModelX,type,...)}
    if(gModelRef=="zero"){
      phiModel = gModelX
      phiModel$V = numeric(phiModel$G+1)
      phiModel$vMat = matrix(0,nrow=nrow(phiModel$cMat),ncol=2*phiModel$G+1)}
  }
  return(phiModel)
}
