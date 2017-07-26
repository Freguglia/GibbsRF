#' Gibbs Random Field simulation.
#'
#' Simulates a realization of a GibbsModel object via Gibbs Sampler.
#' @param gModel A GibbsModel object specifying wich model to sample from.
#' @param macrosteps The number of times each pixel will be simulated from it's conditional distribution.
#' @param initial A matrix indicating the initial object to start sampling from. If NULL, start with each pixel drawn from an independent discrete uniform distribution.
#' @param dim Used when the initial matrix is NULL. Indicates the dimension of the field to be sampled (x-axis first).
#' @return A matrix, realization from the specified Gibbs Model.
#' @author Victor Freguglia Souza
#' @examples
#' rGibbsRF(example.GibbsModel)
#' @export

rGibbsRF = function(gModel,macrosteps = 50,initial=NULL,dim=c(150,100)){
  #Check model
  if(class(gModel)!="GibbsModel"){stop(
    "modelSpec must be a GibbsModel object. Use GibbsModel() function to create it.")}
  cMat = gModel$cMat
  vMat = gModel$vMat
  V = gModel$V
  G = gModel$G
  if(is.null(initial)){X0 = RandomMatrixCpp(dim,G)}
  else{X0 = initial}
  cMat2 = complete_cMat(cMat)
  vMat2 = complete_vMat(vMat)

  return(rGRF(cMat2,vMat2,V,G,macrosteps,X0))
}
