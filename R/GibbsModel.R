#' GibbsModel class object creation
#'
#' Creates a GibbsModel object, given specifications of a Gibbs Random Field model.
#' @param G An integer greater than 0, the maximum color value
#' @param cMat a (C x 2) matrix where each line corresponds to relative position of neighbours.
#' @param V a (G+1) vector corresponding to single site potentials (from 0 to G). Must sum 0. Can be NULL if vMat is also NULL.
#' @param vMat a (C x 2G+1) matrix specifying potential values for each possible difference (from -G to G) and each type of neighbour.
#' @return A GibbsModel object with the specified values.
#' @examples
#' GibbsModel(1,example.cMat,example.V,example.vMat)
#' @author Victor Freguglia Souza
#' @export
#' @useDynLib GibbsRF
#' @importFrom Rcpp evalCpp

GibbsModel = function(G,cMat,V=NULL,vMat=NULL){
  if((G%%1 != 0) || (G<=0)){stop("G must be a positive integer.")}
  if((ncol(cMat)!=2) || (!is.matrix(cMat))){stop(
    "cMat must be a (C x 2)-dimensional matrix.")}
  C = nrow(cMat)
  if(!is.null(vMat) && !is.null(V)){
    if((length(V) != (G+1)) || (!is.numeric(V))){stop(
      "V must be a numeric vector of length (G+1).")}
    if(sum(V)!=0){stop("V elements must have 0 sum.")}
    if((!is.matrix(vMat)) || (ncol(vMat) !=(2*G + 1)) || (nrow(vMat) != C)){stop(
      "vMat must be a (C x 2G+1) matrix.")}
    #if(any(abs(rowSums(vMat))>0.001)){stop("vMat rows must have 0 sum.")}
  }

  gm = list(G=G,cMat=cMat,V=V,vMat=vMat)
  class(gm) = "GibbsModel"
  return(gm)
}
