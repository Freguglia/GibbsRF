#' Gray-level difference histogram
#'
#' Creates the histogram of gray-level difference for each relative position in a given field.
#' @param X the observed field (matrix object).
#' @param positions either the (C x 2) matrix of the relative positions considered or a GibbsModel object.
#' @return a matrix containing the difference histogram (counts from -G to G) for each relative position in each row.
#' @export
#' @author Victor Freguglia Souza

DifHistogram = function(X,positions){
  cMat=NULL
  if(class(positions)=="matrix"){cMat = positions;G=max(as.vector(X))}
  if(class(positions)=="GibbsModel"){cMat = positions$cMat;G = positions$G}
  if(is.null(positions)){stop("positions object doesn't have the required class (matrix or GibbsModel).")}
  return(DifHistogramcpp(X,cMat,G));
}
