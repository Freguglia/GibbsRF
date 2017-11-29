#' Interaction structure selection based on partial auto-correlations
#'
#' @param IMG An image in a matrix object.
#' @param window_size Size of the window for candidates to be included in the interaction structure.
#' @param thr A thresholding parameter. When all absolute conditional auto-correlations not included in the model are belo \code{thr}, the selection stops.
#' @author Victor Freguglia Souza
#' @return A matrix with selected relative positions, ready to be used as \code{cMat} of GibbsModel.
#' @export

PartialCorSelection = function(IMG,window_size = 10,thr = 0.1){
  iMap1 = PartialCorMap(IMG,window_size)
  topcor = max(iMap1$Distance)
  C = NULL
  while(topcor>thr){
    adic = iMap1[iMap1$Distance==max(iMap1$Distance),(1:2)]
    if(is.null(C)){C = rbind(adic,-adic)}
    else{C = rbind(C,adic,-adic)}
    iMap1 = PartialCorMap(IMG,window_size,included_list = C)
    topcor = max(iMap1$Distance)
  }
  C = as.matrix(C)
  C = C[(1:nrow(C))%%2==1,]
  return(C)
}
