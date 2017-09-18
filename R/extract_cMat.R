#' Extract the cMat structure for a Thresholded Interaction Map
#'
#' @param iMap a (thresholded) Interaction Map.
#' @return a Matrix with the interacting pairs.
#'
#' @author Victor Freguglia Souza
#' @export

extract_cMat = function(iMap){
  cMat = iMap[iMap$Distance>0,1:2] %>% as.matrix
  cMat = cMat[(nrow(cMat)/2+1):(nrow(cMat)),]
  return(cMat)
}
