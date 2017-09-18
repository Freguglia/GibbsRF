#' Threshold Interaction Map
#'
#' Removes from Interaction Map all relative positions with lower Distance than a threshold based on k.
#'
#' @param iMap a DataFrame with an Interaction Map (see function InteractionMap() for details).
#' @param k the thresholding parameter. The cut value is set as MD + k*SD. MD being the mean distance in the map, and SD their sample deviation.
#' @return a DataFrame with the new Interaction Map with values thresholded.
#' @author Victor Freguglia Souza
#' @export


ThreshIMap = function(iMap,k){
  MD = mean(iMap$Distance)
  sigma = sd(iMap$Distance)
  iMap[iMap$Distance<(MD +k*sigma),3] = 0
  return(iMap)
}
