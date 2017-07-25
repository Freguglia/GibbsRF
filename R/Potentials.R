#' GibbsModel to DataFrame conversion.
#' Data Frame representation of GibbsModel object for better visualization.
#'
#' @param gModel A GibbsModel object with non-NULL potentials.
#' @return A data.frame object with potentials values for each difference (column) and each relative position (row).

Potentials = function(gModel){
  C = nrow(gModel$cMat)
  V = gModel$V
  vMat = gModel$vMat
  cMat = gModel$cMat
  if(is.null(gModel$V) || is.null(gModel$vMat)){
    stop("This GibbsModel object has no potentials specified")}
  positions = apply(gModel$cMat,1, function(x){
    return(paste("(",x[1],",",x[2],")",sep=""))})
  positions = c("Single",positions)
  Vrow = c(rep(NA,gModel$G),V)
  difs = seq(from=-gModel$G,gModel$G,by=1)
  difs = paste("V(",difs,")",sep="")
  df = data.frame(rbind(Vrow,vMat),row.names = positions)
  names(df) = c(difs)
  return(df)
}
