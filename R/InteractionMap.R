#' Interaction Map
#'
#' Calculates distances between observed (normalized) GLDH and IRF expected probabilities.
#'
#' @author Victor Freguglia Souza
#' @export

InteractionMap = function(X,window_size=10){
  w = window_size
  G = max(X)
  direc = expand.grid(-w:w,-w:w) %>% as.matrix
  direc = direc[direc[,1]!=0 | direc[,2]!=0,]
  H = DifHistogram(X,direc)
  H = apply(H,1,function(x){return(x/sum(x))}) %>% t
  MP = -G:G
  MP = (1+G-abs(MP))/(1+G)^2
  MP = matrix(rep(MP,each=nrow(H)),nrow=nrow(H))
  distances = (MP-H)^2 %>% rowSums
  iMap = data.frame(cbind(direc,distances))
  names(iMap) = c("X","Y","Distance")
  return(iMap)
}
