#' Correlation Map
#'
#' Computes Correlation Map
#'
#' @param IMG An image in a matrix object.
#' @param window_size Size of the windows of relative positions com compute correlation.
#' @return A data.frame object with relative positions and absolute auto-correlation of each relative position
#' @author Victor Freguglia Souza
#' @export

CorMap = function(IMG,window_size = 10){
  N = nrow(IMG);M = ncol(IMG)
  positions = expand.grid(X = c(-window_size:window_size),Y = c(-window_size:window_size)) %>% tbl_df
  positions = positions[!((positions$X==0)&positions$Y==0),]
  w = window_size
  posCor = function(x,y){
    Xi = IMG[(w+1):(N-w),]
    Xi = Xi[,(w+1):(M-w)]
    Xic = IMG[(w+1+x):(N-w+x),]
    Xic = Xic[,(w+1+y):(M-w+y)]
    return(cor(as.vector(Xi),as.vector(Xic)))
  }
  positions = positions %>% rowwise() %>% mutate(Distance = abs(posCor(X,Y)))
  return(positions)
}
