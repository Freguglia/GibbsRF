#' Partial Correlation Map
#'
#' Computes conditional auto-correlations given a set of positions to condition in
#' @param IMG An image in a matrix object.
#' @param window_size Size of the windows to compute partial auto-correlations
#' @param included_list a data.frame with columns named \code{X} and \code{Y} and relative positions to condition on when computing the auto-correlations.
#' @return A data.frame with absolute conditional auto-correlations for each relative position.
#' @author Victor Freguglia Souza
#' @export

PartialCorMap = function(IMG,window_size=10,
                         included_list = data.frame(X = numeric(0),Y=numeric(0))){
  if(nrow(included_list)==0){return(CorMap(IMG,window_size))}
  else{
    N = nrow(IMG);M = ncol(IMG)
    positions = expand.grid(X = c(-window_size:window_size),Y = c(-window_size:window_size)) %>% tbl_df
    positions = positions[!((positions$X==0)&positions$Y==0),]
    w = window_size
    PartposCor = function(x,y){
      Xi = IMG[(w+1):(N-w),]
      Xi = Xi[,(w+1):(M-w)]
      Xi = as.vector(Xi)
      desi = matrix(nrow=length(Xi),ncol = nrow(included_list))
      for(i in 1:nrow(included_list)){
        a = included_list$X[i]
        b = included_list$Y[i]
        Xcol = IMG[(w+1+a):(N-w+a),]
        Xcol = Xcol[,(w+1+b):(M-w+b)]
        desi[,i] = as.vector(Xcol)
      }
      Xic = IMG[(w+1+x):(N-w+x),]
      Xic = Xic[,(w+1+y):(M-w+y)]
      Xp = lm(Xi~desi)$residuals
      return(cor(Xp,as.vector(Xic)))
    }
    positions = positions %>% rowwise() %>% mutate(Distance = abs(PartposCor(X,Y)))
    return(positions)
  }
}
