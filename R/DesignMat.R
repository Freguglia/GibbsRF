#'Creates design Matrix for a selected set of fourier basis.
#'@param coeftab A filtered list of fourier basis. Result from betaPLS function.
#'@param dim Dimensions of the image.
#'@return a design matrix
#'@export
#'@author Victor Freguglia Souza
DesignMat = function(coeftab,dime){
  tblY = matrix(NA,nrow=dime[1],ncol=dime[2]) %>% as.cimg %>% cimg2df
  X = tblY
  N = dime[1]
  M = dime[2]
  d = nrow(coeftab)
  head(X)
  for(i in 1:d){
    fun1 = ifelse(coeftab$f1[i]=="cos",cos,sin)
    fun2 = ifelse(coeftab$f2[i]=="cos",cos,sin)
    Xcolumn = fun1(2*pi*X$x*coeftab$n[i]/N)*fun2(2*pi*X$y*coeftab$m[i]/M)
    Xcolumn = Xcolumn/(sum(Xcolumn^2))
    X[,(i+3)]= Xcolumn
  }
  D = X[,-c(1:3)] %>% as.matrix
  return(D)
}
