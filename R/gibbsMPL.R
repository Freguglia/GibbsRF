#' Maximum Pseudo-Likelihood potentials estimate.
#'
#' Potentials estimate via Maximum Pseudo-Lilelihood.
#' @param X a matrix, realization of a GibbsModel.
#' @param gModel A GibbsModel object, vMat and cMat value should be NULL. If the object has non-NULL vMat or V values, they'll be ignored.
#' @param type Wich type of Energy function to be considered. "general" indicates all potentials are free. "symetric" indicates differences of the same absolute value have the same potential value. "equal" indicates all non-zero differences potentials have the same value.
#' @param ... arguments to be passed to the optim() function during maximization.
#' @return A GibbsModel object with the estimated vMat and V potentials values.
#' @author Victor Freguglia Souza
#' @examples
#' gibbsMPL(example.X,example.GibbsModel) %>% Potentials
#' #Check if close to real values example.GibbsModel %>% Potentials
#' @export

gibbsMPL = function(X,gModel,type = "symetric", ...){
  if(!type %in% c("symetric","equal","general")){
    stop("specified type must be symetric, equal or general.")}
  cMat = gModel$cMat
  G = gModel$G
  vecMPL = function(vecPars){
    gmobj = vec2gModel(vecPars,cMat,G,type)
    vMat = gmobj$vMat
    V = gmobj$V
    return(log_plik(X,cMat,vMat,V,G))
  }
  zeros = gModel2vec(gModel,type)*0

  optres = optim(zeros,vecMPL, ..., control = list(fnscale=-1))
  optgm = vec2gModel(optres$par,cMat,G,type)
  return(optgm)
}
