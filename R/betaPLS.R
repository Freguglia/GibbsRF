#' Calculates penalized least squares estimates for orthogonal design matrix, based on ordinary least squares estimates.
#' @param df.beta a data.frame representation of the OLS coefficients (result of buildFrequencyDf()).
#' @param type the type of penalty. "L0" corresponds to hard thresholding, "L1" corresponds to LASSO or "SCAD".
#' @param fnlambda a function that maps each pair of frequencies to a constant between 0 and 1, corresponding to how the penalty relates to n and m.
#' @param lambda a constant.
#' @return a data.frame with rows corresponding to combinations of frequencies and functions and their respective real-valued fourier transform coefficient.
#' @export
#' @author Victor Freguglia Souza
#'

betaPLS = function(coeftab,lambda,penalty.type="L0",fnlambda = function(b1,b2) exp(pmax(b1,b2)/2),a=3.7){
  coeftab = coeftab  %>% mutate(lambdas = lambda*fnlambda(n,m))
  if(penalty.type=="L0"){
    coeftab = coeftab %>% mutate(coefs = coefs*(abs(coefs)>lambdas)) %>%
      filter(abs(coefs)>0)
    return(coeftab)
  }

  if(penalty.type=="L1"){
    coeftab = coeftab %>% mutate(coefs = coefs*pmax(0,(1-lambdas/abs(coefs)))) %>%
      filter(abs(coefs)>0)
    return(coeftab)
  }

  if(penalty.type=="SCAD"){
    SCAD = function(b,lambda){
      num = length(b)

      ifelse(abs(b)<=(2*lambda),
             pmax(0,(abs(b) - lambda))*sign(b),
             ifelse(
               (((2*lambda) < abs(b)) && (abs(b) <= a*lambda )),
               ((a-1)*b - sign(b)*a*lambda)/(a-2) ,
                    b))
    }
    SCAD2 = Vectorize(SCAD)
    coeftab = coeftab %>% rowwise() %>% mutate(coefs = SCAD2(coefs,lambdas)) %>%
      filter(abs(coefs)>0)
    return(coeftab)

  }
  message("Penalty type not valid.")
}
