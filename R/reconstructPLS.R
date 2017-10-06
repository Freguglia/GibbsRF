#' Reconstruct an image based on frequencies coefficients.
#'
#' @param beta_df the data.frame representation of frequency coefficients. Result of betaPLS() or buildFrequencyDf().
#' @param img_dim Dimensions of the image.
#' @author Victor Freguglia Souza
#' @export

reconstructPLS = function(beta_df,img_dim){
  N = img_dim[1]
  M = img_dim[2]
  dfY2 = expand.grid(x=c(0:(N-1)),y=c(0:(M-1))) %>% tbl_df %>% mutate(value = 0)

  for(i in 1:nrow(beta_df)){
    fun1 = ifelse(beta_df$f1[i]=="cos",cos,sin)
    fun2 = ifelse(beta_df$f2[i]=="cos",cos,sin)
    Xcolumn = fun1(2*pi*dfY2$x*beta_df$n[i]/N)*fun2(2*pi*dfY2$y*beta_df$m[i]/M)
    Xcolumn = Xcolumn/(sum(Xcolumn^2))
    dfY2$value = dfY2$value + Xcolumn*beta_df$coefs[i]
  }
  return(dfY2)
}
