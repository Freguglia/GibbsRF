#' Calculates penalized least squares for orthogonal design matrix, based on ordinary least squares estimates.
#' @param df.beta a data.frame representation of the OLS coefficients (result of buildFrequencyDf()).
#' @param type the type of penalty. "L0" corresponds to hard thresholding, "L1" corresponds to LASSO or "SCAD".
#' @param fnlambda a function that maps each pair of frequencies to a constant between 0 and 1, corresponding to how the penalty relates to n and m.
#' @param lambda a constant.
#' @value a data.frame with rows corresponding to combinations of frequencies and functions and their respective real-valued fourier transform coefficient.
#' @export
#' @author Victor Freguglia Souza
#'

betaPLS = function(df.beta,lambda,penalty.type="L0",fnlambda){

}
