#' @export
low_pass <- function(cimg, ...){
   cimg %>% 
     cimg2df() %>%
     buildFrequencyDf() %>%
     betaPLS(., ...) %>%
    reconstructPLS(., img_dim = dim(cimg)) %>%
    df2cimg()     
}
