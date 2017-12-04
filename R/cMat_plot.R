#' Plots interaction structure
#'
#' @export
#' @author Victor Freguglia Souza

cMatPlot = function(cMat){
  x = complete_cMat(cMat)[,1]
  y = complete_cMat(cMat)[,2]
  df = data.frame(x,y)
  p = ggplot(df,aes(x=x,y=y,xmin=-5,xmax=5)) +
    geom_tile(fill="gray90",colour="black") +
    scale_x_continuous(limits=c(-5,5)) +
    scale_y_continuous(limits=c(-5,5)) +
    theme_void() +
    annotate("text", x = 0, y = 0, label = "i")
  plot(p)
}
