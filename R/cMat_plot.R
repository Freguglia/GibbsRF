#' Plots interaction structure
#'
#' @export
#' @author Victor Freguglia Souza

cMatPlot = function(cMat){
  x = complete_cMat(cMat)[,1]
  y = complete_cMat(cMat)[,2]
  df = data.frame(x,y)
  maxw = max(abs(df$x),abs(df$y)) + 0.5
  p = ggplot(df,aes(x=x,y=y,xmin=-maxw,xmax=maxw)) +
    geom_tile(fill="gray90",colour="black") +
    scale_x_continuous(limits=c(-maxw,maxw)) +
    scale_y_continuous(limits=c(-maxw,maxw)) +
    theme_void() +
    annotate("text", x = 0, y = 0, label = "i")
  plot(p)
}
