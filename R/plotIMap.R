#' Visualization for Interaction Maps
#'
#' Plots the interaction map.
#'
#' @param iMap a DataFrame object with an Interaction Map (result of the InteractionMap() function)
#' @author Victor Freguglia Souza
#' @export

plotIMap = function(iMap){
  n = max(iMap$X) + .5
  m = max(iMap$Y) + .5
  pl = ggplot(iMap,aes(x=X,y=Y,z=Distance)) +
    stat_summary_2d(colour="black",breaks = list(x=seq(-n,n,1),y=seq(-m,m,1))) +
    scale_fill_gradient(low="white", high="black") +
    theme_bw() +
    theme(axis.line = element_line(colour = "black"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          panel.background = element_blank())
  plot(pl)
}
