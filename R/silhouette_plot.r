silhouette_plot<-function(ALPERC_results){
  ALPERC_results$clust_choice%>%ggplot(aes(x=clusters, y=y, group = 1))+
    geom_point()+
    geom_line()+
    geom_vline(aes(xintercept=ALPERC_results$clust_choice%>%filter(y==max(y))%>%pull(clusters)), linetype = "dashed")+
    ylab("Average silhouette width")+
    xlab("Number of clusters")+
    theme_bw()
}