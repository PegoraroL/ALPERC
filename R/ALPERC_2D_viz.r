ALPERC_2D_viz<-function(D_add_j, LAMBDA, model, viz_factors, fixed_factors, fixed_factor_levels, what_plot){
  names_factors <- names(D_add_j)[-1]
  grid<-expand.grid(seq(0,1,length.out=100),
                    seq(0,1,length.out=100))
  names(grid)<-viz_factors
  empty_grid<-matrix(data=NA, nrow=nrow(grid), ncol=(length(names_factors)-length(names(grid))))
  for(i in 1:ncol(empty_grid)){
    empty_grid[,i]<-fixed_factor_levels[i]
  }
  empty_grid<-as_tibble(empty_grid)
  names(empty_grid)<-fixed_factors


  pred_grid<-as_tibble(cbind(grid,empty_grid))
  pred_grid<-pred_grid[names_factors]

  if(class(model)[1]=="hetGP"){
    if(what_plot=="mean"){
      #this only works for hetgp models!
      X<-as.matrix(pred_grid)
      y<-predict(x = X, object = model)$mean
      #################################
    }

    if(what_plot=="uncertainty"){
      #this only works for hetgp models!
      X<-as.matrix(pred_grid)
      y<-(sqrt(predict(x = X, object = model)$nugs +
                 predict(x = X, object = model)$sd2))
      #################################
    }
  }else if(model["method"]=="ranger"){#this only works for ranger models trained via caret!
    if(what_plot=="mean"){
      y<-predict(model$finalModel, pred_grid,type = "response")$predictions
      #################################
    }

    if(what_plot=="uncertainty"){
      infjack<-NULL
      infjack<-predict(model$finalModel, pred_grid,type = "se",
                       se.method = "infjack")$se

      jack<-NULL
      jack<-predict(model$finalModel, pred_grid,type = "se",
                    se.method = "jack")$se

      y<-((infjack+jack)/2)
      #################################
    }
  }


  s<-NULL
  s<-interp(pred_grid%>%pull(viz_factors[1]), pred_grid%>%pull(viz_factors[2]), y, dupl="mean")

  #non usiamo s_norm. Just in case dovessimo ampliare funzionalità funzione.
  s_norm<-NULL
  s_norm<-list(x=(s$x-min(s$x))/(max(s$x)-min(s$x)),y=(s$y-min(s$y))/(max(s$y)-min(s$y)),
               z=(s$z-min(s$z))/(max(s$z)-min(s$z)))


  grp_factors<- lapply(names_factors, as.symbol)
  grp_factors_viz<-lapply(viz_factors, as.symbol)
  runs_to_do<-D_add_j%>%group_by(!!!grp_factors)%>%count()%>%rename(n_repl=n)%>%ungroup()


  proposed_points_Xs<-NULL
  proposed_points_Xs<-LAMBDA%>%mutate(cluster=as.factor(cluster))
  proposed_points_afterclusters_Xs<-proposed_points_Xs


  df_runs_to_do<-df_prop_point<-NULL
  df_runs_to_do<-right_join(proposed_points_afterclusters_Xs, runs_to_do, by = names_factors)
  df_prop_point<-left_join(proposed_points_afterclusters_Xs, runs_to_do, by = names_factors)%>%select(-n_repl)

  pal <- setNames(c25[1:nlevels(df_prop_point$cluster)], levels(df_prop_point$cluster))
  df_runs_to_do$colore<-pal[df_runs_to_do$cluster]
  df_runs_to_do$colore<-col2hex(df_runs_to_do$colore)
  df_runs_to_do<-df_runs_to_do%>%group_by(!!!grp_factors_viz)%>%mutate(n=sum(n_repl))

  df_prop_point$colore<-pal[df_prop_point$cluster]
  df_prop_point$colore<-col2hex(df_prop_point$colore)

  #numero di sfumature di blu per info su rank
  pal_blues<-colorRampPalette(brewer.pal(9,"Blues"))(nlevels(as.factor(proposed_points_afterclusters_Xs$rank)))
  pal_blues<-setNames(pal_blues[1:nlevels(as.factor(proposed_points_afterclusters_Xs$rank))], levels(as.factor(proposed_points_afterclusters_Xs$rank)))
  proposed_points_afterclusters_Xs$colore<-pal_blues[as.factor(proposed_points_afterclusters_Xs$rank)]


  sub_2d_plot_element<-c()
  for(i in seq_along(fixed_factors)){
    sub_tmp<-paste0(fixed_factors[i], "=",fixed_factor_levels[i])
    sub_2d_plot_element<-c(sub_2d_plot_element,sub_tmp)
  }
  sub_2d_plot<-paste(sub_2d_plot_element,collapse=", ")


  #2D plot
  # png(filename=paste0(plot.path,which_seed,"_",fnct_name,"_plot_2D_it",which_iteration,".png"),
  #     width = 600, height = 600)

  # plot.new() ## clean up device
  image(s,col=topo.colors(128), main="", xlab = viz_factors[1], ylab=viz_factors[2],
        sub=sub_2d_plot)
  points(proposed_points_afterclusters_Xs%>%select(all_of(viz_factors)),pch=as.numeric(as.character(proposed_points_afterclusters_Xs$cluster)),cex=2,
         col=proposed_points_afterclusters_Xs$colore)
  text(runs_to_do%>%ungroup()%>%select(all_of(viz_factors))%>%unique,
       labels=runs_to_do%>%ungroup()%>%select(all_of(viz_factors), n_repl)%>%group_by(!!!grp_factors_viz)%>%mutate(n=sum(n_repl))%>%
         select(-n_repl)%>%unique%>%pull(n),
       cex=.9, col="red", font=2)
  # dev.off()
  p <- recordPlot()
  return(p)
}
