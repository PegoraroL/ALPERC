ALPERC_2D_viz_multiv<-function(D_add_j, LAMBDA, model, viz_factors, fixed_factors, fixed_factor_levels, what_plot){
  names_factors <- names(D_add_j)[-1]
  grid<-expand.grid(seq(0,1,length.out=100),
                    seq(0,1,length.out=100))
  names(grid)<-viz_factors
  empty_grid<-matrix(data=NA, nrow=nrow(grid), ncol=(length(names_factors)-length(names(grid))))
  for(i in 1:ncol(empty_grid)){
    empty_grid[,i]<-fixed_factor_levels[i]
  }
  empty_grid<-as_tibble(empty_grid, .name_repair = 'unique')
  names(empty_grid)<-fixed_factors


  pred_grid<-as_tibble(cbind(grid,empty_grid), .name_repair = 'unique')
  pred_grid<-pred_grid[names_factors]

  #several models!
  s_list<-list()
  s_norm_list<-list()
  data_obs_pred_test_list<-list()
  lista_z_norm<-list()
  indicatore<-0


  if(what_plot=="mean"){
    for(mdl in model){
      if(class(mdl)[1]=="hetGP"){
        y<-SE_pred_y<-NULL
        indicatore<-indicatore+1
        #this only works for hetgp models!
        X<-as.matrix(pred_grid)
        y<-predict(x = X, object = mdl)$mean

        #################################

        s<-NULL
        s<-interp(pred_grid%>%pull(viz_factors[1]), pred_grid%>%pull(viz_factors[2]), y, duplicate="mean")
        s_list[[indicatore]]<-s

        #non usiamo s_norm. Just in case dovessimo ampliare funzionalità funzione.
        s_norm<-NULL
        s_norm<-list(x=(s$x-min(s$x))/(max(s$x)-min(s$x)),y=(s$y-min(s$y))/(max(s$y)-min(s$y)),
                     z=(s$z-min(s$z))/(max(s$z)-min(s$z)))
        s_norm_list[[indicatore]]<-s_norm
        lista_z_norm[[indicatore]]<-s_norm$z
      }else if(mdl["method"]=="ranger"){
        y<-SE_pred_y<-NULL
        indicatore<-indicatore+1

        y<-predict(mdl$finalModel, pred_grid,type = "response")$predictions

        #################################

        s<-NULL
        s<-interp(pred_grid%>%pull(viz_factors[1]), pred_grid%>%pull(viz_factors[2]), y, dupl="mean")
        s_list[[indicatore]]<-s

        #non usiamo s_norm. Just in case dovessimo ampliare funzionalità funzione.
        s_norm<-NULL
        s_norm<-list(x=(s$x-min(s$x))/(max(s$x)-min(s$x)),y=(s$y-min(s$y))/(max(s$y)-min(s$y)),
                     z=(s$z-min(s$z))/(max(s$z)-min(s$z)))
        s_norm_list[[indicatore]]<-s_norm
        lista_z_norm[[indicatore]]<-s_norm$z
      }
    }
  }

  if(what_plot=="uncertainty"){
    for(mdl in model){
      if(class(mdl)[1]=="hetGP"){
        y<-SE_pred_y<-NULL
        indicatore<-indicatore+1
        #this only works for hetgp models!
        X<-as.matrix(pred_grid)
        y<-(sqrt(predict(x = X, object = mdl)$nugs +
                   predict(x = X, object = mdl)$sd2))
        SE_pred_y<-(sqrt(predict(x = X, object = mdl)$nugs +
                           predict(x = X, object = mdl)$sd2))

        pred_grid_y<-as_tibble(cbind(pred_grid,y, SE_pred_y), .name_repair = 'unique')
        #################################

        s<-NULL
        s<-interp(pred_grid%>%pull(viz_factors[1]), pred_grid%>%pull(viz_factors[2]), y, dupl="mean")
        s_list[[indicatore]]<-s

        #non usiamo s_norm. Just in case dovessimo ampliare funzionalità funzione.
        s_norm<-NULL
        s_norm<-list(x=(s$x-min(s$x))/(max(s$x)-min(s$x)),y=(s$y-min(s$y))/(max(s$y)-min(s$y)),
                     z=(s$z-min(s$z))/(max(s$z)-min(s$z)))
        s_norm_list[[indicatore]]<-s_norm
        lista_z_norm[[indicatore]]<-s_norm$z
      }else if(mdl["method"]=="ranger"){#this only works for ranger models trained via caret!
        y<-SE_pred_y<-NULL
        indicatore<-indicatore+1

        infjack<-NULL
        infjack<-predict(mdl$finalModel, pred_grid,type = "se",
                         se.method = "infjack")$se

        jack<-NULL
        jack<-predict(mdl$finalModel, pred_grid,type = "se",
                      se.method = "jack")$se

        y<-((infjack+jack)/2)
        SE_pred_y<-y

        pred_grid_y<-as_tibble(cbind(pred_grid,y, SE_pred_y), .name_repair = 'unique')
        #################################

        s<-NULL
        s<-interp(pred_grid%>%pull(viz_factors[1]), pred_grid%>%pull(viz_factors[2]), y, dupl="mean")
        s_list[[indicatore]]<-s

        #non usiamo s_norm. Just in case dovessimo ampliare funzionalità funzione.
        s_norm<-NULL
        s_norm<-list(x=(s$x-min(s$x))/(max(s$x)-min(s$x)),y=(s$y-min(s$y))/(max(s$y)-min(s$y)),
                     z=(s$z-min(s$z))/(max(s$z)-min(s$z)))
        s_norm_list[[indicatore]]<-s_norm
        lista_z_norm[[indicatore]]<-s_norm$z
      }
    }
  }

  grp_factors<- lapply(names_factors, as.symbol)
  grp_factors_viz<-lapply(viz_factors, as.symbol)
  runs_to_do<-D_add_j%>%group_by(!!!grp_factors)%>%count()%>%rename(n_repl=n)%>%ungroup()


  proposed_points_Xs<-NULL
  proposed_points_Xs<-LAMBDA%>%mutate(cluster=as.factor(cluster))
  proposed_points_afterclusters_Xs<-proposed_points_Xs


  df_runs_to_do<-df_prop_point<-NULL
  df_runs_to_do<-right_join(proposed_points_afterclusters_Xs, runs_to_do, by = names_factors)
  df_prop_point<-left_join(proposed_points_afterclusters_Xs, runs_to_do, by = names_factors)%>%select(-n_repl)

  #here, we consider the case in which we have from 2 to 10 responses.
  if(length(s_norm_list)>10){ message("ERROR: The plot is done for up to 10 responses.")
    stop()}
  if(length(s_norm_list)==2){
    s_norm_mean<-list(x=pmap_dbl(list(s_norm_list[[1]]$x,s_norm_list[[2]]$x), compose(partial(mean, na.rm = T), c)),
                      # mean(s1_norm$x,s2_norm$x,s3_norm$x,s4_norm$x),
                      y=pmap_dbl(list(s_norm_list[[1]]$y,s_norm_list[[2]]$y), compose(partial(mean, na.rm = T), c)),
                      # mean(s1_norm$y,s2_norm$y,s3_norm$y,s4_norm$y),
                      z=Reduce('+',lista_z_norm)/length(lista_z_norm)
                      # mean(s1_norm$z,s2_norm$z,s3_norm$z,s4_norm$z)
    )
  }
  if(length(s_norm_list)==3){
    s_norm_mean<-list(x=pmap_dbl(list(s_norm_list[[1]]$x,s_norm_list[[2]]$x,s_norm_list[[3]]$x), compose(partial(mean, na.rm = T), c)),
                      # mean(s1_norm$x,s2_norm$x,s3_norm$x,s4_norm$x),
                      y=pmap_dbl(list(s_norm_list[[1]]$y,s_norm_list[[2]]$y,s_norm_list[[3]]$y), compose(partial(mean, na.rm = T), c)),
                      # mean(s1_norm$y,s2_norm$y,s3_norm$y,s4_norm$y),
                      z=Reduce('+',lista_z_norm)/length(lista_z_norm)
                      # mean(s1_norm$z,s2_norm$z,s3_norm$z,s4_norm$z)
    )
  }
  if(length(s_norm_list)==4){
    s_norm_mean<-list(x=pmap_dbl(list(s_norm_list[[1]]$x,s_norm_list[[2]]$x,s_norm_list[[3]]$x,s_norm_list[[4]]$x), compose(partial(mean, na.rm = T), c)),
                      # mean(s1_norm$x,s2_norm$x,s3_norm$x,s4_norm$x),
                      y=pmap_dbl(list(s_norm_list[[1]]$y,s_norm_list[[2]]$y,s_norm_list[[3]]$y,s_norm_list[[4]]$y), compose(partial(mean, na.rm = T), c)),
                      # mean(s1_norm$y,s2_norm$y,s3_norm$y,s4_norm$y),
                      z=Reduce('+',lista_z_norm)/length(lista_z_norm)
                      # mean(s1_norm$z,s2_norm$z,s3_norm$z,s4_norm$z)
    )
  }
  if(length(s_norm_list)==5){
    s_norm_mean<-list(x=pmap_dbl(list(s_norm_list[[1]]$x,s_norm_list[[2]]$x,s_norm_list[[3]]$x,s_norm_list[[4]]$x,s_norm_list[[5]]$x), compose(partial(mean, na.rm = T), c)),
                      # mean(s1_norm$x,s2_norm$x,s3_norm$x,s4_norm$x),
                      y=pmap_dbl(list(s_norm_list[[1]]$y,s_norm_list[[2]]$y,s_norm_list[[3]]$y,s_norm_list[[4]]$y,s_norm_list[[5]]$y), compose(partial(mean, na.rm = T), c)),
                      # mean(s1_norm$y,s2_norm$y,s3_norm$y,s4_norm$y),
                      z=Reduce('+',lista_z_norm)/length(lista_z_norm)
                      # mean(s1_norm$z,s2_norm$z,s3_norm$z,s4_norm$z)
    )
  }
  if(length(s_norm_list)==6){
    s_norm_mean<-list(x=pmap_dbl(list(s_norm_list[[1]]$x,s_norm_list[[2]]$x,s_norm_list[[3]]$x,s_norm_list[[4]]$x,s_norm_list[[5]]$x,s_norm_list[[6]]$x), compose(partial(mean, na.rm = T), c)),
                      # mean(s1_norm$x,s2_norm$x,s3_norm$x,s4_norm$x),
                      y=pmap_dbl(list(s_norm_list[[1]]$y,s_norm_list[[2]]$y,s_norm_list[[3]]$y,s_norm_list[[4]]$y,s_norm_list[[5]]$y,s_norm_list[[6]]$y), compose(partial(mean, na.rm = T), c)),
                      # mean(s1_norm$y,s2_norm$y,s3_norm$y,s4_norm$y),
                      z=Reduce('+',lista_z_norm)/length(lista_z_norm)
                      # mean(s1_norm$z,s2_norm$z,s3_norm$z,s4_norm$z)
    )
  }
  if(length(s_norm_list)==7){
    s_norm_mean<-list(x=pmap_dbl(list(s_norm_list[[1]]$x,s_norm_list[[2]]$x,s_norm_list[[3]]$x,
                                      s_norm_list[[4]]$x,s_norm_list[[5]]$x,s_norm_list[[6]]$x,s_norm_list[[7]]$x), compose(partial(mean, na.rm = T), c)),
                      # mean(s1_norm$x,s2_norm$x,s3_norm$x,s4_norm$x),
                      y=pmap_dbl(list(s_norm_list[[1]]$y,s_norm_list[[2]]$y,s_norm_list[[3]]$y,
                                      s_norm_list[[4]]$y,s_norm_list[[5]]$y,s_norm_list[[6]]$y,s_norm_list[[7]]$y), compose(partial(mean, na.rm = T), c)),
                      # mean(s1_norm$y,s2_norm$y,s3_norm$y,s4_norm$y),
                      z=Reduce('+',lista_z_norm)/length(lista_z_norm)
                      # mean(s1_norm$z,s2_norm$z,s3_norm$z,s4_norm$z)
    )
  }
  if(length(s_norm_list)==8){
    s_norm_mean<-list(x=pmap_dbl(list(s_norm_list[[1]]$x,s_norm_list[[2]]$x,s_norm_list[[3]]$x,
                                      s_norm_list[[4]]$x,s_norm_list[[5]]$x,s_norm_list[[6]]$x,
                                      s_norm_list[[7]]$x,s_norm_list[[8]]$x), compose(partial(mean, na.rm = T), c)),
                      # mean(s1_norm$x,s2_norm$x,s3_norm$x,s4_norm$x),
                      y=pmap_dbl(list(s_norm_list[[1]]$y,s_norm_list[[2]]$y,s_norm_list[[3]]$y,
                                      s_norm_list[[4]]$y,s_norm_list[[5]]$y,s_norm_list[[6]]$y,
                                      s_norm_list[[7]]$y,s_norm_list[[8]]$y), compose(partial(mean, na.rm = T), c)),
                      # mean(s1_norm$y,s2_norm$y,s3_norm$y,s4_norm$y),
                      z=Reduce('+',lista_z_norm)/length(lista_z_norm)
                      # mean(s1_norm$z,s2_norm$z,s3_norm$z,s4_norm$z)
    )
  }
  if(length(s_norm_list)==9){
    s_norm_mean<-list(x=pmap_dbl(list(s_norm_list[[1]]$x,s_norm_list[[2]]$x,s_norm_list[[3]]$x,
                                      s_norm_list[[4]]$x,s_norm_list[[5]]$x,s_norm_list[[6]]$x,
                                      s_norm_list[[7]]$x,s_norm_list[[8]]$x,s_norm_list[[9]]$x), compose(partial(mean, na.rm = T), c)),
                      # mean(s1_norm$x,s2_norm$x,s3_norm$x,s4_norm$x),
                      y=pmap_dbl(list(s_norm_list[[1]]$y,s_norm_list[[2]]$y,s_norm_list[[3]]$y,
                                      s_norm_list[[4]]$y,s_norm_list[[5]]$y,s_norm_list[[6]]$y,
                                      s_norm_list[[7]]$y,s_norm_list[[8]]$y,s_norm_list[[9]]$y), compose(partial(mean, na.rm = T), c)),
                      # mean(s1_norm$y,s2_norm$y,s3_norm$y,s4_norm$y),
                      z=Reduce('+',lista_z_norm)/length(lista_z_norm)
                      # mean(s1_norm$z,s2_norm$z,s3_norm$z,s4_norm$z)
    )
  }
  if(length(s_norm_list)==10){
    s_norm_mean<-list(x=pmap_dbl(list(s_norm_list[[1]]$x,s_norm_list[[2]]$x,s_norm_list[[3]]$x,
                                      s_norm_list[[4]]$x,s_norm_list[[5]]$x,s_norm_list[[6]]$x,
                                      s_norm_list[[7]]$x,s_norm_list[[8]]$x,s_norm_list[[9]]$x,s_norm_list[[10]]$x), compose(partial(mean, na.rm = T), c)),
                      # mean(s1_norm$x,s2_norm$x,s3_norm$x,s4_norm$x),
                      y=pmap_dbl(list(s_norm_list[[1]]$y,s_norm_list[[2]]$y,s_norm_list[[3]]$y,
                                      s_norm_list[[4]]$y,s_norm_list[[5]]$y,s_norm_list[[6]]$y,
                                      s_norm_list[[7]]$y,s_norm_list[[8]]$y,s_norm_list[[9]]$y,s_norm_list[[10]]$y), compose(partial(mean, na.rm = T), c)),
                      # mean(s1_norm$y,s2_norm$y,s3_norm$y,s4_norm$y),
                      z=Reduce('+',lista_z_norm)/length(lista_z_norm)
                      # mean(s1_norm$z,s2_norm$z,s3_norm$z,s4_norm$z)
    )
  }

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
  image(s_norm_mean,col=topo.colors(128), main="", xlab = viz_factors[1], ylab=viz_factors[2],
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
