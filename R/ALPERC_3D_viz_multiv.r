ALPERC_3D_viz_multiv<-function(D_add_j, LAMBDA, model, viz_factors, fixed_factors, fixed_factor_levels, what_plot, n_pred_points=0){
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
  s_list<-z_list<-list()
  s_norm_list<-list()
  data_obs_pred_test_list<-list()
  lista_z_norm<-list()
  pred_grid_y_list<-y_list<-SE_pred_y_list<-list()
  indicatore<-0


  if(what_plot=="mean"){
    for(mdl in model){
      if(class(mdl)[1]=="hetGP"){
        y<-SE_pred_y<-NULL
        indicatore<-indicatore+1
        #this only works for hetgp models!
        X<-as.matrix(pred_grid)
        y<-predict(x = X, object = mdl)$mean
        SE_pred_y<-(sqrt(predict(x = X, object = mdl)$nugs +
                           predict(x = X, object = mdl)$sd2))

        pred_grid_y<-as_tibble(cbind(pred_grid,y, SE_pred_y), .name_repair = 'unique')
        pred_grid_y_list[[indicatore]]<-pred_grid_y
        y_list[[indicatore]]<-y
        SE_pred_y_list[[indicatore]]<-SE_pred_y
        #################################

        s<-NULL
        s<-interp(pred_grid%>%pull(viz_factors[1]), pred_grid%>%pull(viz_factors[2]), y, dupl="mean")
        s_list[[indicatore]]<-s
        z_list[[indicatore]]<-s$z

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

        infjack<-NULL
        infjack<-predict(mdl$finalModel, pred_grid,type = "se",
                         se.method = "infjack")$se

        jack<-NULL
        jack<-predict(mdl$finalModel, pred_grid,type = "se",
                      se.method = "jack")$se

        SE_pred_y<-((infjack+jack)/2)

        pred_grid_y<-as_tibble(cbind(pred_grid,y, SE_pred_y), .name_repair = 'unique')
        pred_grid_y_list[[indicatore]]<-pred_grid_y
        y_list[[indicatore]]<-y
        SE_pred_y_list[[indicatore]]<-SE_pred_y
        #################################

        s<-NULL
        s<-interp(pred_grid%>%pull(viz_factors[1]), pred_grid%>%pull(viz_factors[2]), y, duplicate="mean")
        s_list[[indicatore]]<-s
        z_list[[indicatore]]<-s$z

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
        pred_grid_y_list[[indicatore]]<-pred_grid_y
        y_list[[indicatore]]<-y
        SE_pred_y_list[[indicatore]]<-SE_pred_y
        #################################

        s<-NULL
        s<-interp(pred_grid%>%pull(viz_factors[1]), pred_grid%>%pull(viz_factors[2]), y, dupl="mean")
        s_list[[indicatore]]<-s
        z_list[[indicatore]]<-s$z

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
        pred_grid_y_list[[indicatore]]<-pred_grid_y
        y_list[[indicatore]]<-y
        SE_pred_y_list[[indicatore]]<-SE_pred_y
        #################################

        s<-NULL
        s<-interp(pred_grid%>%pull(viz_factors[1]), pred_grid%>%pull(viz_factors[2]), y, dupl="mean")
        s_list[[indicatore]]<-s
        z_list[[indicatore]]<-s$z

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

  #here, we consider the case in which we have from 2 to 10 responses. Normalized version
  # if(length(s_norm_list)>10){ message("ERROR: The plot is done for up to 10 responses.")
  #   stop()}
  # if(length(s_norm_list)==2){
  #   s_norm_mean<-list(x=pmap_dbl(list(s_norm_list[[1]]$x,s_norm_list[[2]]$x), compose(partial(mean, na.rm = T), c)),
  #                     # mean(s1_norm$x,s2_norm$x,s3_norm$x,s4_norm$x),
  #                     y=pmap_dbl(list(s_norm_list[[1]]$y,s_norm_list[[2]]$y), compose(partial(mean, na.rm = T), c)),
  #                     # mean(s1_norm$y,s2_norm$y,s3_norm$y,s4_norm$y),
  #                     z=Reduce('+',lista_z_norm)/length(lista_z_norm)
  #                     # mean(s1_norm$z,s2_norm$z,s3_norm$z,s4_norm$z)
  #   )
  # }
  # if(length(s_norm_list)==3){
  #   s_norm_mean<-list(x=pmap_dbl(list(s_norm_list[[1]]$x,s_norm_list[[2]]$x,s_norm_list[[3]]$x), compose(partial(mean, na.rm = T), c)),
  #                     # mean(s1_norm$x,s2_norm$x,s3_norm$x,s4_norm$x),
  #                     y=pmap_dbl(list(s_norm_list[[1]]$y,s_norm_list[[2]]$y,s_norm_list[[3]]$y), compose(partial(mean, na.rm = T), c)),
  #                     # mean(s1_norm$y,s2_norm$y,s3_norm$y,s4_norm$y),
  #                     z=Reduce('+',lista_z_norm)/length(lista_z_norm)
  #                     # mean(s1_norm$z,s2_norm$z,s3_norm$z,s4_norm$z)
  #   )
  # }
  # if(length(s_norm_list)==4){
  #   s_norm_mean<-list(x=pmap_dbl(list(s_norm_list[[1]]$x,s_norm_list[[2]]$x,s_norm_list[[3]]$x,s_norm_list[[4]]$x), compose(partial(mean, na.rm = T), c)),
  #                     # mean(s1_norm$x,s2_norm$x,s3_norm$x,s4_norm$x),
  #                     y=pmap_dbl(list(s_norm_list[[1]]$y,s_norm_list[[2]]$y,s_norm_list[[3]]$y,s_norm_list[[4]]$y), compose(partial(mean, na.rm = T), c)),
  #                     # mean(s1_norm$y,s2_norm$y,s3_norm$y,s4_norm$y),
  #                     z=Reduce('+',lista_z_norm)/length(lista_z_norm)
  #                     # mean(s1_norm$z,s2_norm$z,s3_norm$z,s4_norm$z)
  #   )
  # }
  # if(length(s_norm_list)==5){
  #   s_norm_mean<-list(x=pmap_dbl(list(s_norm_list[[1]]$x,s_norm_list[[2]]$x,s_norm_list[[3]]$x,s_norm_list[[4]]$x,s_norm_list[[5]]$x), compose(partial(mean, na.rm = T), c)),
  #                     # mean(s1_norm$x,s2_norm$x,s3_norm$x,s4_norm$x),
  #                     y=pmap_dbl(list(s_norm_list[[1]]$y,s_norm_list[[2]]$y,s_norm_list[[3]]$y,s_norm_list[[4]]$y,s_norm_list[[5]]$y), compose(partial(mean, na.rm = T), c)),
  #                     # mean(s1_norm$y,s2_norm$y,s3_norm$y,s4_norm$y),
  #                     z=Reduce('+',lista_z_norm)/length(lista_z_norm)
  #                     # mean(s1_norm$z,s2_norm$z,s3_norm$z,s4_norm$z)
  #   )
  # }
  # if(length(s_norm_list)==6){
  #   s_norm_mean<-list(x=pmap_dbl(list(s_norm_list[[1]]$x,s_norm_list[[2]]$x,s_norm_list[[3]]$x,s_norm_list[[4]]$x,s_norm_list[[5]]$x,s_norm_list[[6]]$x), compose(partial(mean, na.rm = T), c)),
  #                     # mean(s1_norm$x,s2_norm$x,s3_norm$x,s4_norm$x),
  #                     y=pmap_dbl(list(s_norm_list[[1]]$y,s_norm_list[[2]]$y,s_norm_list[[3]]$y,s_norm_list[[4]]$y,s_norm_list[[5]]$y,s_norm_list[[6]]$y), compose(partial(mean, na.rm = T), c)),
  #                     # mean(s1_norm$y,s2_norm$y,s3_norm$y,s4_norm$y),
  #                     z=Reduce('+',lista_z_norm)/length(lista_z_norm)
  #                     # mean(s1_norm$z,s2_norm$z,s3_norm$z,s4_norm$z)
  #   )
  # }
  # if(length(s_norm_list)==7){
  #   s_norm_mean<-list(x=pmap_dbl(list(s_norm_list[[1]]$x,s_norm_list[[2]]$x,s_norm_list[[3]]$x,
  #                                     s_norm_list[[4]]$x,s_norm_list[[5]]$x,s_norm_list[[6]]$x,s_norm_list[[7]]$x), compose(partial(mean, na.rm = T), c)),
  #                     # mean(s1_norm$x,s2_norm$x,s3_norm$x,s4_norm$x),
  #                     y=pmap_dbl(list(s_norm_list[[1]]$y,s_norm_list[[2]]$y,s_norm_list[[3]]$y,
  #                                     s_norm_list[[4]]$y,s_norm_list[[5]]$y,s_norm_list[[6]]$y,s_norm_list[[7]]$y), compose(partial(mean, na.rm = T), c)),
  #                     # mean(s1_norm$y,s2_norm$y,s3_norm$y,s4_norm$y),
  #                     z=Reduce('+',lista_z_norm)/length(lista_z_norm)
  #                     # mean(s1_norm$z,s2_norm$z,s3_norm$z,s4_norm$z)
  #   )
  # }
  # if(length(s_norm_list)==8){
  #   s_norm_mean<-list(x=pmap_dbl(list(s_norm_list[[1]]$x,s_norm_list[[2]]$x,s_norm_list[[3]]$x,
  #                                     s_norm_list[[4]]$x,s_norm_list[[5]]$x,s_norm_list[[6]]$x,
  #                                     s_norm_list[[7]]$x,s_norm_list[[8]]$x), compose(partial(mean, na.rm = T), c)),
  #                     # mean(s1_norm$x,s2_norm$x,s3_norm$x,s4_norm$x),
  #                     y=pmap_dbl(list(s_norm_list[[1]]$y,s_norm_list[[2]]$y,s_norm_list[[3]]$y,
  #                                     s_norm_list[[4]]$y,s_norm_list[[5]]$y,s_norm_list[[6]]$y,
  #                                     s_norm_list[[7]]$y,s_norm_list[[8]]$y), compose(partial(mean, na.rm = T), c)),
  #                     # mean(s1_norm$y,s2_norm$y,s3_norm$y,s4_norm$y),
  #                     z=Reduce('+',lista_z_norm)/length(lista_z_norm)
  #                     # mean(s1_norm$z,s2_norm$z,s3_norm$z,s4_norm$z)
  #   )
  # }
  # if(length(s_norm_list)==9){
  #   s_norm_mean<-list(x=pmap_dbl(list(s_norm_list[[1]]$x,s_norm_list[[2]]$x,s_norm_list[[3]]$x,
  #                                     s_norm_list[[4]]$x,s_norm_list[[5]]$x,s_norm_list[[6]]$x,
  #                                     s_norm_list[[7]]$x,s_norm_list[[8]]$x,s_norm_list[[9]]$x), compose(partial(mean, na.rm = T), c)),
  #                     # mean(s1_norm$x,s2_norm$x,s3_norm$x,s4_norm$x),
  #                     y=pmap_dbl(list(s_norm_list[[1]]$y,s_norm_list[[2]]$y,s_norm_list[[3]]$y,
  #                                     s_norm_list[[4]]$y,s_norm_list[[5]]$y,s_norm_list[[6]]$y,
  #                                     s_norm_list[[7]]$y,s_norm_list[[8]]$y,s_norm_list[[9]]$y), compose(partial(mean, na.rm = T), c)),
  #                     # mean(s1_norm$y,s2_norm$y,s3_norm$y,s4_norm$y),
  #                     z=Reduce('+',lista_z_norm)/length(lista_z_norm)
  #                     # mean(s1_norm$z,s2_norm$z,s3_norm$z,s4_norm$z)
  #   )
  # }
  # if(length(s_norm_list)==10){
  #   s_norm_mean<-list(x=pmap_dbl(list(s_norm_list[[1]]$x,s_norm_list[[2]]$x,s_norm_list[[3]]$x,
  #                                     s_norm_list[[4]]$x,s_norm_list[[5]]$x,s_norm_list[[6]]$x,
  #                                     s_norm_list[[7]]$x,s_norm_list[[8]]$x,s_norm_list[[9]]$x,s_norm_list[[10]]$x), compose(partial(mean, na.rm = T), c)),
  #                     # mean(s1_norm$x,s2_norm$x,s3_norm$x,s4_norm$x),
  #                     y=pmap_dbl(list(s_norm_list[[1]]$y,s_norm_list[[2]]$y,s_norm_list[[3]]$y,
  #                                     s_norm_list[[4]]$y,s_norm_list[[5]]$y,s_norm_list[[6]]$y,
  #                                     s_norm_list[[7]]$y,s_norm_list[[8]]$y,s_norm_list[[9]]$y,s_norm_list[[10]]$y), compose(partial(mean, na.rm = T), c)),
  #                     # mean(s1_norm$y,s2_norm$y,s3_norm$y,s4_norm$y),
  #                     z=Reduce('+',lista_z_norm)/length(lista_z_norm)
  #                     # mean(s1_norm$z,s2_norm$z,s3_norm$z,s4_norm$z)
  #   )
  # }

  # here, we consider the case in which we have from 2 to 10 responses. NON Normalized version
  if(length(s_list)>10){ message("ERROR: The plot is done for up to 10 responses.")
    stop()}
  if(length(s_list)==2){
    s_norm_mean<-list(x=pmap_dbl(list(s_list[[1]]$x,s_list[[2]]$x), compose(partial(mean, na.rm = T), c)),
                      # mean(s1_norm$x,s2_norm$x,s3_norm$x,s4_norm$x),
                      y=pmap_dbl(list(s_list[[1]]$y,s_list[[2]]$y), compose(partial(mean, na.rm = T), c)),
                      # mean(s1_norm$y,s2_norm$y,s3_norm$y,s4_norm$y),
                      z=Reduce('+',z_list)/length(z_list)
                      # mean(s1_norm$z,s2_norm$z,s3_norm$z,s4_norm$z)
    )
  }
  if(length(s_list)==3){
    s_norm_mean<-list(x=pmap_dbl(list(s_list[[1]]$x,s_list[[2]]$x,s_list[[3]]$x), compose(partial(mean, na.rm = T), c)),
                      # mean(s1_norm$x,s2_norm$x,s3_norm$x,s4_norm$x),
                      y=pmap_dbl(list(s_list[[1]]$y,s_list[[2]]$y,s_list[[3]]$y), compose(partial(mean, na.rm = T), c)),
                      # mean(s1_norm$y,s2_norm$y,s3_norm$y,s4_norm$y),
                      z=Reduce('+',z_list)/length(z_list)
                      # mean(s1_norm$z,s2_norm$z,s3_norm$z,s4_norm$z)
    )
  }
  if(length(s_list)==4){
    s_norm_mean<-list(x=pmap_dbl(list(s_list[[1]]$x,s_list[[2]]$x,s_list[[3]]$x,s_list[[4]]$x), compose(partial(mean, na.rm = T), c)),
                      # mean(s1_norm$x,s2_norm$x,s3_norm$x,s4_norm$x),
                      y=pmap_dbl(list(s_list[[1]]$y,s_list[[2]]$y,s_list[[3]]$y,s_list[[4]]$y), compose(partial(mean, na.rm = T), c)),
                      # mean(s1_norm$y,s2_norm$y,s3_norm$y,s4_norm$y),
                      z=Reduce('+',z_list)/length(z_list)
                      # mean(s1_norm$z,s2_norm$z,s3_norm$z,s4_norm$z)
    )
  }
  if(length(s_list)==5){
    s_norm_mean<-list(x=pmap_dbl(list(s_list[[1]]$x,s_list[[2]]$x,s_list[[3]]$x,s_list[[4]]$x,s_list[[5]]$x), compose(partial(mean, na.rm = T), c)),
                      # mean(s1_norm$x,s2_norm$x,s3_norm$x,s4_norm$x),
                      y=pmap_dbl(list(s_list[[1]]$y,s_list[[2]]$y,s_list[[3]]$y,s_list[[4]]$y,s_list[[5]]$y), compose(partial(mean, na.rm = T), c)),
                      # mean(s1_norm$y,s2_norm$y,s3_norm$y,s4_norm$y),
                      z=Reduce('+',z_list)/length(z_list)
                      # mean(s1_norm$z,s2_norm$z,s3_norm$z,s4_norm$z)
    )
  }
  if(length(s_list)==6){
    s_norm_mean<-list(x=pmap_dbl(list(s_list[[1]]$x,s_list[[2]]$x,s_list[[3]]$x,s_list[[4]]$x,s_list[[5]]$x,s_list[[6]]$x), compose(partial(mean, na.rm = T), c)),
                      # mean(s1_norm$x,s2_norm$x,s3_norm$x,s4_norm$x),
                      y=pmap_dbl(list(s_list[[1]]$y,s_list[[2]]$y,s_list[[3]]$y,s_list[[4]]$y,s_list[[5]]$y,s_list[[6]]$y), compose(partial(mean, na.rm = T), c)),
                      # mean(s1_norm$y,s2_norm$y,s3_norm$y,s4_norm$y),
                      z=Reduce('+',z_list)/length(z_list)
                      # mean(s1_norm$z,s2_norm$z,s3_norm$z,s4_norm$z)
    )
  }
  if(length(s_list)==7){
    s_norm_mean<-list(x=pmap_dbl(list(s_list[[1]]$x,s_list[[2]]$x,s_list[[3]]$x,
                                      s_list[[4]]$x,s_list[[5]]$x,s_list[[6]]$x,s_list[[7]]$x), compose(partial(mean, na.rm = T), c)),
                      # mean(s1_norm$x,s2_norm$x,s3_norm$x,s4_norm$x),
                      y=pmap_dbl(list(s_list[[1]]$y,s_list[[2]]$y,s_list[[3]]$y,
                                      s_list[[4]]$y,s_list[[5]]$y,s_list[[6]]$y,s_list[[7]]$y), compose(partial(mean, na.rm = T), c)),
                      # mean(s1_norm$y,s2_norm$y,s3_norm$y,s4_norm$y),
                      z=Reduce('+',z_list)/length(z_list)
                      # mean(s1_norm$z,s2_norm$z,s3_norm$z,s4_norm$z)
    )
  }
  if(length(s_list)==8){
    s_norm_mean<-list(x=pmap_dbl(list(s_list[[1]]$x,s_list[[2]]$x,s_list[[3]]$x,
                                      s_list[[4]]$x,s_list[[5]]$x,s_list[[6]]$x,
                                      s_list[[7]]$x,s_list[[8]]$x), compose(partial(mean, na.rm = T), c)),
                      # mean(s1_norm$x,s2_norm$x,s3_norm$x,s4_norm$x),
                      y=pmap_dbl(list(s_list[[1]]$y,s_list[[2]]$y,s_list[[3]]$y,
                                      s_list[[4]]$y,s_list[[5]]$y,s_list[[6]]$y,
                                      s_list[[7]]$y,s_list[[8]]$y), compose(partial(mean, na.rm = T), c)),
                      # mean(s1_norm$y,s2_norm$y,s3_norm$y,s4_norm$y),
                      z=Reduce('+',z_list)/length(z_list)
                      # mean(s1_norm$z,s2_norm$z,s3_norm$z,s4_norm$z)
    )
  }
  if(length(s_list)==9){
    s_norm_mean<-list(x=pmap_dbl(list(s_list[[1]]$x,s_list[[2]]$x,s_list[[3]]$x,
                                      s_list[[4]]$x,s_list[[5]]$x,s_list[[6]]$x,
                                      s_list[[7]]$x,s_list[[8]]$x,s_list[[9]]$x), compose(partial(mean, na.rm = T), c)),
                      # mean(s1_norm$x,s2_norm$x,s3_norm$x,s4_norm$x),
                      y=pmap_dbl(list(s_list[[1]]$y,s_list[[2]]$y,s_list[[3]]$y,
                                      s_list[[4]]$y,s_list[[5]]$y,s_list[[6]]$y,
                                      s_list[[7]]$y,s_list[[8]]$y,s_list[[9]]$y), compose(partial(mean, na.rm = T), c)),
                      # mean(s1_norm$y,s2_norm$y,s3_norm$y,s4_norm$y),
                      z=Reduce('+',z_list)/length(z_list)
                      # mean(s1_norm$z,s2_norm$z,s3_norm$z,s4_norm$z)
    )
  }
  if(length(s_list)==10){
    s_norm_mean<-list(x=pmap_dbl(list(s_list[[1]]$x,s_list[[2]]$x,s_list[[3]]$x,
                                      s_list[[4]]$x,s_list[[5]]$x,s_list[[6]]$x,
                                      s_list[[7]]$x,s_list[[8]]$x,s_list[[9]]$x,s_list[[10]]$x), compose(partial(mean, na.rm = T), c)),
                      # mean(s1_norm$x,s2_norm$x,s3_norm$x,s4_norm$x),
                      y=pmap_dbl(list(s_list[[1]]$y,s_list[[2]]$y,s_list[[3]]$y,
                                      s_list[[4]]$y,s_list[[5]]$y,s_list[[6]]$y,
                                      s_list[[7]]$y,s_list[[8]]$y,s_list[[9]]$y,s_list[[10]]$y), compose(partial(mean, na.rm = T), c)),
                      # mean(s1_norm$y,s2_norm$y,s3_norm$y,s4_norm$y),
                      z=Reduce('+',z_list)/length(z_list)
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

  #compute means of all ys and stderrs
  y_mean<-pmap_dbl(y_list, compose(partial(mean, na.rm = T), c))
  SE_pred_y_mean<-pmap_dbl(SE_pred_y_list, compose(partial(mean, na.rm = T), c))
  pred_grid_y_mean<-pred_grid_y%>%mutate(y=y_mean, SE_pred_y=SE_pred_y_mean)

  pred_grid_y_sample<-pred_grid_y_mean[sample(1:nrow(pred_grid_y_mean),n_pred_points),]
  # #3D PLOT
  plot_3d<-plot_ly(x=s_norm_mean$x,y=s_norm_mean$y,z=s_norm_mean$z, showscale=F) %>%
    add_surface(contours = list(
      z = list(
        show=TRUE,
        usecolormap=TRUE,
        # highlightcolor="#ff0000",
        project=list(z=TRUE)
      )
    ))%>%
    add_trace(data = pred_grid_y_sample, x = pred_grid_y_sample%>%pull(viz_factors[2]), y = pred_grid_y_sample%>%pull(viz_factors[1]), z = pred_grid_y_sample$y, mode = "markers", type = "scatter3d",
              opacity = 0.3,
              error_z = ~list(array = pred_grid_y_sample%>%pull(SE_pred_y),
                              color = 'red'),
              marker = list(size = 3, color = "red", symbol = 104)
    )%>%
    add_trace(data = df_runs_to_do, x = df_runs_to_do%>%pull(viz_factors[2]), y = df_runs_to_do%>%pull(viz_factors[1]), z = 0, mode = "markers+text", type = "scatter3d",
              text = ~n,
              marker = list(color=~colore,size = 5))%>%
    add_trace(data = df_prop_point, x = df_prop_point%>%pull(viz_factors[2]), y = df_prop_point%>%pull(viz_factors[1]), z = 0, mode = "markers", type = "scatter3d",
              marker = list(color=~colore,size = 3))%>%
    layout(title = "",
           # scene=list(zaxis="Y"),
           annotations= list(x = 1, y = 0, text = paste0("x: ", viz_factors[2], ", y: ",viz_factors[1], "; <br>", sub_2d_plot),
                             showarrow=FALSE),
           showlegend = FALSE)
  return(plot_3d)
}
