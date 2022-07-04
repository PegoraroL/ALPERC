ALPERC<-function(n_add,D_cand,sigma_cand,S=NULL,strategy="exploration",varimp_distance,n_clust=NULL,n_boot=1000,force_n_boot=TRUE, alpha_rank=0.1, seed_rank=2105, paral=FALSE){
  #check on ID column

  if(missing(n_add)){
    stop("Argument \"n_add\" is missing, with no default")
  }
  if(missing(D_cand)){
    stop("Argument \"D_cand\" is missing, with no default")
  }
  if(missing(sigma_cand)){
    stop("Argument \"sigma_cand\" is missing, with no default")
  }
  if(missing(strategy)){
    stop("Argument \"strategy\" is missing, with no default")
  }
  if(missing(varimp_distance)){
    stop("Argument \"varimp_distance\" is missing, with no default")
  }
  if(missing(n_boot)){
    stop("Argument \"n_boot\" is missing, with no default")
  }
  if(missing(alpha_rank)){
    stop("Argument \"alpha_rank\" is missing, with no default")
  }
  if(missing(S)){
    stop("Argument \"S\" is missing, with no default. If you don't want to supply S, set S=NULL")
  }
  if(is.null(seed_rank)){
    message("ERROR: argument \"seed_rank\" set to NULL.")
    stop()
  }
  if(names(D_cand)[1]!=names(sigma_cand)[1]){
    stop("The first column of D_cand and sigma_cand must have same name (ID column)")
  }
  #check on strategy name
  if(!strategy %in% c("exploration" , "exploitation")){
    stop("Strategy must be \"exploration\" or \"exploitation\" ")
  }
  if(is.null(S)){
    varimp_distance=FALSE
  }else{
    #check if column names of  D_cand, == first row of S
    if(any(names(D_cand)[-1]!=S%>%pull(1))){
      stop("Column names of D_cand (except from ID) must be equal to first row of S")
    }
  }
  if(force_n_boot==FALSE && factorial(length(names(sigma_cand)[-1])*2)/(factorial(length(names(sigma_cand)[-1]))*factorial(length(names(sigma_cand)[-1])*2-length(names(sigma_cand)[-1])))
  <200){
    message(paste0("Warning: with ",length(names(sigma_cand)[-1])," responses, the p-value of the pairwise comparisons used in the ranking is bounded to be greater than ",round(1/(factorial(length(names(sigma_cand)[-1])*2)/(factorial(length(names(sigma_cand)[-1]))*factorial(length(names(sigma_cand)[-1])*2-length(names(sigma_cand)[-1])))),digits=5) ,
                        ". We suggest to set force_n_boot==TRUE."))
  }


  #calculate mean variable importance
  # VARIMP
  if(is.null(S)){
    S_varimp=NULL
  }else{
    S_varimp<-cbind(S, meanImportance=rowMeans(S[,-1]))
    if(sum(S_varimp$meanImportance)==0){
      message("Warning: mean importance of all factors is 0. varimp_distance is automatically set to FALSE")
	  varimp_distance=FALSE
    }
  }

  #hierarchical clustering
  hier_clust<-NULL
  if(varimp_distance==FALSE){
    #hierarchical clustering
    distances<-NULL
    distances<-dist(as.matrix(D_cand[,-1]), method = 'euclidean')
    clust_opt<-fviz_nbclust(as.matrix(D_cand[,-1]), hcut, k.max = nrow(D_cand)/2, method = "silhouette")
    if(is.null(n_clust)){
      best_nclust=as.numeric(as.character(clust_opt$data[which(clust_opt$data$y==max(clust_opt$data$y)),1]))
      nclust_criterion="silhouette"
    }else{
      best_nclust=n_clust
      nclust_criterion="set by analyst"
    }
    hier_clust<-hclust(distances, method = 'centroid')
  }
  if(varimp_distance==TRUE){
    distances<-NULL
    distances<-dist(as.matrix(D_cand[,-1]), method = 'euclidean')
    normalized_varimp<-NULL
    normalized_varimp<-S_varimp$meanImportance/sum(S_varimp$meanImportance)

    weighted.euc.dist <- function(x1, x2) sqrt(sum(normalized_varimp*(x1 - x2) ^ 2))
    distance_custom<-dist_make((as.matrix(D_cand[,-1])), weighted.euc.dist)

    #potentially print best number of clusters and plot
    clust_opt<-fviz_nbclust(as.matrix(D_cand[,-1]), hcut, k.max = nrow(D_cand)/2, method = "silhouette")
    #if nclust not supplied, then calculate best cluster and use that. Else use supplied nclust
    if(is.null(n_clust)){
      best_nclust=as.numeric(as.character(clust_opt$data[which(clust_opt$data$y==max(clust_opt$data$y)),1]))
      nclust_criterion="silhouette"
    }else{
      best_nclust=n_clust
      nclust_criterion="set by analyst"
    }
    hier_clust<-hclust(distance_custom, method = 'centroid')
  }


  #compute NPC rank
  sigma_cand_long<-sigma_cand%>%pivot_longer(!names(sigma_cand)[1], names_to = "Ys", values_to = "sigma_pred")
  #source("test_stat_NPC_r.R")
  risultati<-NULL
  t0<-Sys.time()
  print("Start NPC ranking")
  print(t0)
  risultati<-pair_comp(as.data.frame(sigma_cand_long%>%select(-Ys)),B=n_boot,alpha.ranking=alpha_rank,seed=seed_rank,st="dm",paral,force_n_boot)
  t1<-Sys.time()
  print(t1-t0)
  print("End NPC ranking")
  risultati_glob<-tibble(names(risultati[[2]]), "rank"=risultati$globale)
  names(risultati_glob)[1]=names(sigma_cand)[1]


  LAMBDA<-NULL
  LAMBDA<-as_tibble(cbind(full_join(D_cand, risultati_glob, by=names(sigma_cand)[1]),cluster=cutree(hier_clust, k=best_nclust))%>%arrange(desc(rank)))


  count_clust_rank<-LAMBDA%>%group_by(rank, cluster)%>%count
  LAMBDA_count<-left_join(LAMBDA, count_clust_rank, by=c("rank","cluster"))
  pred_se_tmp<-NULL
  pred_se_tmp<-sigma_cand%>%mutate(mean_sigma_pred=rowMeans(sigma_cand[,-1]))
  LAMBDA_unique<-left_join(LAMBDA_count,pred_se_tmp%>%select(names(sigma_cand)[1],mean_sigma_pred), by=names(sigma_cand)[1])%>%
    group_by(rank, cluster, n)%>%mutate(MAX_mean_sigma_pred=max(mean_sigma_pred))%>%
    filter(MAX_mean_sigma_pred==mean_sigma_pred)%>%arrange(desc(rank), desc(mean_sigma_pred))%>%
    relocate(names(sigma_cand)[1], .before=NULL)



  if(n_add>nrow(LAMBDA_unique)){
    message("Warning: n_add is larger than the number of unique candidates in LAMBDA")
  }

  n_batch<-min(n_add, nrow(LAMBDA_unique))

  D_add_j<-NULL
  #do not perform replicates if more configurations belong to same cluster and have same rank
  if(strategy=="exploration"){
    D_add_j<-D_cand[(LAMBDA_unique%>%pull(names(sigma_cand)[1]))[1:n_batch],]
  }
  #do perform replicates if more configurations belong to same cluster and have same rank
  if(strategy=="exploitation"){
    D_add_j<-D_cand[rep(LAMBDA_unique%>%pull(names(sigma_cand)[1]), LAMBDA_unique%>%pull(n))[1:n_batch],]
  }
  return(list(strategy=strategy,seed_rank=seed_rank,n_permut_comparisons=risultati$N_permut_unique,varimp_distance=varimp_distance,clust_choice=as_tibble(clust_opt$data),nclust_criterion=nclust_criterion,best_nclust=best_nclust, LAMBDA = LAMBDA , D_add_j = D_add_j))
}
