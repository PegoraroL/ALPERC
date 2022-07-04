#funzione statistica test
test.stat <- function(y, stat, alternative, B, force_n_boot){
  #source("t2p.r")

  #riscrivo label come 1 e 2
  label.orig <- y[,1]
  size <- table(label.orig)
  label = vector(length=length(label.orig))
  label[label.orig==unique(label.orig)[1]] <- 1 #rep(1:2, size)
  label[label.orig==unique(label.orig)[2]] <- 2

  #escludo colonna gruppi da miei dati
  obs <- y[,2:ncol(y)]
  if((ncol(y)-1)==1) dim(obs) <- c(nrow(y),(ncol(y)-1)) #se ? un vettore lo voglio comunque a 2 dimensioni

  #compute total nr of possible combinations for difference in means
  N_combinat=factorial(nrow(y))/(factorial(min(table(label)))*factorial(nrow(y)-min(table(label))))

  if(B > N_combinat && force_n_boot == FALSE){
    N_permut_comparisons=N_combinat
    combi<-combinat::combn(nrow(y),min(table(label)))
    all_combin_tmp1<-rbind(combi[,1:(dim(combi)[2]/2)],combi[,dim(combi)[2]:(1+(dim(combi)[2]/2))])
    all_combin_tmp2<-rbind(combi[,dim(combi)[2]:(1+(dim(combi)[2]/2))],combi[,1:(dim(combi)[2]/2)])

    all_combin<-cbind(all_combin_tmp1,all_combin_tmp2)

    T <- array(dim=c((N_combinat),ncol(obs)))
    dim(T) <- c((N_combinat),ncol(obs))
    #differenza in media
    T[1,] <- apply(obs, 2, function(x) mean(x[label==2], na.rm=TRUE) - mean(x[label==1], na.rm=TRUE))

    for(bb in 2:(N_combinat)){ #up until N_combinat
      u <- all_combin[,bb] #take all exact combinations instead of random permuting
      obs.star <- obs[u,1:dim(obs)[2]]#cbind(label, obs[u,1:dim(obs)[2]])
      if((ncol(y)-1)==1) dim(obs.star) = c(nrow(obs),ncol(obs))
      T[bb,] <-apply(obs.star, 2, function(x) mean(x[label==2],na.rm=TRUE) - mean(x[label==1],na.rm=TRUE))
    }
    #create vector of NA to avoid problem with dimensions
    NA_vec<-rep(NA,(B+1)-N_combinat)
    dim(NA_vec)<-c(length(NA_vec),1)
    P <- rbind(t2p(T),NA_vec)
    P.min <- rbind(t2p(-T),NA_vec)
    P.comb<-rbind(t2p(comb(t2p(T),"F")),NA_vec)
    P.comb.min<-rbind(t2p(comb(t2p(-T),"F")),NA_vec)
  }

  if(B <= N_combinat || force_n_boot == TRUE){
    N_permut_comparisons=B
    T <- array(dim=c((B+1),ncol(obs)))
    dim(T) <- c((B+1),ncol(obs))
    #differenza in media
    T[1,] <- apply(obs, 2, function(x) mean(x[label==2], na.rm=TRUE) - mean(x[label==1], na.rm=TRUE))

    for(bb in 2:(B+1)){
      u <- sample(1:dim(obs)[1], dim(obs)[1])
      obs.star <- obs[u,1:dim(obs)[2]]#cbind(label, obs[u,1:dim(obs)[2]])
      if((ncol(y)-1)==1) dim(obs.star) = c(nrow(obs),ncol(obs))
      T[bb,] <-apply(obs.star, 2, function(x) mean(x[label==2],na.rm=TRUE) - mean(x[label==1],na.rm=TRUE))
    }
    P <- t2p(T)
    P.min <- t2p(-T)
    P.comb<-t2p(comb(P,"F"))
    P.comb.min<-t2p(comb(P.min,"F"))
  }

  return(list(pv = P, pv.min = P.min, pv.comb = P.comb, pv.comb.min = P.comb.min, N_permut_comparisons=N_permut_comparisons))
}
