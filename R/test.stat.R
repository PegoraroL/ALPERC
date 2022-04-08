#funzione statistica test
test.stat <- function(y, stat, alternative, B){
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

  return(list(pv = P, pv.min = P.min, pv.comb = P.comb, pv.comb.min = P.comb.min))
}
