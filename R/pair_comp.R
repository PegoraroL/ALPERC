#funzione principale. Ritorna lista con ranking parziali e globale. Regolare l'alpha.ranking per rendere più o meno sensibile il ranking finale.
pair_comp<-function(dati,B=2000,alpha.ranking=0.05,seed=125,st="dm", paral=FALSE){
  set.seed(seed)
  cell<-dati[,1] #colonna con i gruppi

  C <- dim(table(cell))

  #The comparisons
  pair.comp = NULL
  for(rr in 1:(C-1)){
    for(cc in (rr+1):C){
      pair.tmp <- c(rr,cc)
      pair.comp <- rbind(pair.comp,pair.tmp) # pairs of products
    }#for cc
  }#for rr
  nc = C*(C-1)/2 #no. comparisons
  P.dependent4dependent <- P.dependent4dependent.min <- array(dim=c((B+1),(ncol(dati)-1), nc))
  P.glob <- P.glob.min <- array(dim=c((B+1),1, nc))

  #n. of dependent variables per dataset
  ndependent <- ncol(dati) - 1

  if(paral==FALSE){
    # Initialize progress bar
    progbar <- txtProgressBar(min = 0, max = nc, style = 3, char = "=", width=getOption("width"))
    for(nn in 1:nc){
      data.pair <- dati[(cell==unique(cell)[pair.comp[nn,1]])|(cell==unique(cell)[pair.comp[nn,2]]),]
      res <- test.stat(y = data.pair,stat = st, alternative=1,B)
      P.dependent4dependent[,,nn] <- res$pv
      P.dependent4dependent.min[,,nn] <- res$pv.min
      P.glob[,1,nn] <- res$pv.comb
      P.glob.min[,1,nn] <- res$pv.comb.min
      #set progress bar
      setTxtProgressBar(progbar, nn)
    }
    close(progbar)
  }

  if(paral==TRUE){
    # Initialize progress bar
    progbar <- txtProgressBar(min = 0, max = nc, style = 3, char = "=", width=getOption("width"))
    print_progress <- function(n) setTxtProgressBar(progbar, n)
    opts_snow <- list(progress = print_progress)
    #here use foreach to speed up computation through parallelization (cluster must have been defined)
    result<-foreach(nn = 1:nc,.packages = 'ALPERC',.options.snow = opts_snow) %dopar%{
      # source("R/test.stat.R")
      # source("R/t2p.R")
      # source("R/comb.R")
      set.seed(seed)

      data.pair <- dati[(cell==unique(cell)[pair.comp[nn,1]])|(cell==unique(cell)[pair.comp[nn,2]]),]
      t_stat<-test.stat(y = data.pair,stat = st, alternative=1,B)
      return(t_stat)
    }
    close(progbar)

    for(nn in 1:nc){
      P.dependent4dependent[,,nn]<-result[[nn]]$pv
      P.dependent4dependent.min[,,nn] <- result[[nn]]$pv.min
      P.glob[,1,nn] <- result[[nn]]$pv.comb
      P.glob.min[,1,nn] <- result[[nn]]$pv.comb.min
    }
  }

  #Partial --> creo matrici dei p-value
  dependent4dependent <- array(dim=c(C,C,ndependent))
  a<-0
  for(i1 in 1:(C-1)){
    for(i2 in (i1+1):C){
      a <- a+1
      dependent4dependent[i1,i2,] <- p.adjust(P.dependent4dependent[1,,a],"BH")  #applico correzione molto lieve
      dependent4dependent[i2,i1,] <- p.adjust(P.dependent4dependent.min[1,,a],"BH")  #applico correzione molto lieve
    }
  }
  dimnames(dependent4dependent) <- list(paste(unique(cell)),paste(unique(cell)),names(dati)[2:(ncol(dati))])

  #Global --> creo matrici dei p-value
  global <- array(dim=c(C,C))
  a<-0
  for(i1 in 1:(C-1)){
    for(i2 in (i1+1):C){
      a <- a+1
      global[i1,i2] <- P.glob[1,,a]
      global[i2,i1] <- P.glob.min[1,,a]
    }
  }
  dimnames(global) <- list(paste(unique(cell)),paste(unique(cell)))


  #Ranking parziali
  ranking.dependent4dependent <- array(dim=c(ndependent,C))
  for(kk2 in 1:ndependent){
    diag(dependent4dependent[,,kk2])<-1
    ranking.dependent4dependent[kk2,] <- rank_populations(dependent4dependent[,,kk2],alpha.ranking)
    colnames(ranking.dependent4dependent) <- paste(unique(cell))
    rownames(ranking.dependent4dependent) <- names(dati)[2:(ncol(dati))]

  }

  #Ranking globali
  diag(global)<-1
  ranking.global <- rank_populations(global,alpha.ranking)
  names(ranking.global) <- paste(unique(cell))

  return(list(parziali = ranking.dependent4dependent , globale = ranking.global))

}

