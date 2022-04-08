#funzione per combinazione
comb<-function(pv,fcomb,tau=FALSE,a=FALSE){
  if(fcomb=="F"){
    #cat('Combination Function: Fisher \n')
    fi<-apply(pv,1,function(x) -2*log(prod(x,na.rm=TRUE)))##Fisher
  }
  if(fcomb=="L"){
    cat('Combination Function: Liptak \n')
    fi<-apply(pv,1,function(x) sum(qnorm(1-x),na.rm=TRUE))##Liptak
  }
  if(fcomb=="T"){
    #cat('Combination Function: Tippet \n')
    fi<-apply(pv,1,function(x) max((1-x),na.rm=TRUE)) ##Tippet
  }

  if(fcomb=="Max_T"){
    fi<-apply(pv,1,function(x) max((x),na.rm=TRUE)) ##Max
  }

  if(fcomb=="Average"){
    fi<-apply(pv,1,function(x) sum((x),na.rm=TRUE))
  }
  if(fcomb=="concave"){
    fi<-apply(pv,1,function(x) sum(qlnorm(x),na.rm=TRUE))
  }
  if(fcomb=="Truncated"){
    fi<-apply(pv,1,function(x) prod(x^(x<tau)))
  }
  if(fcomb=="max_truncated"){
    ind<-1:dim(pv)[2]
    ind_max<-ind[pv[1,]==max(pv[1,])]
    if(length(ind_max)>1){
      fi<-pv[,ind_max[1]]
    }
    if(length(ind_max)==1) fi<-pv[,ind_max]
  }
  if(fcomb=="weight"){
    fi<-apply(pv,1,function(x) sum(x^a))
  }
  if(fcomb=="sum_prod"){
    fi<-apply(pv,1,function(x) {sum(x)-prod(x)})
  }
  return(fi)

}
