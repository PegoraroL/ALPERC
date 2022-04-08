#funzione per calcolo ranking
rank_populations<-function(P,alpha,alpha_coeff=2) { # P è la matrice CxC dei p-value dei confronti a coppie direzionali
  C<-dim(P)[1]
  zeros_and_ones<-ifelse(P<alpha/alpha_coeff,1,0)
  r_u<-rank((C-apply(zeros_and_ones,1,sum)),ties.method="min")
  r_d<-(1+apply(zeros_and_ones,2,sum))
  r<-rank((apply(rbind(r_u,r_d),2,mean)),ties.method="min")
  return(r)
}
