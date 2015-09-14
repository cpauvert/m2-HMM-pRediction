require(seqinr)

#' Fonction de calcul de l'estimation de la matrice de transition
#' d'une chaîne de Markov d'ordre $k$.
estimMarkovK<-function(fasta, ordre, alphabet ) {
  Nij<-matrix(count(fasta, ordre+1), byrow = T, ncol=length(alphabet))
  
  if(ordre != 0){
    Ni<-count(fasta[1: (length(fasta)-ordre)], ordre)
  }
  else {
    Ni<-length(fasta)
  }
  matT<-Nij/as.vector(Ni)
  return(matT)
}

LogVrais<-function(matTrans,fasta, ordre, alphabet) {
  # Matrice des transitions observées dans la séquence
  NAlphaBeta<-matrix(count(fasta, ordre+1), byrow = T, ncol=length(alphabet))
  return(sum(NAlphaBeta*(log(matTrans))))
}
LogVraisCorr<-function(matTrans,fasta, ordre, alphabet) {
  # log vraisemblance corrigée pour éviter log(0) = -Inf
  # si on a pas toutes les transitions 
  
  # Matrice des transitions observées dans la séquence
  NAlphaBeta<-matrix(count(fasta, ordre+1), byrow = T, ncol=length(alphabet))
  return(sum(NAlphaBeta*(log1p(matTrans))))
}

rapportLogVrais<-function(alphabet, ordre, L0, L1) {
  # Attention ici on travaille avec le nombre max de paramètres
  #   estimées, pas de L(L-1) (pour le moment)
  # 
  statDeTest<-(-2*(L0-L1))
  # Nb paramètres Modèle 1 - Nb paramètres Modèle 0
  L<-(length(alphabet))
  DDL<-L*(L^(ordre+1)-L^(ordre))
  
  # P( X > statDeTest ) ~ pvalue
  pvalue<-pchisq(statDeTest, df = DDL, lower.tail = F)
  return(pvalue)
}

# Fonctions de convenances
AIC<-function(L0, K) -2*L0+2*K
BIC<-function(L0, K, n) -2*L0+K*log(n)



