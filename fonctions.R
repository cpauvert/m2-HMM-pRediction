require(seqinr)

#' Fonction de calcul de l'estimation de la matrice de transition
#' d'une chaîne de Markov d'ordre $k$.

estimMarkovK<-function(fasta, ordre, alphabet ) {
  #ordre=1
  tabCount<-count(fasta[-length(fasta)], ordre)
  # attention à ne pas prendre en compte la dernière lettre
  tabCountPlus<-count(fasta, ordre+1)
  Nij<-matrix(tabCountPlus, byrow = T, ncol=length(alphabet), dimnames = list( alphabet, names(tabCountPlus)))
  matT<-Nij/as.vector(tabCount)
  return(matT)
}