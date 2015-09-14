require(seqinr)

####### Partie Markov #######
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

####### Partie Markov Caché #######
# Implémentation algorithme Viterbi


viterbi_log<-function(fasta, matT,matE,alphabet, depart) {

	# matT : matrice de transitions des états cachés
	# matE : matrice des probabilités d'émissions des symboles
       	#  	  en fonction des états cachés
  
  mV<-matrix(ncol = 4,nrow = length(fasta))
# TODO : dimnames lists matrix  
  #' col 1 : état 1
  #' col 2 : état 2
  #' col 3 : état k qui maximise (nu_t-1*a) etat 1
  
  etatCacheInitialisation<-c(1,2)
  mV[1,]<-c(log(depart)+log(matE[,fasta[1]]),etatCacheInitialisation)
  #' depart : loi initiale pour les deux états
  #' b : matrice des probas des lettres (M0) selon le modèle
  for( t in 2:(length(fasta))  ){
    #' b_y2^j : b[j,which(alphabetabet == y2)]
    #' 
    #' nu_1 : logV[1]
    #' nu_2 : logV[2]
    #' 
    #' a_1j : A[1,j]
    #' 
    #     logV<-rbind(logV,log1p(b[,which(alphabetabet==Y[1])]))[,2]
    
    #' Stockage de l'état précédent qui maximise logV
    #' A[, mV[i-1,3] ] : correspond à la colonne de la matrice
   
	# Exploration des possibilités.
	#  l'ajout de la matrice de transitions étend le champ de l'observation
	#  des valeurs de vraisemblance des deux états (colonnes 1 et 2)
	#  aux transitions possibles dans la séquence :
	#	   C->N 
	#	   C->C 
	#	   N->N 
	#	   N->C 
    MMax<-mV[t-1,1:2] +log(matT)
    

  	# On cherche l'état caché suivant qui maximise la vraisemblance
  	#  en fonction que l'on provienne de l'un ou l'autre des états cachés;
    EtatMax<-apply(MMax,2,which.max)

    	# Les valeurs des vraisemblances sont identifiées 
    ValMax<-apply(MMax,2,max)
    
	# Les états et les vraisemblances correspondantes
	#  sont ajoutées aux précédentes.
    mV[t,1:2]<-log(matE[,fasta[t]]) +  ValMax
    mV[t,3:4]<-EtatMax
  }
  
  # Setoile correspond au vecteur de taille des données
  #  illustrant les états cachés déduits.
  Setoile<-vector("numeric", length(fasta))

  # Initialisation du dernier état caché qui maximise la vraisemblance
  Setoile[length(fasta)]<-which.max(mV[length(fasta),1:2])
  
  # Remontée de la séquence en sens inverse.
  for(t in seq(length(fasta),2)){
    
    # tau(t) : mv[, 3:4]
    ST<-2+Setoile[t] # S étoile en t, 2 pour choisir colonne 3 ou 4
    Setoile[t-1]<-mV[ t , ST ]
    #which.max(mV[i,1:2])
  }
  tmp<-data.frame(cbind(mV,Setoile))
  names(tmp)<-c("LL1", "LL2","S1","S2","S")
  return(tmp)
  
}
