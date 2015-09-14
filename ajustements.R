# Librairie seqinr
require(seqinr)

# Chargement des fonctions définis
source("fonctions.R")

# Question A
# Pour des modèles d'ordre de 0 à m




ajustementRapportVrais<-function(mMax, fasta, alphabet, seuil=0.05){
  # En prenant le test de rapport de
  #  vraisemblance comme critère d'ajustement, 
  #  on cherche l'ordre de chaîne de Markov
  #  pour lequel on ne peut plus rejetter H0.
  #
  
  # Tableau trace 
  # R commence à 0; et m+1 nécessaire pour comparer
  Trace<-data.frame(matrix(nrow = mMax+1, 
                           ncol = 1, 
                           dimnames=list(0:mMax,"pvalue")))
#   print(Trace)
  m<-0
  repeat { 
    L0<-LogVraisCorr(matTrans = estimMarkovK(fasta, m, alphabet),
                      fasta = fasta,
                      ordre = m,
                      alphabet)
    L1<-LogVraisCorr(matTrans = estimMarkovK(fasta, m+1, alphabet),
                      fasta = fasta,
                      ordre = m+1,
                      alphabet)
    
    # Calcul pvalue pour test du rapport de vraisemblance :
    # H0: modèle de markov ordre m
    # H1: modèle de markov ordre m+1
    ##################################
    # /!\ le Trace[m+1] correspond à l'ordre m (artifice de numérotation R)
    ##################################
    Trace[m+1, 1]<-rapportLogVrais(alphabet, ordre = m, L0, L1)
    
    if(is.finite(Trace[m+1,1]) ) {
      # Test si le rapport est NaN.
      
#       cat(paste0("DBG: m=",m,"  trace[m+1]=",Trace[m+1,1],"\n"))
      if ( (m == mMax) ){
        # On arrête l'ajustement lorsque l'on a
        #  atteint l'itération maximale, ou
        #  que l'on a pas incrémenter l'ordre m
        cat(paste0("Itération max. Ordre de la chaîne : ",m,"\n"))
        break
      }
      if( Trace[m+1,1] > seuil ){
        # On arrête l'ajustement lorsque l'on a
        #  l'on rejette H1. (Pvalue > seuil)
        
        # On conclut donc à l'ordre m si m = 0.
        #  et m-1 sinon.
        if ( m == 0){
          cat(paste0("Ordre chaîne Markov ajusté ",m,"\n"))
        }  else {
          cat(paste0("Ordre chaîne Markov ajusté ",m-1,"\n"))
        }
        break
      }
    } else{
      # On arrête l'ajustement si la valeur calculée est NaN
      cat(paste0("NaN calculée\nOrdre chaîne Markov ajusté ",m-1,"\n"))
      break
    }
    
   m<-m+1
  }
    
  return(Trace)
}

# ajustementRapportVrais(mMax = 4, fasta = s2c(seqTest), ALPH)

ajustementCritere<-function(critere, mMax, fasta, alphabet, eps=0.1){
  # En prenant critère d'ajustement à minimiser (AIC/BIC par ex), 
  #  on cherche l'ordre de chaîne de Markov optimal
  #  dans la limite de mMax itérations.
  #
  
  # Tableau trace 
  # R commence à 0; et m+1 nécessaire pour comparer
  Trace<-data.frame(matrix(nrow = mMax+1, 
                           ncol = 1, 
                           dimnames=list(0:mMax,critere)))
  m<-0
  repeat { 
    L0<-LogVraisCorr(matTrans = estimMarkovK(fasta, m, alphabet),
                      fasta = fasta,
                      ordre = m,
                      alphabet)
    # Estimation (grossière) du nombre de paramètres (non optimisé)
    Nparametre<-length(alphabet)*length(alphabet)^m
    
    # Valeur du critère pour l'ordre m
    Trace[m+1, 1]<-switch(critere,
                          "AIC" = AIC(L0, K = Nparametre),
                          "BIC" = BIC(L0, K = Nparametre,n=length(fasta)))
    
    # Sauvegarde du critère précédent
    #  initialisé à Inf pour m = 0.
    critere.old<-ifelse(m == 0,Inf, Trace[m,1] )
    
    cat(paste0("m: ",m,"  LL: ",L0,"  Trace:", Trace[m+1,1]," OLD: ", critere.old,"\n"))
   if(is.finite(Trace[m+1,1]) ) {
     if ( m == mMax ) {
       # On arrête l'ajustement lorsque 
       #  si l'on a atteint l'itération maximale. 
       cat(paste0("Itération max. Ordre de la chaîne avec ", critere," : ",m,"\n"))
       break
     }
      if ( Trace[m+1,1] > critere.old ){
        cat(paste0("Ordre chaîne Markov ajusté avec ",critere," : ",m-1,"\n"))
        break
      }
   } else{
     # On arrête l'ajustement si la valeur calculée est NaN
     cat(paste0("NaN calculée\nOrdre chaîne Markov ajustée ",critere," : ",m-1,"\n"))
     break
   }
    
   m<-m+1
  }
  
  return(Trace)
}

# Alphabet nucléique standard
ALPH<-c("a", "c", "g", "t")

# Données
pfu50<-read.fasta("headN50_Pfu_DSM3638.fasta")[[1]]
pfu<-read.fasta("complete_genome_Pfu_DSM3638.fasta")[[1]]


# Ajustements
pfuRapport<-ajustementRapportVrais(mMax = 7, fasta = pfu,alphabet = ALPH)
pfuAIC<-ajustementCritere("AIC",mMax = 7, fasta = pfu,ALPH)
pfuBIC<-ajustementCritere("BIC",mMax = 7, fasta = pfu,ALPH)


#tabAjustement<-data.frame(cbind(rownames(pfuRapport), pfuRapport,pfuAIC,pfuBIC))
# voir tab.tex
