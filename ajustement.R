# Chargement des fonctions définis
source("fonctions.R")

# Question A
# Pour des modèles d'ordre de 0 à m




ajustementRapportVrais<-function(mMax, fasta, alphabet){
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
  print(Trace)
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
    Trace[m+1, 1]<-rapportLogVrais(alphabet, ordre = m, L0, L1)
    
   if(m+1 == mMax){ break }
   m<-m+1
  }
  # On arrête l'ajustement lorsque l'on rejette le modèle de complexité supérieure
  #  ou que l'on a atteint l'itération maximale.
    
  cat(paste0("Ordre chaîne Markov ajustée :",m," \n", sep=" "))
  return(Trace)
}

# ajustementRapportVrais(mMax = 4, fasta = s2c(seqTest), ALPH)

ajustementAIC<-function(mMax, fasta, alphabet, eps=0.1){
  # En prenant l'AIC comme critère d'ajustement, 
  #  on cherche l'ordre de chaîne de Markov
  #  pour lequel on atteint un limite de mMax
  #
  
  # Tableau trace 
  # R commence à 0; et m+1 nécessaire pour comparer
  Trace<-data.frame(matrix(nrow = mMax+1, 
                           ncol = 1, 
                           dimnames=list(0:mMax,"AIC")))
  m<-0
  repeat { 
    L0<-LogVraisCorr(matTrans = estimMarkovK(fasta, m, alphabet),
                      fasta = fasta,
                      ordre = m,
                      alphabet)
    Nparametre<-length(alphabet)*length(alphabet)^m
    Trace[m+1, 1]<-AIC(L0, K = Nparametre)
    
   if(is.finite(Trace[m+1,1])) {
     if ( (m == mMax) || min((Trace[m+1,1] < Trace[m,1])){     
      break
     }
   } 
   else{
     print("NaN Error !")
     break
   }
    
   m<-m+1
  }
  # On arrête l'ajustement lorsque l'on a atteint l'itération maximale.
  cat(paste0("Ordre chaîne Markov ajustée AIC :",m," \n", sep=" "))
  return(Trace)
}

ajustementAIC(mMax = 4, fasta = s2c(seqTest),ALPH)
