# Test Viterbi
require(HMM)

mT<-matrix(c( 0.7805608, 0.2194392, 0.2219266, 0.7780734), byrow=TRUE,nrow=2)
mE<-matrix(c(0.4321041, 0.06694529, 0.34011263, 0.1608380, 0.1586833, 0.34191773, 0.06642913, 0.4329698), byrow=TRUE,nrow=2)
hmmPostEstim<-initHMM(States = c("C","N"),
                      Symbols = ALPH,
                      startProbs = statio,
                      transProbs = mT,
                      emissionProbs = mE)
tabEstim<-viterbi_log(fasta = pfu, matT = hmmPostEstim$transProbs, matE = hmmPostEstim$emissionProbs, alphabet= ALPH, depart =statio)
tabEstimHMM<-viterbi(hmmPostEstim,pfu)
# On examine la distribution de la longueur des gènes déterminés pour
#  évaluer l'estimation markovienne


tailleSeq<-function(vectorViterbi,indexSymbole){
  symboleOppose<-ifelse(indexSymbole == 1,2,1)
  vectorWhich<-which(vectorViterbi == symboleOppose)
  resultant<-c(vectorWhich,NA) - c(0,vectorWhich)
  resultant<-resultant[ which( resultant != 1) ]
  return(resultant-1)
}

longueurSeqC<-tailleSeq(vectorViterbi = tabEstim$S, 1)

# La distribution des longueurs n'est pas réaliste.
# Aucun des gènes de dépassent une taille de 100 base, ce qui est improbable.
hist(longueurSeqC, breaks = 1000)

###
# On essaie de restreindre le modèle pour que l'état codant soit
#  plus difficile à quitter.
# On définit ainsi la transition C->C à 0.99

mTstringent<-matrix(c( 0.99, 0.01,0.87, 0.13), byrow=TRUE,nrow=2)
mE<-matrix(c(0.4321041, 0.06694529, 0.34011263, 0.1608380, 0.1586833, 0.34191773, 0.06642913, 0.4329698), byrow=TRUE,nrow=2)
hmmPostEstimStringent<-initHMM(States = c("C","N"),
                      Symbols = ALPH,
                      startProbs = statio,
                      transProbs = mTstringent,
                      emissionProbs = mE)
tabEstimStringent<-viterbi_log(fasta = pfu,
                               matT = hmmPostEstimStringent$transProbs,
                               matE = hmmPostEstimStringent$emissionProbs,
                               alphabet= ALPH, depart =statio)


longueurSeqCStringent<-tailleSeq(vectorViterbi = tabEstimStringent$S, 1)
# Cependant c'est un modèle bien trop stringent. Seul l'état codant est réprésenté
hist(longueurSeqCStringent, breaks = 1000)


mTstringent2<-matrix(c( 0.87, 0.13,0.87, 0.13), byrow=TRUE,nrow=2)
mE<-matrix(c(0.4321041, 0.06694529, 0.34011263, 0.1608380, 0.1586833, 0.34191773, 0.06642913, 0.4329698), byrow=TRUE,nrow=2)
hmmPostEstimStringent2<-initHMM(States = c("C","N"),
                               Symbols = ALPH,
                               startProbs = statio,
                               transProbs = mTstringent2,
                               emissionProbs = mE)
tabEstimStringent2<-viterbi_log(fasta = pfu,
                               matT = hmmPostEstimStringent2$transProbs,
                               matE = hmmPostEstimStringent2$emissionProbs,
                               alphabet= ALPH, depart =statio)
summary(tabEstimStringent2$S)
# toujours trop stringent

mTsemiAleat<-matrix(c( 0.87, 0.13,0.5, 0.5), byrow=TRUE,nrow=2)
mE<-matrix(c(0.4321041, 0.06694529, 0.34011263, 0.1608380, 0.1586833, 0.34191773, 0.06642913, 0.4329698), byrow=TRUE,nrow=2)
hmmPostEstimSemiAleat<-initHMM(States = c("C","N"),
                                Symbols = ALPH,
                                startProbs = statio,
                                transProbs = mTsemiAleat,
                                emissionProbs = mE)
tabEstimSemiAleat<-viterbi_log(fasta = pfu,
                                matT = hmmPostEstimSemiAleat$transProbs,
                                matE = hmmPostEstimSemiAleat$emissionProbs,
                                alphabet= ALPH, depart =statio)
summary(tabEstimSemiAleat$S)


longueurSeqCSemiAleat<-tailleSeq(vectorViterbi = tabEstimSemiAleat$S, 1)
hist(longueurSeqCSemiAleat, breaks = 1000)
# C'est pas trop mal mais 200 en gène le plus long pas assez.
# La sortie de régions codantes n'est pas assez

###############
mTPenalise<-matrix(c( 0.99, 0.01,0.5, 0.5), byrow=TRUE,nrow=2)
mE<-matrix(c(0.4321041, 0.06694529, 0.34011263, 0.1608380, 0.1586833, 0.34191773, 0.06642913, 0.4329698), byrow=TRUE,nrow=2)
hmmPostEstimPenalise<-initHMM(States = c("C","N"),
                               Symbols = ALPH,
                               startProbs = statio,
                               transProbs = mTPenalise,
                               emissionProbs = mE)
tabEstimPenalise<-viterbi_log(fasta = pfu,
                               matT = hmmPostEstimPenalise$transProbs,
                               matE = hmmPostEstimPenalise$emissionProbs,
                               alphabet= ALPH, depart =statio)
summary(tabEstimPenalise$S)


longueurSeqCPenalise<-tailleSeq(vectorViterbi = tabEstimPenalise$S, 1)
hist(longueurSeqCPenalise, breaks = 1000)


########### On pénalise désormais les zones intergéniques
###############
mTPenaliseNC<-matrix(c( 0.99, 0.01,0.60, 0.4), byrow=TRUE,nrow=2)
mE<-matrix(c(0.4321041, 0.06694529, 0.34011263, 0.1608380, 0.1586833, 0.34191773, 0.06642913, 0.4329698), byrow=TRUE,nrow=2)
hmmPostEstimPenaliseNC<-initHMM(States = c("C","N"),
                              Symbols = ALPH,
                              startProbs = statio,
                              transProbs = mTPenaliseNC,
                              emissionProbs = mE)
tabEstimPenaliseNC<-viterbi_log(fasta = pfu,
                              matT = hmmPostEstimPenaliseNC$transProbs,
                              matE = hmmPostEstimPenaliseNC$emissionProbs,
                              alphabet= ALPH, depart =statio)
summary(tabEstimPenaliseNC$S)


longueurSeqCPenaliseNC<-tailleSeq(vectorViterbi = tabEstimPenaliseNC$S, 1)
hist(longueurSeqCPenaliseNC, breaks = 1000)

########### On pénalise désormais les zones intergéniques (mais un pu moins)
###############
mTPenaliseNCless<-matrix(c( 0.99, 0.01,0.65, 0.35), byrow=TRUE,nrow=2)
mE<-matrix(c(0.4321041, 0.06694529, 0.34011263, 0.1608380, 0.1586833, 0.34191773, 0.06642913, 0.4329698), byrow=TRUE,nrow=2)
hmmPostEstimPenaliseNCless<-initHMM(States = c("C","N"),
                                Symbols = ALPH,
                                startProbs = statio,
                                transProbs = mTPenaliseNCless,
                                emissionProbs = mE)
tabEstimPenaliseNCless<-viterbi_log(fasta = pfu,
                                matT = hmmPostEstimPenaliseNCless$transProbs,
                                matE = hmmPostEstimPenaliseNCless$emissionProbs,
                                alphabet= ALPH, depart =statio)
summary(tabEstimPenaliseNCless$S)


longueurSeqCPenaliseNCless<-tailleSeq(vectorViterbi = tabEstimPenaliseNCless$S, 1)
hist(longueurSeqCPenaliseNCless, breaks = 1000)

############## 15/09/2015

mTPenaliseNCajust<-matrix(c( 0.99, 0.01,0.61, 0.39), byrow=TRUE,nrow=2)
mE<-matrix(c(0.4321041, 0.06694529, 0.34011263, 0.1608380, 0.1586833, 0.34191773, 0.06642913, 0.4329698), byrow=TRUE,nrow=2)
hmmPostEstimPenaliseNCajust<-initHMM(States = c("C","N"),
                                    Symbols = ALPH,
                                    startProbs = statio,
                                    transProbs = mTPenaliseNCajust,
                                    emissionProbs = mE)
tabEstimPenaliseNCajust<-viterbi_log(fasta = pfu,
                                    matT = hmmPostEstimPenaliseNCajust$transProbs,
                                    matE = hmmPostEstimPenaliseNCajust$emissionProbs,
                                    alphabet= ALPH, depart =statio)
summary(tabEstimPenaliseNCajust$S)
longueurSeqCPenaliseNCajust<-tailleSeq(vectorViterbi = tabEstimPenaliseNCajust$S, 1)
hist(longueurSeqCPenaliseNCajust, breaks = 1000)


## #Changmeent de la matrice d'émissions

mTPenaliseNCemis<-matrix(c( 0.99, 0.01,0.61, 0.39), byrow=TRUE,nrow=2)
hmmPostEstimPenaliseNCemis<-initHMM(States = c("C","N"),
                                     Symbols = ALPH,
                                     startProbs = statio,
                                     transProbs = mTPenaliseNCemis)
tabEstimPenaliseNCemis<-viterbi_log(fasta = pfu,
                                     matT = hmmPostEstimPenaliseNCemis$transProbs,
                                     matE = hmmPostEstimPenaliseNCemis$emissionProbs,
                                     alphabet= ALPH, depart =statio)
summary(tabEstimPenaliseNCemis$S)
longueurSeqCPenaliseNCemis<-tailleSeq(vectorViterbi = tabEstimPenaliseNCemis$S, 1)
hist(longueurSeqCPenaliseNCemis, breaks = 1000)



mTcodant<-matrix(c( 0.99, 0.01,0.61, 0.39), byrow=TRUE,nrow=2)
mE<-matrix(rep(c(0.3,0.2,0.2,0.3),2), byrow=TRUE,nrow=2)

hmmPostEstimCodant<-initHMM(States = c("C","N"),
                                    Symbols = ALPH,
                                    startProbs = statio,
                                    transProbs = mTcodant,emissionProbs = mE)
tabEstimCodant<-viterbi_log(fasta = pfu,
                                    matT = hmmPostEstimCodant$transProbs,
                                    matE = hmmPostEstimCodant$emissionProbs,
                                    alphabet= ALPH, depart =statio)
summary(tabEstimCodant$S)
longueurSeqCodant<-tailleSeq(vectorViterbi = tabEstimCodant$S, 1)
hist(longueurSeqCodant, breaks = 1000)