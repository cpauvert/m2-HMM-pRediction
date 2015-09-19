#' Test de l'estimation des chaînes de markov via 
#' le package \verb+HMM+.
require(HMM)
source('fonctions.R')

# Alphabet nucléique standard
ALPH<-c("a", "c", "g", "t")
# Données
pfu50<-read.fasta("headN50_Pfu_DSM3638.fasta")[[1]]
pfu<-read.fasta("complete_genome_Pfu_DSM3638.fasta")[[1]]

#' Estimation de la  loi stationnaire du modèle.
#' Probabilité d'être dans l'état C ``codant'' estimé par
#' la longueur moyenne des gènes du génome rapporté à la taille totale du génome.
statio<-c(0.87, 0.13)
names(statio)<-c('C','N')

#' Initialisation d'une chaîne de Markov
#' Les matrices suivantes proviennent de l'estimation réalisé sur genotoul.
#'  avec 
#'   > hmm1<-initHMM(States = c("C","N"), Symbols = ALPH,startProbs = statio)
#'   > estimHMM1<-baumWelch(hmm = hmm1,observation = pfu)
#'   
mT<-matrix(c( 0.7805608, 0.2194392,
              0.2219266, 0.7780734),
           byrow=TRUE,nrow=2)
mE<-matrix(c(0.4321041, 0.06694529, 0.34011263, 0.1608380,
             0.1586833, 0.34191773, 0.06642913, 0.4329698),
           byrow=TRUE,nrow=2)
hmmPostEstim<-initHMM(States = c("C","N"),
                      Symbols = ALPH,
                      startProbs = statio,
                      transProbs = mT,
                      emissionProbs = mE)
tabEstim<-viterbi_log(fasta = pfu,
                      matT = hmmPostEstim$transProbs,
                      matE = hmmPostEstim$emissionProbs,
                      alphabet= ALPH,
                      depart = statio)

#' Longueur des séquences codantes.
longueurSeqC<-tailleSeq(vectorViterbi = tabEstim$S, 1)

# La distribution des longueurs n'est pas réaliste.
# Aucun des gènes de dépassent une taille de 100 base, ce qui est improbable.
hist(longueurSeqC, breaks = 1000)

LogVraisHMM(tabEstim) # -2.777997e+12


### Deuxième estimation ###
mTfinal<-matrix(c(0.8795310, 0.1204690,
                  0.2681361, 0.7318639), byrow=TRUE,nrow=2)
mEfinal<-matrix(c(0.2173744, 0.28504395, 0.1179185, 0.3796631,
                  0.4715304, 0.02250706, 0.3957343, 0.1102283),
                byrow=TRUE,nrow=2)

hmmPostEstimfinal<-initHMM(States = c("C","N"),
                           Symbols = ALPH,
                           startProbs = statio,
                           transProbs = mTfinal,
                           emissionProbs = mEfinal)
tabEstimfinal<-viterbi_log(fasta = pfu,
                           matT = hmmPostEstimfinal$transProbs,
                           matE = hmmPostEstimfinal$emissionProbs,
                           alphabet= ALPH, depart =statio)

longueurSeqfinal<-tailleSeq(vectorViterbi = tabEstimfinal$S, 1)
hist(longueurSeqfinal, breaks = 1000)

LogVraisHMM(tabEstimfinal) # -2.707342e+12 

# Pour générer les histogrammes directement avec le fichier LaTeX.
saveRDS(tabEstim,'tEstim.rds')
saveRDS(tabEstimfinal,'tFinal.rds')
