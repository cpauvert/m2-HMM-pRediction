#' Test de l'estimation des chaînes de markov via 
#' le package \verb+HMM+.
#+ data,echo=TRUE
require(HMM)
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
hmm1<-initHMM(States = c("C","N"), Symbols = ALPH,startProbs = statio)
estimHMM1<-baumWelch(hmm = hmm1,observation = pfu)
