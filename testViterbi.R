# Test Viterbi
require(HMM)

mT<-matrix(c( 0.7805608, 0.2194392, 0.2219266, 0.7780734), byrow=TRUE,nrow=2)
mE<-matrix(c(0.4321041, 0.06694529, 0.34011263, 0.1608380, 0.1586833, 0.34191773, 0.06642913, 0.4329698), byrow=TRUE,nrow=2)
hmmPostEstim<-initHMM(States = c("C","N"),
                      Symbols = ALPH,
                      startProbs = statio,
                      transProbs = mT,
                      emissionProbs = mE)
tabEstim<-viterbi_log(fasta = pfu50, matT = hmmPostEstim$transProbs, matE = hmmPostEstim$emissionProbs, alphabet= ALPH, depart =statio)
