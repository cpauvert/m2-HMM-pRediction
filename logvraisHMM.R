

LogVraisHMM<-function(estimationHMM){
  
  tmp<-sum(estimationHMM[which(estimationHMM$S==1),1]) +
  sum(estimationHMM[which(estimationHMM$S==2),2])
 return(tmp)
}

LogVraisHMM(tabEstim)
LogVraisHMM(tabEstimStringent)
LogVraisHMM(tabEstimStringent2)
LogVraisHMM(tabEstimSemiAleat)
LogVraisHMM(tabEstimPenalise)
LogVraisHMM(tabEstimPenaliseNC)
