#install.packages('HiddenMarkov')
require(HiddenMarkov)

#' Il faut définir un objet contenant notre modèle de Markov.
#' Est ce que l'on est en temps discret  ?
#' Oui puisque notre v.a $Y$ est défini comme une des lettres de 
#' l'alphabet $\mathcal{A} = \{ a, c, g, t\}$ de la séquence $\mathbb{X}$.
#' 
#' On définit pour l'état caché de notre chaîne de Markov
#'  la variable $S$ où $S=\{C,N \}$ pour modéliser l'état codant ou
#' non.


dthmm(pfu50)
Pi <- matrix(c(0.8, 0.1, 0.1,
               0.1, 0.6, 0.3,
               0.2, 0.3, 0.5),
             byrow=TRUE, nrow=2)

#' On ne sait pas quelle matrice de transition utiliser pour
#' initialiser l'estimation des paramètres.
#' On va donc définir une fonction pour déterminer des matrices 
#' de transition aléatoirement, *avec* la contrainte que la somme 
#' des lignes de la matrice vaut $1$.

aleatMatrixA<-function(){
  x<-sample(seq(0,1,length.out = 1000), 1)
  y<-sample(seq(0,1,length.out = 1000), 1)
  return(matrix(c(x,1-x,
                  1-y,y),
                byrow=TRUE,nrow=2))  

}

matriceA<-c(0.99,0.01,0.1,0.9)
A<-matrix(matriceA,ncol=2,byrow=T)


matriceb<-c(0.4,0.1,0.1,0.4,0.05,0.4,0.5,0.05)
b<-matrix(matriceb,ncol=4,byrow=T)

x<-dthmm(NULL, A, delta = c(1,0), discrete = T, distn = "norm",pm=list(A=b[,1],C=b[,2],G=b[,3],T=b[,4]))
x<-simulate(x, nsim = 1000)


###############
require(HMM)

statio<-c(0.87, 0.13)
names(statio)<-c('C','N')

hmm1<-initHMM(States = c("C","N"), Symbols = ALPH,startProbs = statio)
estimHMM1<-baumWelch(hmm = hmm1,observation = pfu)
