# m2-HMM-pRediction
Projet M2 de prédictions de gènes via HMM

L'organisme modèle choisie est [_Pyrococcus furiosus_](http://www.ncbi.nlm.nih.gov/genome/?term=Pyrococcus%20furiosus). Cette archeae présente l'avantage d'avoir un génome séquencé complet et en un contig. 
De plus entre 2000 et 2012, deux prédictions protéiques différentes ont été publiés fournissant ainsi matière à comparaison pour notre projet.

## Compilation du rapport

```
Rscript -e "library(knitr);knit('CR.Rnw')" 
pdflatex CR.tex
```
