install.packages('zoo')
install.packages('tseries')
install.packages('fUnitRoots')
install.packages('questionr')
install.packages('corrplot')
install.packages('readr')
install.packages('tidyverse')
install.packages('dplyr')
install.packages('Hmisc')
install.packages('lmtest')
install.packages('margins')
install.packages('psych')
install.packages('forecast')
require(zoo) #format de serie temporelle pratique et facile d'utilisation (mais plus volumineux)
require(tseries) #diverses fonctions sur les series temporelles
require(fUnitRoots)
require(forecast)
require(car)
require(readr)
require(questionr)
library(readr)
library(tidyverse)
library(dplyr)
library(questionr)
library(corrplot)
library(Hmisc)
library(lmtest)
library(margins)
library(psych)


###### Partie 1 : Les donnees ######

#Extraction des donnees : 
path <- "C:/Users/berti/Documents/ENSAE - 2A/Time series" 
setwd(path) 


#Mise en forme : 
datafile <- "serietemp/patates.csv"
data <- read.csv('patates.csv',sep=";")
#data<- as.data.frame(data)


### Q1 : Representation de la serie 
dates_char <- as.character(data$dates)
dates_char[1];dates_char[length(dates_char)] #affiche la premi`ere et la derni`ere date
dates <- as.yearmon(seq(from=2006+3/12,to=2023+1/12,by=1/12)) #index des dates pour spread
serie <- zoo(data$spread,order.by=dates)
serie_diff <- diff(spread,1) #diff´erence premi`ere
plot(serie,xlab="Dates",ylab="Indice de production industrielle",main="Indice")

plot(cbind(serie,serie_diff))

summary(lm(serie~dates))
#regression lineaire simple de la serie chronologique s en fonction de sa position temporelle seq(1,n)
#calcul des coefficients de regression et des statistiques associees a cette relation lineaire

monthplot(serie)

acf(serie)
pacf(serie)

fit1<-decompose(serie)
plot(fit1)



### Tests de stationnarite :

kpss.test(spread,null="Trend")

#On verifie l'autocorrelation des residus jusqu'a l'ordre k 


# Philippe perron
pp.test(spread) 
#pp test sur la serie (cas general : avec constante et tendance)

# dickey fuller
adf.test(spread)
# On ne rejette pas H0 

# test KPSS
kpss.test(spread,null="Trend") 
#on rejette H0 (H0: la serie est stationnaire) a  5% et 10%
#la serie n'est donc pas stationnaire


#### Q2 #### Stationnariser la serie
s_diff <- diff(s)
plot(s_diff, xaxt="s")


Qtests<-function(series,k,fitdf=0){
  pvals<-apply(matrix(1:k),1,FUN=function(l)
  {pval<-
    if(l<=fitdf)NA 
  else Box.test(series,lag=l,type="Ljung-Box",fitdf=fitdf)$p.value
  return(c("lag"=l,"pval"=pval))})
  return(t(pvals))}

#test de stationnarité H0: n'est pas stationnaire H1: statio
adf.test(s_diff)
pp.test(s_diff)
kpss.test(s_diff,null="Trend") 

#Enlever la moyenne
s_centre <- s_diff - mean(s_diff)

#représentation des deux séries
plot(cbind(s,s_diff),main="Representation des deux series")

acf(as.numeric(s_diff) , main="ACF de s_diff")
pacf(as.numeric(s_diff), main = "PACF de s_diff")
#les ordres maximaux sont p*=4 et q* = 1

par(mfrow = c(1,2))

acf(s_centre,24);
pacf(s_centre,24)

# statio centré donc on peut mettre ARMA
# ARMA d'ordres 4,1 à vérifier 


###### Partie 2 : Modèles ARMA #######


# On vérifie avec le test de LjungBox

# Validation du modèle

Qtests <- function(series, k, fitdf=0) {
  pvals <- apply(matrix(1:k), 1, FUN=function(l) {
    pval <- if (l<=fitdf) NA else Box.test(series, lag=l, type="Ljung-Box", fitdf=fitdf)$p.value
    return(c("lag"=l,"pval"=pval))
  })
  return(t(pvals))
}
#tests de LB pour les ordres 1 a 24

# L’absence d’autocorrelation n’est jamais rejetee a 95% jusqu’a 24 retards. Le modele est donc valide.
# dans la fonction Qtestes: parametre numero 1: series
# le test de Ljung-Box ne rejette pas l'absence d'autocorrelation des residus a l'ordre 6
# matrice de pvals de k lignes et 1 colonnes: remplie de 1, pour chaque element on applique la fonction (if... on met NA, sinon ... quand on sort de la fonction on met return)
#on trouve que pour 1 a 5: NA: le k est plus petit 
#dans ce cas on veut accepter H0 (residus non correles) donc on veut des p-valeurs elevees: on est ravis ici



#les coeff de l'ar3, on fait le ratio de 0.1748/0.1604 = 1.08 < 1.96 --> coef non significatif, on peut simplifier le modele
#On teste tous les modeles pour trouver les modeles bien ajustes et valides, on les compare tous 



summary(lm(spread ~ dates))
#on fait une regression : les p-valeurs sont tres faibles, coef significatif au seuil de 1% (on le sait aussi car il y a 3 etoiles)
#on fait le test de DF avec c et t : comme les deux coeff sont significatifs au seuil de 1% il sont non nuls donc on les inclus. 
#deltaXt = c + bt + betaX(t-1) + somme(1?k) Phi(l)DeltaXt-l + Epsilon

#install.packages("fUnitRoots")#tests de racine unitaire plus modulables
#library(fUnitRoots)

arimafit <- function(estim){
  adjust <- round(signif(estim),3)
  pvals <- Qtests(estim$residuals,24,length(estim$coef)-1)
  pvals <- matrix(apply(matrix(1:24,nrow=6),2,function(c) round(pvals[c,],3)),nrow=6)
  colnames(pvals) <- rep(c("lag", "pval"),4)
  cat("tests de nullite des coefficients :\n")
  print(adjust)
  cat("\n tests d'absence d'autocorrelation des residus : \n")
  print(pvals)
}


#### Q4 ####
#on cr?e une fonction pour automatiser le calcul de ratios
#fonction de test des significations individuelles des coefficients
#on r?cup?re les coeffs
#"on r?cup?re les std errors"
#"on calcule le ratio"
signif <- function(estim){ 
  coef <- estim$coef 
  se <- sqrt(diag(estim$var.coef))
  t <- coef/se 
  pval <- (1-pnorm(abs(t)))*2
  return(rbind(coef,se,pval))
}

##
arimafit <- function(estim){
  adjust <- round(signif(estim),3)
  pvals <- Qtests(estim$residuals,24,length(estim$coef)-1)
  pvals <- matrix(apply(matrix(1:24,nrow=6),2,function(c) round(pvals[c,],3)),nrow=6)  #on reformate le tableau pour avoir 6 lignes
  colnames(pvals) <- rep(c("lag", "pval"),4)
  cat("tests de nullite des coefficients :\n")
  print(adjust)
  cat("\n tests d'absence d'autocorrelation des residus : \n")
  print(pvals)
}

y<- s_centre

estim <- arima(y,c(4,0,1)); arimafit(estim)
#pas bien ajuste

estim <- arima(y,c(3,0,1)); arimafit(estim)
arma31 <- estim
#OK 

estim <- arima(y,c(2,0,1)); 
arimafit(estim)
#pas bien ajusté 

estim <- arima(y,c(1,0,1)); 
arimafit(estim)
#valide et bien ajuste 
arma11 <- estim

estim <- arima(y,c(0,0,1)); arimafit(estim)
# bien ajusté et pas valide

estim <- arima(y,c(1,0,0)); arimafit(estim)
# bien ajuste et pas du tout valide

estim <- arima(y,c(2,0,0)); arimafit(estim)
# bien ajuste et pas valide

estim <- arima(y,c(3,0,0)); arimafit(estim)
# bien ajuste et pas valide

estim <- arima(y,c(4,0,0)); arimafit(estim)
ar4 <-estim

# on a selectionne 3 modeles ar3, ma2 et arma21 valides et bien ajustes. On regarde les tests BIC et AIC pour selectionner le meilleur modele.
models <- c("arma11","ar4","arma31"); names(models) <- models
apply(as.matrix(models),1, function(m) c("AIC"=AIC(get(m)), "BIC"=BIC(get(m))))
#on regarde en 1 AIC et BIC cr??s pour comparer des modeles; ensuite on regarde le R
#Si AIC et BIC selectionnent 2 modeles, ils restent deux modeles candidats que l'on departage par le R

#on prend le arma11 car plus petits AIC et BIC. on regarde ensuite R pour confirmer notre intuition


##### Partie3:Previsions #####

#### Question 6

tsdiag(arma11)
qqnorm(arma11$residuals)

#jarque.bera.test(arima11$residuals)
plot(density(arma11$residuals,lwd=0.5),xlim=c(-10,10),main="Densite des residus",xlab="Valeurs",ylab="Densite")

mu<-mean(arma11$residuals)
sigma<-sd(arma11$residuals)
x<-seq(-10,10)
y<-dnorm(x,mu,sigma)
lines(x,y,lwd=0.5,col="blue")

#Extraction des coefs du modele et de la variance des residus
arma11$coef
phi_1<-as.numeric(arma11$coef[1])
phi_2<-as.numeric(arma11$coef[2])
sigma2<-as.numeric(arma11$sigma2)

phi_1
phi_2
sigma2
