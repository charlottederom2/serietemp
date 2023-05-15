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
install.packages('ellipsis')
install.packages('ellipse')
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
library(base)
require(base)
require(tseries)
library(base)
require(base)

rm(list=ls())


###### Partie 1 : Les donnees ######

#Extraction des donnees : 
datafile <- "serietemp/patates.csv"
data <- read.csv('patates_bisbis.csv',sep=";")

#mise en forme des donnees
dates_char <- as.character(data$dates)
dates_char[1];dates_char[length(dates_char)] #affiche la premiere et la derniere date
dates <- as.yearmon(seq(from=2023+2/12,to=2005,by=-1/12)) #index des dates pour l'indice (colonne nommee 'spread' dans le csv)
serie <- zoo(data$spread,order.by=dates)
T <- length(serie)
data <- serie[1:(T-4)] #on enlève les quatre dernieres valeurs pour effectuer la prediction finale
serie_diff <- diff(serie,1) #difference premiere

#### Q1 : Representation de la serie 
plot(serie,xlab="Dates",ylab="Indice de production industrielle",main="Indice")#on represente la serie

plot(cbind(serie,serie_diff)) #on represente la serie et la serie differenciee pour avoir une idee

fit1<-decompose(serie)
plot(fit1)#on represente la serie decomposee

#2 

summary(lm(serie~dates))
#regression lineaire simple de la serie chronologique s en fonction de sa position temporelle seq(1,n)
#calcul des coefficients de regression et des statistiques associees a cette relation lineaire
#on voit que les coefficients associes a la tendance et a la constante sont significativement non nuls



### Tests de stationnarite :

## dickey fuller
adf <- adfTest(serie, lag=0, type="ct") #on met le type 'ct' pour avoir un test ADF dans le cas avec constante et tendance

#On verifie l'autocorrelation des residus jusqu'a l'ordre k 
Qtests <- function(series, k, fitdf=0) {
  pvals <- apply(matrix(1:k), 1, FUN=function(l) {
    pval <- if (l<=fitdf) NA else Box.test(series, lag=l, type="Ljung-Box", fitdf=fitdf)$p.value
    return(c("lag"=l,"pval"=pval))
  })
  return(t(pvals))
}

adfTest_valid <- function(series,kmax,type){ #tests ADF jusqu’`a des r´esidus non autocorr´el´es
  k <- 0
  noautocorr <- 0
  while (noautocorr==0){
    cat(paste0("ADF with ",k, " lags: residuals OK? "))
    adf <- adfTest(series,lags=k,type=type)
    pvals <- Qtests(adf@test$lm$residuals,24,fitdf=length(adf@test$lm$coefficients))[,2]
    if (sum(pvals<0.05,na.rm=T) == 0) {
      noautocorr <- 1; cat("OK \n")}
    else cat("nope \n")
    k <- k + 1
  }
  return(adf)
}
adf <- adfTest_valid(serie,24,"ct")
adf
#l'hypothese de racine unitaire n'est pas rejetee


## Philippe perron
pp.test(serie) 
#pp test sur la serie (cas general : avec constante et tendance)


## test KPSS
kpss.test(serie,null="Trend") #prend en compte à la fois une tendance et une constante
#on rejette H0 (H0: la serie est stationnaire) a  5% et 10%
#la serie n'est donc pas stationnaire


#### Q2 #### Stationnariser la serie
serie_diff <- diff(serie) #on differencie notre serie
plot(serie_diff, xaxt="s", type = "l") #

summary(lm(serie_diff~dates[-1])) #on enleve la première date car la série est differenciee
#on trouve une serie avec tendance et constante nulles

#test de stationnarite (H0: n'est pas stationnaire ; H1: statio)
adf <- adfTest_valid(serie_diff,24,"nc")
adf #On rejette H0 sur tous les niveaux de confiance usuels 

pp.test(serie_diff) #Idem

kpss.test(serie_diff,null="Trend") #On ne rejette pas H0 (= Est stationnaire)

#### Q3 #### représentation de la serie brute et de la serie differenciee
plot(cbind(serie,serie_diff),main="Representation des deux series") 


###### Partie 2 : Modèles ARMA #######

#### Q4 ####

#PACF et ACF pour avoir les ordres maximaux de p et q
acf(as.numeric(serie_diff) , main="ACF de s_diff")
pacf(as.numeric(serie_diff), main = "PACF de s_diff")
#les ordres maximaux sont p*=4 et q* = 1

# On vérifie avec le test de LjungBox

# Validation du modele

Qtests <- function(series, k, fitdf=0) {
  pvals <- apply(matrix(1:k), 1, FUN=function(l) {
    pval <- if (l<=fitdf) NA else Box.test(series, lag=l, type="Ljung-Box", fitdf=fitdf)$p.value
    return(c("lag"=l,"pval"=pval))
  })
  return(t(pvals))
  
}

#on utilise une fonction pour automatiser le calcul de ratios
#fonction de test des significations individuelles des coefficients

signif <- function(estim){ 
  coef <- estim$coef 
  se <- sqrt(diag(estim$var.coef))
  t <- coef/se 
  pval <- (1-pnorm(abs(t)))*2
  return(rbind(coef,se,pval))
}

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

y<-serie_diff

estim <- arima(y,c(4,0,1)); arimafit(estim)
#pas bien ajuste

estim <- arima(y,c(3,0,1)); arimafit(estim)
estim
arma31 <- arima(y,c(3,0,1),include.mean = F)
#bien ajuste et valide ##OK

estim <- arima(y,c(2,0,1)); arimafit(estim)
#pas bien ajusté 

estim <- arima(y,c(1,0,1)); arimafit(estim)
#bien ajuste et valide ##OK
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
# bien ajusté mais pas valide 


# on a selectionne 2 modeles arma11 et arma31 valides et bien ajustes. On regarde les tests BIC et AIC pour selectionner le meilleur modele.
models <- c("arma11","arma31"); names(models) <- models
apply(as.matrix(models),1, function(m) c("AIC"=AIC(get(m)), "BIC"=BIC(get(m))))

#AIC et BIC selectionnent chacun un modele different, il reste deux modeles candidats que l'on departage par le R

#Calcul du R^2 ajuste
adj_r2 <- function(model){
  ss_res <- sum(model$residuals^2)
  p <- model$arma[1]
  q <- model$arma[2]
  ss_tot <- sum(y[-c(1:max(p,q))]^2)
  n <- model$nobs-max(p,q)
  adj_r2 <- 1-(ss_res/(n-p-q-1))/(ss_tot/(n-1))
  return(adj_r2)
}
adj_r2(arma11)
adj_r2(arma31)

#on selectionne le arma31 avec le plus grand R^2


#### Q5 ####
arima311<-arima(serie,c(3,1,1),include.mean=F)
arima311



##### Partie3:Previsions #####

#### Question 7 #### on teste les hypotheses qui nous ont permis de trouver une region de confiance

####hypothese sur les residus:

#on represente la densite de nos residus
plot(density(arma31$residuals,lwd=0.5),xlim=c(-10,10),main="Densite des residus",xlab="Valeurs",ylab="Densite")

#on represente une gaussienne de meme moyenne et variance que nos residus
mu<-mean(arma31$residuals)
sigma<-sd(arma31$residuals)
x<-seq(-10,10)
y<-dnorm(x,mu,sigma)
lines(x,y,lwd=0.5,col="blue")

####hypothese sur le modele:

#extraction des coefficients du modele pour verifier les racines des polynomes phi et theta de l'arma31
arma31$coef

# Obtention des coefficients AR et MA
phi <- coef(arma31)[1:3]
theta <- coef(arma31)[4]

# Calcul des racines du polynome AR
ar_roots <- polyroot(c(1, -phi))

# Calcul des racines du polynome MA
ma_roots <- polyroot(c(1, theta))

# Affichage des racines et de leur module : les racines sont bien en dehors du cercle unite
print(ar_roots)
Mod(ar_roots)

print(ma_roots)

#les racines ont toutes un module supérieur à 1, elles sont toutes en dehors du cercle unité


#### Question 8 : representation de la region de confiance

#On trace la région de confiance pour la serie a 95%

library(forecast)
fore=forecast(arima311,h=4,level=95) #h correspond au nombre de valeurs a predire, ici h=4 comme le nombre retire au debut
par(mfrow=c(1,1))
plot(fore,col=1,fcol=2,xlim=c(2020,2024),shaded=TRUE,xlab="Temps",ylab="Valeur",main="Prevision pour un ARIMA(3,1,1)")

#On peut representer la region de confiance bivariee a 95% avec une ellipse

XT1=predict(arma31,n.ahead=2)$pred[1]
XT2=predict(arma31,n.ahead=2)$pred[2]
XT1
XT2

require(ellipsis)
require(car)
require(ellipse)
library(ellipse)
arma=arima0(serie_diff,order=c(3,0,1))

arma31<-arima(serie_diff,c(3,0,1),include.mean=F)

#on a déjà extrait les coefficients de l'ARMA31 dans phi et theta
sigma2<-as.numeric(arma31$sigma2)
sigma2
phi
theta

Sigma<-matrix(c(sigma2,(phi[1]+theta[1])*sigma2,(phi[1]+theta[1])*sigma2,(1+(phi[1]+theta[1])^2)*sigma2),ncol=2)
Sigma
inv_Sigma<-solve(Sigma)


plot(XT1,XT2,xlim=c(-17,15),ylim=c(-19,17),xlab="Prevision de X(T+1)",ylab="Prevision de X(T+2)",main="Region de confiance bivariee a 95%")

lines(ellipse(Sigma,centre=c(XT1,XT2)),type="l",col="red",xlab="Xt+1",ylab="Xt+2",main="Ellipse de confiance pour(Xt+1,Xt+2)")
abline(h=XT1,v=XT2)

