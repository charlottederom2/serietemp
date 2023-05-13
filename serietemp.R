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

rm(list=ls())
###### Partie 1 : Les donnees ######

#Extraction des donnees : 
path <- "C:/Users/berti/Documents/ENSAE - 2A/Time series" 
setwd(path) 

#Mise en forme : 
datafile <- "serietemp/patates.csv"
data <- read.csv('patates_bisbis.csv',sep=";")

require(zoo)
require(tseries)
library(base)
require(base)

data.source <- zoo(data[[1]], order.by=data$dates) #convertit le 1er element de data en "zoo"
T <- length(data.source)
data <- data.source[1:(length(data.source)-4)] #supprime les 4 dernieres valeurs
data_df <- data.frame(date = as.Date(index(data.source), format = "%Y-%m"), spread = coredata(data.source))

dates_char <- as.character(data_df$date)
is.unique(data_df$date)

data_df$date <- as.Date(data_df$date, format="%Y-%m") # changer le format de la date
data_df_unique <- data_df[!duplicated(data_df$date),]
data_df_unique$spread <- as.numeric(data_df_unique$spread)

plot(serie, xlab = "Dates", ylab = "Indice de production industrielle", main = "Indice")

serie_diff <- diff(serie, 1)

plot(cbind(serie, serie_diff), xlab = "Dates", ylab = "Indice de production industrielle", main = "Indice et différences premières")

plot(serie, xlab = "Dates", ylab = "Indice de production industrielle", main = "Indice")

serie_diff <- diff(serie, 1)
plot(cbind(serie, serie_diff), xlab = "Dates", ylab = "Indice de production industrielle", main = "Indice et différences premières")

### Q1 : Representation de la serie 
dates_char <- as.character(data$dates)
dates_char[1];dates_char[length(dates_char)] #affiche la premi`ere et la derni`ere date
dates <- as.yearmon(seq(from=2005+2/12,to=2022+10/12,by=1/12)) #index des dates pour spread
serie <- zoo(data,order.by=dates)
serie_diff <- diff(data$spread,1) #difference premiere
plot(serie,xlab="Dates",ylab="Indice de production industrielle",main="Indice")

plot(cbind(serie,serie_diff))

monthplot(serie)

acf(serie)
pacf(serie)

fit1<-decompose(serie)
plot(fit1)

#2 

summary(lm(serie~dates))
#regression lineaire simple de la serie chronologique s en fonction de sa position temporelle seq(1,n)
#calcul des coefficients de regression et des statistiques associees a cette relation lineaire



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

Qtests(adf@test$lm$residuals24,length(adf@test$lm$coefficients))

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
#la racine unitaire est rejetée au seuil de 1%


## Philippe perron
pp.test(serie) 
#pp test sur la serie (cas general : avec constante et tendance)


## test KPSS
kpss.test(serie,null="Trend") #prend en compte à la fois une tendance et une constante
#on rejette H0 (H0: la serie est stationnaire) a  5% et 10%
#la serie n'est donc pas stationnaire


#### Q2 #### Stationnariser la serie
serie_diff <- diff(serie)
plot(serie_diff, xaxt="s", type = "l")

summary(lm(serie_diff~dates[-1])) #on enlève la première date car la série est différenciée
#on trouve a nouveau une série avec tendance et constante non nulles (les deux coefficients sont significatifs)

#test de stationnarité (H0: n'est pas stationnaire ; H1: statio)
adf <- adfTest_valid(serie_diff,24,"ct")
adf #On rejette H0 sur tous les niveaux de confiance usuels 

pp.test(serie_diff) #Idem
kpss.test(serie_diff,null="Trend") 
#On ne rejette pas H0 (= Est stationnaire)

#Enlever la moyenne
s_centre <- serie_diff - mean(serie_diff)

#représentation des deux séries
plot(cbind(serie,serie_diff),main="Representation des deux series") 

acf(as.numeric(serie_diff) , main="ACF de s_diff")
pacf(as.numeric(serie_diff), main = "PACF de s_diff")
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




#### Q4 ####
#on cree une fonction pour automatiser le calcul de ratios
#fonction de test des significations individuelles des coefficients
#on recupere les coeffs
#"on recupere les std errors"
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

arima111<-arima(serie,c(1,1,1),include.mean=F)
arima111

##### Partie3:Previsions #####

#### Question 6

plot(density(arma11$residuals,lwd=0.5),xlim=c(-10,10),main="Densite des residus",xlab="Valeurs",ylab="Densite")

mu<-mean(arma11$residuals)
sigma<-sd(arma11$residuals)
x<-seq(-10,10)
y<-dnorm(x,mu,sigma)
lines(x,y,lwd=0.5,col="blue")

### Question 8 : Traçons la région de confiance pour la serie à 95%

library(forecast)
fore=forecast(arima111,h=3,level=95)
par(mfrow=c(1,1))
plot(fore,col=1,fcol=2,shaded=TRUE,xlab="Temps",ylab="Valeur",main="Prevision pour un ARIMA(1,1,1) avec une moyenne nulle")

#Ensuite,onrepresentelaregiondeconfiancebivarieea95%.require(ellipse)

XT1=predict(arma11,n.ahead=2)$pred[1]
XT2=predict(arma11,n.ahead=2)$pred[2]
XT1
XT2

require(ellipsis)
require(car)
require(ellipse)
library(ellipse)
arma=arima0(serie_diff,order=c(1,0,1))

arma11<-arima(serie_diff,c(1,0,1),include.mean=F)
arma11$coef
phi_1<-as.numeric(arma11$coef[1])
theta_1<-as.numeric(arma11$coef[2])
sigma2<-as.numeric(arma11$sigma2)

phi_1
theta_1
sigma2


Sigma<-matrix(c(sigma2,(phi_1+theta_1)*sigma2,(phi_1+theta_1)*sigma2,(1+(phi_1+theta_1)^2)*sigma2),ncol=2)

inv_Sigma<-solve(Sigma)

par(mar=c(5, 5, 4, 2) + 0.1)
plot(XT1, XT2, xlim=c(-10,10), ylim=c(-10,10), xlab="PrevisiondeX(T+1)", ylab="PrevisiondeX(T+2)", main="Regiondeconfiancebivariee 95%")

# Calcul des bornes de la region de confiance
conf <- ellipse(
  x = XT1, y = XT2, level = 0.95, scale = c(sqrt(sigma2), sqrt(sigma2)), 
  draw = TRUE, col = "red"
)
lines(conf, col="red")
