install.packages('zoo')
install.packages('tseries')
install.packages('fUnitRoots')
require(zoo) #format de serie temporelle pratique et facile d'utilisation (mais plus volumineux)
require(tseries) #diverses fonctions sur les series temporelles
require(fUnitRoots)

###### Partie 1 : Les donnees ######

#Extraction des donnees : 
path <- "C:/Users/charl/Documents/ENSAE/S2/Séries temporelles" 
setwd(path) 

#Mise en forme : 
datafile <- "huiles.csv" 
data <- read.csv(datafile,sep=";")

xm <- zoo(data[[2]]) #convertit le premier element de data en serie temporelle de type "zoo"
T <- length(xm)

#xm <- xm.source[1:(T-4)] #supprime les 4 dernieres valeurs
#class <- class(data$indice)
#head(xm, n=10)

# Q1 : Representation de la serie 

plot(xm, xaxt="s") 
axis(side=1,at=seq(1990,2022,2)) 

#acf(xm) #on retire peut etre après 

#lag.plot(xm,lags=12,layout=c(3,4),do.lines = FALSE )
#graphique qui montre les corrÃ©lations entre la sÃ©rie xm et 
#sa propre version dÃ©calÃ©e (lagged version) jusqu'Ã  un maximum de 12 dÃ©calages, 
#disposÃ©s en une grille de 3 lignes et 4 colonnes, sans afficher les lignes de corrÃ©lation.


#la tendance
summary(lm(xm~seq(1,n)))

## test de stationnarité
# Philippe peron
pp.test(xm) #pp test sur la série (cas g´en´eral : avec constante et tendance)

# dickey fuller
adf.test(xm)
# On ne rejette pas H0 

#testKPSS
kpss.test(xm,null="Trend") #on rejette H0 (H0: la sÃ©rie est stationnaire) Ã  5% et 10%
#la sÃ©rie n'est donc pas stationnaire


#### Q2 #### Stationnariser la serie
xm_diff <- diff(xm)
plot(xm_diff, xaxt="s")

#test de stationnarité H0: n'est pas stationnaire H1: statio
adf.test(xm_diff)
pp.test(xm_diff)
kpss.test(xm_diff,null="Trend")

#Enlever la moyenne
xm_centre <- xm_diff - mean(xm_diff)

# statio centré donc on peut mettre ARMA
par(mfrow = c(1,2))
acf(xm_centre,24);pacf(xm_centre,24)

# ARMA d'ordres 4,1 à vérifier 


#Partie 2 : Modèles ARMA


# On vérifie avec le test de JungBox

# Validation du modèle

Qtests <- function(series, k, fitdf=0) {
  pvals <- apply(matrix(1:k), 1, FUN=function(l) {
    pval <- if (l<=fitdf) NA else Box.test(series, lag=l, type="Ljung-Box", fitdf=fitdf)$p.value
    return(c("lag"=l,"pval"=pval))
  })
  return(t(pvals))
}
#tests de LB pour les ordres 1 a 24

# L’absence d’autocorrelation n’est jamais rejetee `a 95% jusqu’`a 24 retards. Le mod`ele est donc valide
# dans la fonction Qtestes: param?tre num?ro 1: series
# le test de Ljung-Box ne rejette pas l'absence d'autocorrelation des r?sidus ? l'ordre 6
# matrice de pvals de k lignes et 1 colonnes: remplie de 1, pour chaque ?l?ment on applique la fonction (if... on met NA, sinon ... quand on sort de la fonction on met return)
#on trouve que pour 1 ? 5: NA: le k est plus petit 
#dans ce cas on veut accepter H0 (r?sidus non corr?l?s) donc on veut des p-valeurs ?lev?es: on est ravis ici



#les coeff de l'ar3, on fait le ratio de 0.1748/0.1604 = 1.08 < 1.96 --> coef non significatif, on peut simplifier le mod?le
#On teste tous les mod?les pour trouver les mod?les bien ajust?s et valides, on les compare tous 



summary(lm(spread ~ dates))
#on fait une regression : les p-valeurs sont tres faibles, coef significatif au seuil de 1% (on le sait aussi car il y a 3 ?toiles)
#onfait le test de DF avec c et t : comme les deux coeff sont significatifs au seuil de 1% il sont non nuls donc on les inclus. 
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


#### Q6 ####
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

estim <- arima(y,c(4,0,1)); arimafit(estim)
# pas bien ajusté

estim <- arima(y,c(3,0,1)); arimafit(estim)
# Apas bien ajusté

estim <- arima(y,c(2,0,1)); arimafit(estim)
# pas bien ajuste 


estim <- arima(y,c(1,0,1)); arimafit(estim)
# bien ajusté et valide : youpi
arma11 <- estim

estim <- arima(y,c(0,0,1)); arimafit(estim)
# bien ajuste et pas du tout valide

estim <- arima(y,c(4,0,0)); arimafit(estim)
# bien ajuste et valide : cool 
ar4 <- estim

estim <- arima(y,c(3,0,0)); arimafit(estim)
# bien ajuste et pas valide

estim <- arima(y,c(2,0,0)); arimafit(estim)
# bien ajuste mais pas valide

estim <- arima(y,c(1,0,0)); arimafit(estim)
# bien ajuste mais pas valide


# on a selectionne 3 modeles ar1ma1 et ar4 valides et bien ajustes. On regarde les test BIC et AIC pour selectionner le meilleur modele.
models <- c("arma11","ar4"); names(models) <- models
apply(as.matrix(models),1, function(m) c("AIC"=AIC(get(m)), "BIC"=BIC(get(m))))
#on regarde en 1 AIC et BIC cr??s pour comparer des mod?les; ensuite on regarde le R?
#Si AIC et BIC s?lectionnent 2 mod?les, ils restent deux mod?les candidats que l'on d?partage par le R?

#on prend le ma2 car plus petits AIC et BIC. on regarde ensuite R? pour confirmer notre intuition

#une fois que l'on a le mod?le, on continue vace le meilleur mod?le. Ici on garde les autres candidats pour la p?dagogie

#### Q7 ####

##cr?ation de s?ries o? ? chaque colonne sera assign?e la pr?diction par un mod?le
models <-  c("ar3","ma2","ar2ma1") #c vecteur
preds <- zoo(matrix(NA,ncol=3,nrow=4),order.by=tail(index(xm.source),4)) #4lignes car on pr?voit sur 4 horizons
colnames(preds) <- models #on met les noms des mod?les: 1er nom de colonne: ar3, 2?me nom de col ma2
desaisonp <- preds #on met une s?rie vierge dans laquelle on stockera les pr?visions
xmp <- preds #

##
for (m in models){
  pred1 <- mean(desaison) + zoo(predict(get(m),4)$pred, order.by=tail(index(xm.source),4)) #dans pred 1 on met la valeur, 
  pred2 <- as.numeric(tail(xm,12))[1:4] + pred1 
  desaisonp[,m] <- pred1
  xmp[,m] <- pred2
}

obs <- tail(xm.source,4) #
cbind(obs,xmp) #
apply(xmp,2, function(x) sqrt(sum((x-obs)^2)/4)/sd(xm.source)) 
# le ma2 est le plus proche du nb d'observations
#on voit que le arma21 performe le mieux car avec la plus petite mean square erreur (une sorte de R? non ajust?)

#### Q8 ####
datafile <- "Donnees2.csv" #definit le fichier de donnees

data <- read.csv(datafile,sep=";") #importe un fichier .csv dans un objet de classe data.frame
xm.source <- zoo(data[[1]]) #convertit le premier element de data en serie temporelle de type "zoo"
T <- length(xm.source)
xm <- xm.source[1:(T-4)] #supprime les 4 dernieres valeurs
dev.off() #reinitialise les parametre de graphique
plot(xm)
### 

trend <- 1:length(xm)
lt <- lm(xm ~ trend) #
summary(lt) #
r <- lt$residuals #
par(mfrow=c(1,2))
plot(r)
acf(r)
### 

pp.test(xm) 
### 

acf(r,24);pacf(r,24) 
### 
### 
pmax=4; qmax=21

### 



## fonction pour estimer un arima et en verifier l'ajustement et la validite
modelchoice <- function(p,q,data=r, k=24){
  estim <- try(arima(data, c(p,0,q),optim.control=list(maxit=20000)))
  if (class(estim)=="try-error") return(c("p"=p,"q"=q,"arsignif"=NA,"masignif"=NA,"resnocorr"=NA, "ok"=NA))
  arsignif <- if (p==0) NA else signif(estim)[3,p]<=0.05
  masignif <- if (q==0) NA else signif(estim)[3,p+q]<=0.05
  resnocorr <- sum(Qtests(estim$residuals,24,length(estim$coef)-1)[,2]<=0.05,na.rm=T)==0
  checks <- c(arsignif,masignif,resnocorr)
  ok <- as.numeric(sum(checks,na.rm=T)==(3-sum(is.na(checks))))
  return(c("p"=p,"q"=q,"arsignif"=arsignif,"masignif"=masignif,"resnocorr"=resnocorr,"ok"=ok))
}

## fonction pour estimer et verifier tous les arima(p,q) avec p<=pmax et q<=max
armamodelchoice <- function(pmax,qmax){
  pqs <- expand.grid(0:pmax,0:qmax)
  t(apply(matrix(1:dim(pqs)[1]),1,function(row) {
    p <- pqs[row,1]; q <- pqs[row,2]
    cat(paste0("Computing ARMA(",p,",",q,") \n"))
    modelchoice(p,q)
  }))
}

armamodels <- armamodelchoice(pmax,qmax) #estime tous les arima (patienter...)


selec <- armamodels[armamodels[,"ok"]==1&!is.na(armamodels[,"ok"]),] #modeles bien ajustes et valides
selec
### On a ? modeles bien ajustes et valides

pqs <- apply(selec,1,function(row) list("p"=as.numeric(row[1]),"q"=as.numeric(row[2]))) #cree une liste des ordres p et q des modeles candidats
names(pqs) <- paste0("arma(",selec[,1],",",selec[,2],")") #renomme les elements de la liste
models <- lapply(pqs, function(pq) arima(r,c(pq[["p"]],0,pq[["q"]]))) #cree une liste des modeles candidats estimes
vapply(models, FUN.VALUE=numeric(2), function(m) c("AIC"=AIC(m),"BIC"=BIC(m))) #calcule les AIC et BIC des modeles candidats
### L'ARMA(?,?) minimise les criteres d'information.

rps <- lapply(models, function(m) as.zoo(predict(m,4)$pred)) #previsions de r
xmps <- lapply(rps, function(rp) rp+cbind(1,c((T-3):T))%*%lt$coefficients) #previsions de xm
rmse <- vapply(xmps, FUN.VALUE=numeric(1), function(xmp) sqrt(sum((as.zoo(xmp)-tail(xm.source,4))^2))) #calcule les rmse out-of-sample
rmse
### L'ARMA(?,?) fait aussi la meilleure prevision