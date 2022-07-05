library(quantmod)
library(ghyp)
library(fGarch)
library(fBasics)
library(evd)
library(ggplot2)
library(copula)

#################################
###     Based on lectures     ###
### https://rpubs.com/iezepov ###
#################################


getSymbols("V",src='yahoo',from="2010-01-01")
V <-as.numeric(V$V.Close);
getSymbols("MA",src='yahoo',from="2010-01-01")
MA <-as.numeric(MA$MA.Close)

T=min(length(MA),length(V))

V <- V[2:T] / V[1:T-1] - 1
MA <- MA[2:T] / MA[1:T-1] - 1
######
# доходности портфеля из двух активов
prt <- cbind(V,MA);
tail(V)
# оценка параметров модели
prt.fit <- fit.hypmv(prt,symmetric=FALSE,silent=TRUE)
aic.mv <- stepAIC.ghyp(prt, dist = "ghyp")
# оценки риска
prt.fit <- fit.ghypmv(prt,symmetric=FALSE,silent=TRUE)
c1<-cor(1:length(MA),MA)
c2<-cor(1:length(V),V)
w <- c(c1/(c1+c2),c2/(c1+c2)) # веса активов в портфеле
####Parametrs####
alpha=.1;
N=1000;

sim <- rghyp(n=N,object=prt.fit)
prt.sim <- w[1]*sim[,1]+w[2]*sim[,2]
prt.sim <- sort(prt.sim)
VaR <- prt.sim[alpha*N]
ES <- mean(prt.sim[1:(alpha*N-1)])
VaR
ES
ES.ghyp<-ES
VaR.ghyp<-VaR
# выбор оптимальных весов активов в портфеле
opt <- portfolio.optimize(prt.fit, 
                          risk.measure="value.at.risk",
                          type="minimum.risk", 
                          level=0.95,silent=TRUE)
opt$opt.weights # искомые веса

# моделирование частных функций распределения
V.fit <- stepAIC.ghyp(V,dist=c("gauss","t","ghyp"),
                        symmetric=NULL,silent=TRUE)$best.model
MA.fit <- stepAIC.ghyp(MA,dist=c("gauss","t","ghyp"),
                        symmetric=NULL,silent=TRUE)$best.model
# расчёт значений F���1 и F���2
V.cdf <- pghyp(V,object=V.fit)
MA.cdf <- pghyp(MA,object=MA.fit)
cdf <- array(c(V.cdf,MA.cdf),dim=c(T-1,2))



# объявление копул
norm.cop <- normalCopula(dim=2,param=0.5,dispstr="un")
stud.cop <- tCopula(dim=2,param=0.5,df=5,
                    df.fixed=TRUE,dispstr="un")
gumb.cop <- gumbelCopula(dim=2,param=2)
clay.cop <- claytonCopula(dim=2,param=2)
# подгонка копулы
norm.fit <- fitCopula(cdf,copula=norm.cop)
stud.fit <- fitCopula(cdf,copula=stud.cop)
gumb.fit <- fitCopula(cdf,copula=gumb.cop)
clay.fit <- fitCopula(cdf,copula=clay.cop)


# значения частных функций распределения
N <- 10^4
stud.sim <- rCopula(n=N,copula=stud.fit@copula)
# доходности активов
V.sim <- qghyp(stud.sim[,1],object=V.fit)
MA.sim <- qghyp(stud.sim[,2],object=MA.fit)
w <- opt$opt.weights
prt.sim <- w[1]*V.sim + w[2]*MA.sim
# измерители риска
alpha <- 0.1
prt.sim <- sort(prt.sim)
VaR <- prt.sim[alpha*N]
ES <- mean(prt.sim[1:(alpha*N-1)])
VaR
ES
ES.copula<-ES
VaR.copula<-VaR

# одномерные GARCH-модели
library(fGarch)
V.gfit <- garchFit(data=V,formula=~garch(1,1),
                     shape=1.25,include.shape=F,cond.dist="ged",trace=F)
MA.gfit <- garchFit(data=MA,formula=~garch(1,1),
                     shape=1.3,include.shape=F,cond.dist="sged",trace=F)

# стандартизированные остатки

z <- matrix(nrow=T-1,ncol=2)
z[,1] <- V.gfit@residuals / V.gfit@sigma.t
z[,2] <- MA.gfit@residuals / MA.gfit@sigma.t
# частные распределения остатков
mean <- c(0,0); sd <- c(1,1); nu <- c(1.25,1.3)
xi <- c(1,MA.gfit@fit$par["skew"])
cdf <- matrix(nrow=T-1,ncol=2)
for (i in 1:2) cdf[,i] <- psged(z[,i],mean=mean[i],
                                sd=sd[i],nu=nu[i],xi=xi[i])

# подгонка копул
norm.fit <- fitCopula(cdf,copula=norm.cop)
stud.fit <- fitCopula(cdf,copula=stud.cop)
gumb.fit <- fitCopula(cdf,copula=gumb.cop)
clay.fit <- fitCopula(cdf,copula=clay.cop)
# метод Монте-Карло
cdf.sim <- rCopula(n=N,copula=stud.fit@copula)
z.sim <- matrix(nrow=N,ncol=2)
for (i in 1:2) z.sim[,i] <- qsged(cdf.sim[,i],
                                  mean=mean[i],sd=sd[i],nu=nu[i],xi=xi[i])
frc1 <- predict(V.gfit,n.ahead=1)
frc2 <- predict(MA.gfit,n.ahead=1)
mu <- c(frc1[,1],frc2[,1])
sigma <- c(frc1[,3],frc2[,3])

# модельные доходности портфеля
prt.sim <- w[1]*(mu[1]+sigma[1]*z.sim[,1]) +
  w[2]*(mu[2]+sigma[2]*z.sim[,2])
# измерители риска
prt.sim <- sort(prt.sim)
VaR <- prt.sim[alpha*N]
ES <- mean(prt.sim[1:(alpha*N-1)])

VaR
ES
ES.gcopula<-ES
VaR.gcopula<-VaR




# загрузка данных
Prtf <- cbind(V,MA)
# расчёт максим
m=T-1
Mn <- rep(0,times=m*2)
dim(Mn) <- c(m,2)
n=10
for (i in 1:2) {
  for (j in 1:(m/n))
    Mn[j,i] <- max(Prtf[((j-1)*n+1):(j*n),i])
}
Mn1<-Mn[1:(m/n),]
rm(Mn)
Mn<-Mn1
rm(Mn1)
# частные распределения на основе GED
fit1 <- fgev(Mn[,1])
fit2 <- fgev(Mn[,2])
# экстремальные копулы
library(copula)
gumb.cop <- gumbelCopula(2)
gal.cop <- galambosCopula(2)
# значения частных функций распределения
cdf1 <- pgev(Mn[,1],loc=fit1$estimate[1],
             scale=fit1$estimate[2],shape=fit1$estimate[3])
cdf2 <- pgev(Mn[,2],loc=fit2$estimate[1],
             scale=fit2$estimate[2],shape=fit2$estimate[3])
cdf <- cbind(cdf1,cdf2)
# подгонка копулы
gumb.fit <- fitCopula(cdf,copula=gumb.cop)
gal.fit <- fitCopula(cdf,copula=gal.cop)

# модельные значения максим
N <- 10^5
cdf.sim <- rCopula(n=N,copula=gal.fit@copula)
sim1 <- qgev(cdf.sim[,1],loc=fit1$estimate[1],
             scale=fit1$estimate[2],shape=fit1$estimate[3])
sim2 <- qgev(cdf.sim[,2],loc=fit2$estimate[1],
             scale=fit2$estimate[2],shape=fit2$estimate[3])

# модельные убытки
w <- opt$opt.weights
loss <- sort(w[1]*sim1+w[2]*sim2)
# расчёт мер риска
k <- 4
alpha <- 1-1/k
VaR <- loss[alpha*N]
ES <- mean(loss[(alpha*N+1):N])
VaR
ES
VaR.Mn <- VaR
ES.Mn <- ES


# выборка значений, превышающих многомерный порог
u <- c(sort(V)[0.9*(T-1)],sort(MA)[0.9*(T-1)])
t.Prtf <- Prtf[(Prtf[,1]>u[1])&(Prtf[,2]>u[2]),]
# частные распределения на основе GED
fit1 <- fpot(t.Prtf[,1],threshold=u[1],
             model="gpd",method="SANN")
fit2 <- fpot(t.Prtf[,2],threshold=u[2],
             model="gpd",method="SANN")
# значения частных функций распределения
cdf1 <- pgpd(t.Prtf[,1],loc=u[1],scale=fit1$par[1],
             shape=fit1$par[2])
cdf2 <- pgpd(t.Prtf[,2],loc=u[2],scale=fit2$par[1],
             shape=fit2$par[2])
cdf <- cbind(cdf1,cdf2)
# подгонка копулы
gumb.fit <- fitCopula(cdf,copula=gumb.cop)
gal.fit <- fitCopula(cdf,copula=gal.cop)

# модельные значения убытков
cdf.sim <- rCopula(n=N,copula=gal.fit@copula)
sim1 <- qgpd(cdf.sim[,1],loc=u[1],scale=fit1$par[1],
             shape=fit1$par[2])
sim2 <- qgpd(cdf.sim[,2],loc=u[2],scale=fit2$par[1],
             shape=fit2$par[2])

# убытки по портфелю
loss <- sort(w[1]*sim1+w[2]*sim2)
# расчёт мер риска
Fu <- nrow(t.Prtf)/(T-1)
alpha <- 1-1/(260*Fu)
VaR <- loss[alpha*N]
ES <- mean(loss[(alpha*N+1):N])

VaR
ES
ES.loss<-ES
VaR.loss<-VaR
##########
###VaRs###
##########
# разделим выборку на обучающую и экзаменующую
T=length(V)
T1 <- T-T2; T2 <- 20
# на пространстве экзаменующей выборки построим набор
# последовательных значений VaR
VaR.l.ghyp <- numeric()
VaR.l.copula <- numeric()
VaR.l.gcopula <- numeric()
VaR.l.Mn <- numeric()
VaR.l.loss <- numeric()

ES.l.ghyp <- numeric()
ES.l.copula <- numeric()
ES.l.gcopula <- numeric()
ES.l.Mn <- numeric()
ES.l.loss <- numeric()

#Объявление Копул
norm.cop <- normalCopula(dim=2,param=0.5,dispstr="un")
stud.cop <- tCopula(dim=2,param=0.5,df=5,
                    df.fixed=TRUE,dispstr="un")
gumb.cop <- gumbelCopula(dim=2,param=2)
clay.cop <- claytonCopula(dim=2,param=2)

###############
###Параметры###
###############
alpha=.1;
h <- 560; # длина обучающей выборки
N <- 10^3
a=10
k <- 4
##########
###Цикл###
##########
for (i in (T1+1):(T1+T2)) 
{
  h.prt <- prt[(i-h):(i-1),]
  h.V <- V[(i-h):(i-1)]
  h.MA <- MA[(i-h):(i-1)]
  prt.fit <- fit.hypmv(h.prt,symmetric=FALSE,silent=TRUE)
  aic.mv <- stepAIC.ghyp(h.prt, dist = "ghyp")
  # выбор оптимальных весов активов в портфеле
  opt <- portfolio.optimize(prt.fit, 
                            risk.measure="value.at.risk",
                            type="minimum.risk", 
                            level=0.95,silent=TRUE)
  w <- opt$opt.weights # искомые веса
  # оценки риска
  prt.fit <- fit.ghypmv(h.prt,symmetric=FALSE,silent=TRUE)
  sim <- rghyp(n=N,object=prt.fit)
  prt.sim <- w[1]*sim[,1]+w[2]*sim[,2]
  prt.sim <- sort(prt.sim)
  VaR.l.ghyp[i-T1] <- prt.sim[alpha*N]
  ES.l.ghyp[i-T1] <- mean(prt.sim[1:(alpha*N-1)])
  # моделирование частных функций распределения
  V.fit <- stepAIC.ghyp(h.V,dist=c("gauss","t","ghyp"),
                        symmetric=NULL,silent=TRUE)$best.model
  MA.fit <- stepAIC.ghyp(h.MA,dist=c("gauss","t","ghyp"),
                         symmetric=NULL,silent=TRUE)$best.model
  # расчёт значений F���1 и F���2
  V.cdf <- pghyp(V,object=V.fit)
  MA.cdf <- pghyp(MA,object=MA.fit)
  cdf <- array(c(V.cdf,MA.cdf),dim=c(h,2))
  # подгонка копулы
  norm.fit <- fitCopula(cdf,copula=norm.cop)
  stud.fit <- fitCopula(cdf,copula=stud.cop)
  gumb.fit <- fitCopula(cdf,copula=gumb.cop)
  clay.fit <- fitCopula(cdf,copula=clay.cop)
  
  
  # значения частных функций распределения
  stud.sim <- rCopula(n=N,copula=stud.fit@copula)
  # доходности активов
  V.sim <- qghyp(stud.sim[,1],object=V.fit)
  MA.sim <- qghyp(stud.sim[,2],object=MA.fit)
  prt.sim <- w[1]*V.sim + w[2]*MA.sim
  # измерители риска
  prt.sim <- sort(prt.sim)
  VaR.l.copula[i-T1] <- prt.sim[alpha*N]
  ES.l.copula[i-T1] <- mean(prt.sim[1:(alpha*N-1)])
  
  
  # одномерные GARCH-модели
  V.gfit <- garchFit(data=h.V,formula=~garch(1,1),
                     shape=1.25,include.shape=F,cond.dist="ged",trace=F)
  MA.gfit <- garchFit(data=h.MA,formula=~garch(1,1),
                      shape=1.3,include.shape=F,cond.dist="sged",trace=F)
  
  # стандартизированные остатки
  
  z <- matrix(nrow=h,ncol=2)
  z[,1] <- V.gfit@residuals / V.gfit@sigma.t
  z[,2] <- MA.gfit@residuals / MA.gfit@sigma.t
  # частные распределения остатков
  mean <- c(0,0); sd <- c(1,1); nu <- c(1.25,1.3)
  xi <- c(1,MA.gfit@fit$par["skew"])
  cdf <- matrix(nrow=h,ncol=2)
  for (i in 1:2) cdf[,i] <- psged(z[,i],mean=mean[i],
                                  sd=sd[i],nu=nu[i],xi=xi[i])
  
  # подгонка копул
  norm.fit <- fitCopula(cdf,copula=norm.cop)
  stud.fit <- fitCopula(cdf,copula=stud.cop)
  gumb.fit <- fitCopula(cdf,copula=gumb.cop)
  clay.fit <- fitCopula(cdf,copula=clay.cop)
  # метод Монте-Карло
  cdf.sim <- rCopula(n=N,copula=stud.fit@copula)
  z.sim <- matrix(nrow=N,ncol=2)
  for (i in 1:2) z.sim[,i] <- qsged(cdf.sim[,i],
                                    mean=mean[i],sd=sd[i],nu=nu[i],xi=xi[i])
  frc1 <- predict(V.gfit,n.ahead=1)
  frc2 <- predict(MA.gfit,n.ahead=1)
  mu <- c(frc1[,1],frc2[,1])
  sigma <- c(frc1[,3],frc2[,3])
  
  # модельные доходности портфеля
  prt.sim <- w[1]*(mu[1]+sigma[1]*z.sim[,1]) +
    w[2]*(mu[2]+sigma[2]*z.sim[,2])
  # измерители риска
  prt.sim <- sort(prt.sim)
  VaR.l.gcopula[i-T1] <- prt.sim[alpha*N]
  ES.l.gcopula[i-T1] <- mean(prt.sim[1:(alpha*N-1)])
  
  # Данные
  Prtf <- cbind(h.V,h.MA)
  # расчёт максим
  Mn <- rep(0,times=h*2)
  dim(Mn) <- c(h,2)
  n=a
  for (i in 1:2) {
    for (j in 1:(h/n))
      Mn[j,i] <- max(-Prtf[((j-1)*n+1):(j*n),i])
  }
  Mn1<-Mn[1:(h/n),]
  rm(Mn)
  Mn<-Mn1
  rm(Mn1)
  # частные распределения на основе GED
  fit1 <- fgev(Mn[,1])
  fit2 <- fgev(Mn[,2])

  # значения частных функций распределения
  cdf1 <- pgev(Mn[,1],loc=fit1$estimate[1],
               scale=fit1$estimate[2],shape=fit1$estimate[3])
  cdf2 <- pgev(Mn[,2],loc=fit2$estimate[1],
               scale=fit2$estimate[2],shape=fit2$estimate[3])
  cdf <- cbind(cdf1,cdf2)
  # подгонка копулы
  gumb.fit <- fitCopula(cdf,copula=gumb.cop)
  gal.fit <- fitCopula(cdf,copula=gal.cop)
  
  # модельные значения максим
  
  cdf.sim <- rCopula(n=N,copula=gal.fit@copula)
  sim1 <- qgev(cdf.sim[,1],loc=fit1$estimate[1],
               scale=fit1$estimate[2],shape=fit1$estimate[3])
  sim2 <- qgev(cdf.sim[,2],loc=fit2$estimate[1],
               scale=fit2$estimate[2],shape=fit2$estimate[3])
  
  # модельные убытки
  loss <- sort(w[1]*sim1+w[2]*sim2)
  # расчёт мер риска
  alpha <- 1-1/k
  VaR.l.Mn[i-T1] <- -loss[alpha*N]
  ES.l.Mn[i-T1] <- mean(-loss[(alpha*N+1):N])
  
  
  
  # выборка значений, превышающих многомерный порог
  u <- c(sort(-h.V)[0.9*h],sort(-h.MA)[0.9*h])
  t.Prtf <- Prtf[(Prtf[,1]>u[1])&(Prtf[,2]>u[2]),]
  # частные распределения на основе GED
  fit1 <- fpot(t.Prtf[,1],threshold=u[1],
               model="gpd",method="SANN")
  fit2 <- fpot(t.Prtf[,2],threshold=u[2],
               model="gpd",method="SANN")
  # значения частных функций распределения
  cdf1 <- pgpd(t.Prtf[,1],loc=u[1],scale=fit1$par[1],
               shape=fit1$par[2])
  cdf2 <- pgpd(t.Prtf[,2],loc=u[2],scale=fit2$par[1],
               shape=fit2$par[2])
  cdf <- cbind(cdf1,cdf2)
  # подгонка копулы
  gumb.fit <- fitCopula(cdf,copula=gumb.cop)
  gal.fit <- fitCopula(cdf,copula=gal.cop)
  
  # модельные значения убытков
  cdf.sim <- rCopula(n=N,copula=gal.fit@copula)
  sim1 <- qgpd(cdf.sim[,1],loc=u[1],scale=fit1$par[1],
               shape=fit1$par[2])
  sim2 <- qgpd(cdf.sim[,2],loc=u[2],scale=fit2$par[1],
               shape=fit2$par[2])
  
  # убытки по портфелю
  loss <- sort(w[1]*sim1+w[2]*sim2)
  # расчёт мер риска
  Fu <- nrow(t.Prtf)/T
  alpha <- 1-1/(260*Fu)
  VaR.l.loss[i-T1] <- -loss[alpha*N]
  ES.l.loss[i-T1] <- mean(-loss[(alpha*N+1):N])
}

  fact <- cbind(V[(T1+1):(T1+T2)],MA[(T1+1):(T1+T2)])
  prt.fit <- fit.hypmv(fact,symmetric=FALSE,silent=TRUE)

# выбор оптимальных весов активов в портфеле
  opt <- portfolio.optimize(prt.fit, 
                          risk.measure="value.at.risk",
                          type="minimum.risk", 
                          level=0.95,silent=TRUE)
  w <- opt$opt.weights # искомые веса
  factR <- w[1]*V[(T1+1):(T1+T2)]+w[2]*MA[(T1+1):(T1+T2)]

#################
###Kupic tests###
#################
#For Ghyp
K <- sum(factR<VaR.l.ghyp); 
alpha0 <- K/T2
S <- -2*log((1-alpha)^(T2-K)*alpha^K)+
  2*log((1-alpha0)^(T2-K)*alpha0^K)
p.value.ghyp <- 1-pchisq(S,df=1)
#For Copula
K <- sum(factR<VaR.l.copula); 
alpha0 <- K/T2
S <- -2*log((1-alpha)^(T2-K)*alpha^K)+
  2*log((1-alpha0)^(T2-K)*alpha0^K)
p.value.copula <- 1-pchisq(S,df=1)
#For gCopula
K <- sum(factR<VaR.l.gcopula); 
alpha0 <- K/T2
S <- -2*log((1-alpha)^(T2-K)*alpha^K)+
  2*log((1-alpha0)^(T2-K)*alpha0^K)
p.value.gcopula <- 1-pchisq(S,df=1)
#For Mn
K <- sum(factR<VaR.l.Mn); 
alpha0 <- K/T2
S <- -2*log((1-alpha)^(T2-K)*alpha^K)+
  2*log((1-alpha0)^(T2-K)*alpha0^K)
p.value.Mn <- 1-pchisq(S,df=1)
#For loss
K <- sum(factR<VaR.l.loss); 
alpha0 <- K/T2
S <- -2*log((1-alpha)^(T2-K)*alpha^K)+
  2*log((1-alpha0)^(T2-K)*alpha0^K)
p.value.loss <- 1-pchisq(S,df=1)

p.value.ghyp
p.value.copula
p.value.gcopula
p.value.Mn
p.value.loss

#############
###Графики###
#############
#For Ghyp
plot(factR,type="l")
lines(VaR.l.ghyp,col="red")
#For Copula
plot(factR,type="l")
lines(VaR.l.copula,col="red")
#For gCopula
plot(factR,type="l")
lines(VaR.l.gcopula,col="red")
#For Mn
plot(factR,type="l")
lines(VaR.l.Mn,col="red")
#For loss
plot(factR,type="l")
lines(VaR.l.loss,col="red")
#############
###Таблица###
#############
  table<-array(c(
    VaR.ghyp, VaR.copula, VaR.gcopula, VaR.Mn, VaR.loss,
    ES.ghyp, ES.copula, ES.gcopula, ES.Mn, ES.loss,
    p.value.ghyp, p.value.copula, p.value.gcopula, p.value.Mn, p.value.loss
    ), dim=c(3,5))
  library(knitr)
  kable(table, format = "html", digits = 4, 
        row.names = c("Value at Risk","ES","P.Value in Kupic test"),
        col.names = c("Ghyp dist.",
                    "Copula based on ghyp",
                    "Copula based on Garch",
                    "Method of maxim",
                    "Method of loss"), align=c("c"), caption =  "Сводная таблица")