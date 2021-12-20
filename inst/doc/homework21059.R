## ----figure-------------------------------------------------------------------
x=rnorm(1000)
hist(x,prob=T,main="Normal Distribution")
curve(dnorm(x),add=T)

## ----text---------------------------------------------------------------------
name<-"Vivian";m<-7;n<-17;k<-27
ls(pat="n")

## ----table--------------------------------------------------------------------
Table1<-data.frame('样本编码'=c('A','B','C'),'2020'=c('1.57','1.65','1.78'),'2021'=c('1.58','1.68','1.80'),check.names=FALSE)
Table1

## ----fig.width=8,fig.height=4-------------------------------------------------
n<-100
u<-runif(n)
sigma1<-1
x1<-{-2*sigma1^2*log(1-u)}^{1/2}
y1<-seq(0,100,.01)
z1<-exp(-(y1^2/(2*sigma1^2)))
hist(x1,prob=TRUE,main=expression(f(x)==e^(-frac(x^2,2*sigma^2))))
lines(y1,z1)
sigma2<-2
x2<-{-2*sigma2^2*log(1-u)}^{1/2}
y2<-seq(0,100,.01)
z2<-exp(-(y1^2/(2*sigma2^2)))
hist(x2,prob=TRUE,main=expression(f(x)==e^(-frac(x^2,2*sigma^2))))
lines(y2,z2)

## -----------------------------------------------------------------------------
n<-1e2
x1<-rnorm(n)
x2<-rnorm(n,3,1)
p1<-0.75
z1<-p1*x1+(1-p1)*x2
hist(z1,prob=TRUE,breaks=35,main=expression("p=0.75"))
p2<-0.5
z2<-p2*x1+(1-p2)*x2
hist(z2,prob=TRUE,breaks=35,main=expression("p=0.5"))
p3<-0.25
z3<-p3*x1+(1-p3)*x2
hist(z3,prob=TRUE,breaks=35,main=expression("p=0.25"))
#根据不同的p值所做直方图，当p=0.5时会出现双峰

## -----------------------------------------------------------------------------
n<-1e2
lambda1<-5
t<-10
Nt1<-rpois(n,lambda1*t)
afa<-1
beta<-1
#此处将Yi的概率密度参数设置为Ga（1，1）
Xt1<-rgamma(n,afa*Nt1,beta)
#由于Yi独立同分布，其和服从Ga（Nt1，1）
hist(Xt1,prob=TRUE)
EX1<-sum(Xt1)*{1/n}
VarX1<-sum((Xt1-EX1)^2)*{1/(n-1)}
#此处计算出由随机数生成的样本观测值的均值和方差
TRUEEX1<-t*lambda1*afa/beta
TRUEVarX1<-t*lambda1*(afa^2+afa)/beta
#此处计算出由公式得出的理论均值和方差
EX1/TRUEEX1
VarX1/TRUEVarX1
#通过相除得出λ=5时，理论与实际比值在1附近
lambda2<-10
Nt2<-rpois(n,lambda2*t)
Xt2<-rgamma(n,afa*Nt2,beta)
hist(Xt2,prob=TRUE)
EX2<-sum(Xt2)*{1/n}
VarX2<-sum((Xt2-EX2)^2)*{1/(n-1)}
TRUEEX2<-t*lambda2*afa/beta
TRUEVarX2<-t*lambda2*(afa^2+afa)/beta
EX2/TRUEEX2
VarX2/TRUEVarX2
#此处取λ=10，比值依旧在1附近
lambda3<-15
Nt3<-rpois(n,lambda3*t)
Xt3<-rgamma(n,afa*Nt3,beta)
hist(Xt3,prob=TRUE)
EX3<-sum(Xt3)*{1/n}
VarX3<-sum((Xt3-EX3)^2)*{1/(n-1)}
TRUEEX3<-t*lambda3*afa/beta
TRUEVarX3<-t*lambda3*(afa^2+afa)/beta
EX3/TRUEEX3
VarX3/TRUEVarX3
#此处取λ=15，比值依旧在1附近

## -----------------------------------------------------------------------------
x <- seq(.1, .9, length = 9)
n <- 100
u <- runif(n)
cdf <- numeric(length(x))
#由于原积分的区间为（0，x），运用变换t=ux将上下限变为（0，1）
for (i in 1:length(x)) {
     g <- 30*x[i]^3*u^2*(1-u*x[i])^2
     cdf[i] <- mean(g)
}
realvalue<- pbeta(x,3,3)
rbind(x, cdf, realvalue)
#结果表明用蒙特卡洛积分所得结果与真实值非常相近

## -----------------------------------------------------------------------------
MC.Phi <- function(x, sigma = 5, R = 1e4, antithetic = TRUE) {
  u <- runif(R/2)
  if (!antithetic) 
    v <- runif(R/2) 
  else
    v <- 1 - u
  u <- c(u, v)
  cdf <- numeric(length(x))
  for (i in 1:length(x)) {
    g <- x[i]^2*u/sigma^2*exp (-x[i]^2*u^2/(2*sigma^2))
    cdf[i] <- mean(g)
  }
  cdf
}

x <- seq(.1, 2.5, length=5)
MC1 <- MC.Phi(x, anti = FALSE)
MC2 <- MC.Phi(x)

m <- 100
MC1 <- MC2 <- numeric(m)
x <- 2
for (i in 1:m) {
MC1[i] <- MC.Phi(x, R = 1e4, anti = FALSE)
MC2[i] <- MC.Phi(x, R = 1e4)
}
print(sd(MC1))
print(sd(MC2))
print((var(MC1) - var(MC2))/var(MC1))
#由此看到方差降低了99%，效果较好

## -----------------------------------------------------------------------------
m <- 100
theta.hat <- se <- numeric(3)
#令t=1/x^2,做积分变换，将上下限变为（0，1），同时定义g（x）
g <- function(x) {
  (1/(2*sqrt(2*pi)*x^(5/2))) * exp(-1/(2*x)) * (x > 0) * (x < 1)
}

x <- runif(m) #f0=1,相当于不使用importance function
f0g <- g(x)
theta.hat[1] <- mean(f0g)
se[1] <- sd(f0g)

u <- runif(m) #f1=exp(-1)/(1 − e^(−1)), inverse transform method
x <- - log(1 - u * (1 - exp(-1)))
f1g <- g(x) / (exp(-x) / (1 - exp(-1)))
theta.hat[2] <- mean(f1g)
se[2] <- sd(f1g)

u <- runif(m) #f2= 4/(pi*(1+x^2)), inverse transform method
x <- tan(pi * u / 4)
f2g <- g(x) / (4 / ((1 + x^2) * pi))
theta.hat[3] <- mean(f2g)
se[3] <- sd(f2g)

rbind(theta.hat, se)
#结果表明方差降低，方法有效

## -----------------------------------------------------------------------------
m <- 100
theta.hat <- se <- numeric(1)
g <- function(x) {
  (1/(2*sqrt(2*pi)*x^(5/2))) * exp(-1/(2*x)) * (x > 0) * (x < 1)
}
u <- runif(m) #f1=exp(-1)/(1 − e^(−1)), inverse transform method，同上题
x <- - log(1 - u * (1 - exp(-1)))
fg <- g(x) / (exp(-x) / (1 - exp(-1)))
theta.hat[1] <- mean(fg)
se[1] <- sd(fg)
rbind(theta.hat, se)

## -----------------------------------------------------------------------------
n<-20
alpha<-0.05
UCL<-replicate(1e4,expr = {
  x<-rchisq(n,df=2)
  mean(x)+qt(1-alpha/2,n-1)*sd(x)/sqrt(n)
})
#UCL全称为upper confidence level,计算的是置信区间上界
set.seed(123)
LowerCL<-replicate(1e4,expr = {
  x<-rchisq(n,df=2)
  mean(x)-qt(1-alpha/2,n-1)*sd(x)/sqrt(n)
})
#LowerCL全称为lower confidence level,计算的是置信区间下界
sum(LowerCL<2 & UCL>2)
mean(LowerCL<2 & UCL>2)
#计算同时符合上界大于2、下界小于2（2在置信区间内）的个数及其概率，即coverage probability
#本题是一个针对均值进行的双侧区间估计，最终结果表明用t检验模拟卡方分布较为稳健，而Example6.4是一个针对方差进行的单侧区间估计，用t检验进行模拟的误差较大。

## -----------------------------------------------------------------------------
n <- 20
alpha <- .05
mu1 <- 1
sigma1 <- sqrt(2)#此处为第一小问，X服从卡方1时期望为1，标准差为根号2
m <- 1e4 
p <- numeric(m) 
for (j in 1:m) {
    x <- rchisq(n, mu1)  #此处生成服从卡方1的随机数
    ttest <- t.test(x, alternative = "two.sided", mu = mu1)#题中要求为双侧检验，下同
    p[j] <- ttest$p.value
    }
p.hat1 <- mean(p < alpha)
se.hat1 <- sqrt(p.hat1 * (1 - p.hat1) / m)
print(c(p.hat1, se.hat1))

mu2 <- 1
sigma2 <- sqrt(1/3)  #此处为第二小问，X服从U(0,2)时，期望为1，标准差为根号下1/3
m <- 1e4 
p <- numeric(m) 
for (j in 1:m) {
    x <- runif(n, 0,2) #此处生成服从U(0,2)的随机数样本
    ttest <- t.test(x, alternative = "two.sided", mu = mu2)
    p[j] <- ttest$p.value
    }
p.hat2 <- mean(p < alpha)
se.hat2 <- sqrt(p.hat2 * (1 - p.hat2) / m)
print(c(p.hat2, se.hat2))

mu3 <- 1
sigma3 <- 1  #此处为第三小问，X服从E(1)时，期望和标准差均为1
m <- 1e4 
p <- numeric(m) 
for (j in 1:m) {
    x <- rexp(n) #此处生成服从E(1)的随机数样本
    ttest <- t.test(x, alternative = "two.sided", mu = mu3)
    p[j] <- ttest$p.value
    }
p.hat3 <- mean(p < alpha)
se.hat3 <- sqrt(p.hat3 * (1 - p.hat3) / m)
print(c(p.hat3, se.hat3))

## -----------------------------------------------------------------------------
library(MASS)
#编写一个函数Mardiatest
Mardiatest<-function(data){
  n=nrow(data)
  c=ncol(data)
  x<-data
  for(i in 1:c){
    x[,i]<-data[,i]-mean(data[,i])
  }
  sigmahat<-t(x)%*%x/n
  a<-x%*%solve(sigmahat)%*%t(x)
  b<-sum(colSums(a^3))/(n^2)#计算出检验统计量
  test<-n*b/6
  cv<-qchisq(0.95,c*(c+1)*(c+2)/6)  #调用卡方分布的临界值
  as.integer(test>cv)  #通过取整函数，得到拒绝原假设的个数
}

set.seed(123)
mu <- c(0,0,0)
sigma <- matrix(c(1,0,0,0,1,0,0,0,1),nrow=3,ncol=3) 
m=1000
n<-c(10, 20, 30, 50, 100, 500)
#m: 重复次数; n: 样本容量
a=numeric(length(n))
for(i in 1:length(n)){
  a[i]=mean(replicate(m, expr={
    data <- mvrnorm(n[i],mu,sigma) #调用服从多元正态分布的随机数，生成矩阵data
    Mardiatest(data)
  }))
}
print(a)

## -----------------------------------------------------------------------------
library(MASS)
#把之前定义的函数再复制一遍
Mardiatest<-function(data){
  n=nrow(data)
  c=ncol(data)
  x<-data
  for(i in 1:c){
    x[,i]<-data[,i]-mean(data[,i])
  }
  sigmahat<-t(x)%*%x/n
  a<-x%*%solve(sigmahat)%*%t(x)
  b<-sum(colSums(a^3))/(n^2)
  test<-n*b/6
  cv<-qchisq(0.95,c*(c+1)*(c+2)/6)
  as.integer(test>cv)
}
#输入两个正态分布的参数
mu1 <- mu2 <- c(0,0,0)
sigma1 <- matrix(c(1,0,0,0,1,0,0,0,1),nrow=3,ncol=3)
sigma2 <- matrix(c(100,0,0,0,100,0,0,0,100),nrow=3,ncol=3)
alpha <- .1
n <- 30
m <- 2500
epsilon <- c(seq(0, .15, .01), seq(.15, 1, .05))
N <- length(epsilon)
pwr <- numeric(N)
for (j in 1:N) { #for each epsilon
    e <- epsilon[j]
    sktests <- numeric(m)
    for (i in 1:m) { #每次重复，以概率1-e和e抽取两个正态分布的随机样本
       index=sample(c(1, 2), replace = TRUE, size = n, prob = c(1-e, e))
       data<-matrix(0,nrow=n,ncol=3)
       for(t in 1:n){
          if(index[t]==1) 
            data[t,]=mvrnorm(1,mu1,sigma1) 
          else 
            data[t,]=mvrnorm(1,mu2,sigma2)
       }  
       sktests[i] <- Mardiatest(data)
    }
    pwr[j] <- mean(sktests)
}
#plot power vs epsilon
plot(epsilon, pwr, type = "b",xlab = bquote(epsilon), ylim = c(0,1))
abline(h = .1, lty = 3)
se <- sqrt(pwr * (1-pwr) / m) #add standard errors
lines(epsilon, pwr+se, lty = 3)
lines(epsilon, pwr-se, lty = 3)

## -----------------------------------------------------------------------------
library(bootstrap)
set.seed(123)
n<-nrow(scor)
lambda_hat<-eigen(cov(scor))$val
theta_hat<-lambda_hat[order(lambda_hat,decreasing=TRUE)[1]]/sum(lambda_hat)
B<-1e2
thetastar<-numeric(B)
for (b in 1:B) {
  i <- sample(1:n, size = n, replace = TRUE)
  data<-scor[i,]
  lambda<-eigen(cov(data))$val
  thetastar[b]<-lambda[order(lambda,decreasing=TRUE)[1]]/sum(lambda)
}
se<- sd(thetastar)
bias<- mean(thetastar - theta_hat)
round(c(se,bias),4)

## -----------------------------------------------------------------------------
library(bootstrap)
set.seed(123)
n = nrow(scor)
lambda_hat = eigen(cov(scor))$val
theta_hat = lambda_hat[order(lambda_hat,decreasing=TRUE)[1]]/sum(lambda_hat)
theta_j = rep(0,n)
for (i in 1:n) {
    x = scor [-i,]
    lambda = eigen(cov(x))$val
    theta_j[i] = lambda[order(lambda,decreasing=TRUE)[1]]/sum(lambda)
}
#estimated bias of theta_hat
bias_jack = (n-1)*(mean(theta_j)-theta_hat)
#estimated se of theta_hat
se_jack = (n-1)*sqrt(var(theta_j)/n)
print(round(c(bias_jack=bias_jack,se_jack=se_jack),4))

## -----------------------------------------------------------------------------
library(boot)
library(bootstrap)
data(scor, package = "bootstrap")
theta.boot<-function(x,i){
  data<-scor[i,]
  lambda<-eigen(cov(data))$val
  theta.b<-lambda[order(lambda,decreasing=TRUE)[1]]/sum(lambda)
}
boot.obj <- boot(scor, statistic = theta.boot, R = 2000)
alpha <- c(.025, .975)
print(quantile(boot.obj$t, alpha, type=6))

## -----------------------------------------------------------------------------
skewness <- function(x){
  xbar <- mean(x)
  m1 <- mean((x-xbar)^3)
  m2 <- mean((x-xbar)^2)
  return( m1 / m2^1.5 )
}
library(boot)
mu <- 0  # the skewness of normal distribution
n <- 20  # the sample size 
m <- 1e2 # replicate time 

boot.skew <- function(x,i) {
  skewness(x[i])
}

ci.norm <- ci.basic <- ci.perc <- matrix(NA,m,2)
for(i in 1:m){
  X <- rnorm(n)
  de <- boot(data = X, statistic = boot.skew, R = 1e3)
  ci <- boot.ci(de,type = c("norm","basic","perc"))
  ci.norm[i,] <- ci$norm[2:3]
  ci.basic[i,] <- ci$basic[4:5]
  ci.perc[i,] <- ci$percent[4:5]
}

cover.prob <- c(mean(ci.norm[,1]<= mu & ci.norm[,2]>= mu),
                mean(ci.basic[,1]<= mu & ci.basic[,2]>= mu),
                mean(ci.perc[,1]<= mu & ci.perc[,2]>= mu))

left.omit <- c(mean(ci.norm[,1]>= mu),
               mean(ci.basic[,1]>= mu),
               mean(ci.perc[,1]>= mu))

right.omit <- c(mean(ci.norm[,2]<= mu),
               mean(ci.basic[,2]<= mu),
               mean(ci.perc[,2]<= mu))
cover.norm <- matrix(data = c(cover.prob,left.omit,right.omit),nrow = 3,byrow = TRUE,)
rownames(cover.norm) <- c("cover probability","miss on the left","miss on the right")
colnames(cover.norm) <- c("standard normal bootstrap confidence interval","basic bootstrap confidence interval","percentile bootstrap confidence interval")
print(cover.norm)

library(boot)
mu <- sqrt(8/5)  # the skewness of chi-squared distribution
n <- 20  # the sample size in bootstrap
m <- 1e3 # replicate time in Monte Carlo experiments


ci.norm <- ci.basic <- ci.perc <- matrix(NA,m,2)

for(i in 1:m){
  X <- rchisq(n,df = 5)
  de <- boot(data = X, statistic = boot.skew, R = 1e3)
  ci <- boot.ci(de,type = c("norm","basic","perc"))
  ci.norm[i,] <- ci$norm[2:3]
  ci.basic[i,] <- ci$basic[4:5]
  ci.perc[i,] <- ci$percent[4:5]
}


cover.prob <- c(mean(ci.norm[,1]<= mu & ci.norm[,2]>= mu),
                mean(ci.basic[,1]<= mu & ci.basic[,2]>= mu),
                mean(ci.perc[,1]<= mu & ci.perc[,2]>= mu))
left.omit <- c(mean(ci.norm[,1]>= mu),
               mean(ci.basic[,1]>= mu),
               mean(ci.perc[,1]>= mu))
right.omit <- c(mean(ci.norm[,2]<= mu),
               mean(ci.basic[,2]<= mu),
               mean(ci.perc[,2]<= mu))

cover.chisq <- matrix(data = c(cover.prob,left.omit,right.omit),nrow = 3,byrow = TRUE,)
rownames(cover.chisq) <- c("cover probability","miss on the left","miss on the right")
colnames(cover.chisq) <- c("standard normal bootstrap confidence interval","basic bootstrap confidence interval","percentile bootstrap confidence interval")
print(cover.chisq)

## -----------------------------------------------------------------------------
x <- c(60, 68, 68, 77, 74, 77, 80)
y <- c(61, 67, 69, 76, 74, 80, 68)
m<-1e3
u <- c(x,y)
zt <- function(x,y){
  n <- length(x)
  x <- rank(x, ties.method = "average")
  y <- rank(y, ties.method = "average")
  v<-n * ((n+1)/2)^2
  for (i in 1:n) {
    r1 = sum(x[i] * y[i]) - v
    r2 = sum(x[i]^2) -v
    r3 = sum(y[i]^2) - v
  }
  r1/sqrt(r2 * r3)
}
zt_initial <- zt(x,y)
t <- t1 <- t2 <- r_permutation <- numeric(m)
for (i in 1:m) {
  t <- sample(u, size = 14,replace = FALSE)
  t1 <- t[1:7]
  t2 <- t[8:14]
  r_permutation[i] <- zt(t1,t2)
}
p_permutation <- mean(r_permutation < zt_initial)
print(p_permutation)

## -----------------------------------------------------------------------------
Gelman.Rubin<-function(zt){
  zt<-as.matrix(zt)
  n<-ncol(zt)
  k<-nrow(zt)
  zt_means<-rowMeans(zt)
  B<-n*var(zt_means)
  zt_w<-apply(zt,1,"var")
  W<-mean(zt_w)
  v.hat<-W*(n-1)/n+(B/n)
  r.hat<-v.hat/W
  return(r.hat)
}

## -----------------------------------------------------------------------------
a_Euc<-function(a){#定义欧式距离
  norm(a,type="F")
}
k_th<-function(a,k){#计算第k项
  d<-length(a)
  (-1)^k/(factorial(k)*2^k)*a_Euc(a)^(2*k+2)/((2*k+1)*(2*k+2))*exp(lgamma((d+1)/2)+lgamma(k+3/2)-lgamma(k+d/2+1))
}
zt<-function(a){
  s<-0
  k<-0
  g<-1
  while(abs(g)>1e-10){##当第k项绝对值小于1e-10时认为级数收敛
    s<-s+k_th(a,k)
    k<-k+1
    g<-k_th(a,k)
  }
  return(s)
}
a<-matrix(c(1,2),2,1)
zt(a)#计算a=（1，2）时的和

## -----------------------------------------------------------------------------
k<-4
C_k<-function(k,a){
  (a^2*k/(k+1-a^2))^0.5
}
gyj<-function(a,k){
  pt(C_k(k,a),k)
}
zt<-uniroot(function(x) {gyj(x,k)-gyj(x,k-1)},lower=1,upper=2)
zt$root

## -----------------------------------------------------------------------------
y<-c(0.54,0.48,0.33,0.43,0.91,0.21,0.85)
n<-10
m<-3
lambda<-1
zt<-function(y,lambda,x){
  y_sum<-sum(y)
  E<-n*log(x)-x*y_sum-x*m*exp(-lambda)
  return(-E)
}
lambda<-optim(lambda,function(x){zt(y,lambda,x)},method="BFGS")

## -----------------------------------------------------------------------------
rsq <- function(mod) summary(mod)$r.squared
formulas <- list(
  mpg ~ disp,
  mpg ~ I(1 / disp),
  mpg ~ disp + wt,
  mpg ~ I(1 / disp) + wt
)
zt1<-lapply(formulas, lm,data=mtcars)
zt2<-lapply(zt1,rsq)
print(zt2)

## -----------------------------------------------------------------------------
x<-data.frame(cbind(x1=seq(2,8,2),x2=c(1:4)))
zt<-vapply(x, sd, FUN.VALUE = c('a'=0))
print(zt)
#结果显示第一列向量（2，4，6，8）的标准差为2.581989，第二列向量（1，2，3，4）的标准差为1.290994

## -----------------------------------------------------------------------------
x<-data.frame(a=1:5,b=runif(5),c=c(TRUE,FALSE,TRUE,FALSE,TRUE))
zt1<-vapply(x, is.numeric, logical(1))
zt2<-vapply(x[zt1], sd, FUN.VALUE = c('a'=0))
print(zt2)

