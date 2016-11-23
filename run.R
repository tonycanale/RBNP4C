#################################################################
##
## SUPPLEMENTARY MATERIAL FOR  
## Robustifying Bayesian nonparametric mixtures for count data   	        
## by Antonio Canale and Igor Pruenster (c) 2016
##
#################################################################

library(RBNP4C)
?darters
data(darters)
darters$x <- apply(darters[,1:3],1,sum)
darters$sumjx <- darters[,1] + 2*darters[,2] + 3*darters[,3]
n <- nrow(darters)
attach(darters)

#--------------------
# PRIOR CALIBRATION
#--------------------
sigma <- c(0, 0.25, 0.5, 0.75)
# choose the prior expectations for the number of cluster
EK1  <- 10
EK2  <- 22
EK3  <- 30
EK4  <- 40
#find thetas in order to match EK_j
fa <- function(a,n) prod(a+ (1:n) - 1)
theta1 <- theta2 <- theta3 <- theta4 <- rep(NA,length(sigma))
GenStirling <- function(sigma, n, k) 
{
res <- 0
for(j in 0:k)
{
res <- res + ((-1)^(j)) * choose(k,j) *  prod( -j*sigma + 0:(n-1) )
}
res/factorial(k)
}
theta1[1] <- uniroot(function(x) sum(x/(x+(1:n)-1)) - EK1, interval=c(0.01,55))$root
theta2[1] <- uniroot(function(x) sum(x/(x+(1:n)-1)) - EK2, interval=c(0.01,55))$root
theta3[1] <- uniroot(function(x) sum(x/(x+(1:n)-1)) - EK3, interval=c(0.01,55))$root
theta4[1] <- uniroot(function(x) sum(x/(x+(1:n)-1)) - EK4, interval=c(50,155))$root
rfa <- function(theta, sigma, n) prod( 1 + sigma/(theta + 0:(n-1)))
for(i in 2:length(sigma))
{
theta1[i] <- try(uniroot(function(x) x/sigma[i]*(rfa(x,sigma[i],n)-1) - EK1, interval=c(-sigma[i]+0.000001,theta1[i-1]) )$root)
theta2[i] <- try(uniroot(function(x) x/sigma[i]*(rfa(x,sigma[i],n)-1) - EK2, interval=c(-sigma[i]+0.000001,theta2[i-1]) )$root)
theta3[i] <- try(uniroot(function(x) x/sigma[i]*(rfa(x,sigma[i],n)-1) - EK3, interval=c(-sigma[i]+0.000001,theta3[i-1]) )$root)
theta4[i] <- try(uniroot(function(x) x/sigma[i]*(rfa(x,sigma[i],n)-1) - EK4, interval=c(-sigma[i]+0.000001,theta4[i-1]) )$root)
}

#--------------------
# MCMC ITERATIONS
#--------------------
nb <- 30000
nrep <- 40000

#-----------------------------------------------------------------------------------------------------------------------
# FIT Rounded Gaussian

mean_mu <- mean(darters$first) + mean(darters$second[darters$second>0])+mean(darters$third[darters$third>0])

res.rg11<- removalsamplingRG(darters$x, darters$sumjx, darters$J, k=53, nrep=nrep, nb=nb, theta=theta1[1], sigma=sigma[1], algo="polya-urn", mixing_type="2PD", basemeasure_hyperprio=TRUE, atau =2, btau=2, lb = 0, ub = 300, plot.it=TRUE, a_pi=20, b_pi=20,  kap=var(darters$first), mu= mean_mu)
res.rg21<- removalsamplingRG(darters$x, darters$sumjx, darters$J, k=53, nrep=nrep, nb=nb, theta=theta1[2], sigma=sigma[2], algo="polya-urn", mixing_type="2PD", basemeasure_hyperprio=TRUE, atau =2, btau=2, lb = 0, ub = 300, plot.it=TRUE, a_pi=20, b_pi=20, kap=var(darters$first), mu= mean_mu)
res.rg31<- removalsamplingRG(darters$x, darters$sumjx, darters$J, k=53, nrep=nrep, nb=nb, theta=theta1[3], sigma=sigma[3], algo="polya-urn", mixing_type="2PD", basemeasure_hyperprio=TRUE, atau=2, btau=2, lb = 0, ub = 300, plot.it=TRUE, a_pi=20, b_pi=20, kap=var(darters$first), mu= mean_mu)
res.rg41<- removalsamplingRG(darters$x, darters$sumjx, darters$J, k=53, nrep=nrep, nb=nb, theta=theta1[4], sigma=sigma[4], algo="polya-urn", mixing_type="2PD", basemeasure_hyperprio=TRUE, atau = 2, btau=2, lb = 0, ub = 300, plot.it=TRUE, a_pi=20, b_pi=20, kap=var(darters$first), mu= mean_mu)
#--------

res.rg12<- removalsamplingRG(darters$x, darters$sumjx, darters$J, k=53, nrep=nrep, nb=nb, theta=theta2[1], sigma=sigma[1], algo="polya-urn", mixing_type="2PD", basemeasure_hyperprio=TRUE, atau =2, btau=2, lb = 0, ub = 300, plot.it=TRUE, a_pi=20, b_pi=20,  kap=var(darters$first), mu= mean_mu)
res.rg22<- removalsamplingRG(darters$x, darters$sumjx, darters$J, k=53, nrep=nrep, nb=nb, theta=theta2[2], sigma=sigma[2], algo="polya-urn", mixing_type="2PD", basemeasure_hyperprio=TRUE, atau =2, btau=2, lb = 0, ub = 300, plot.it=TRUE, a_pi=20, b_pi=20, kap=var(darters$first), mu= mean_mu)
res.rg32<- removalsamplingRG(darters$x, darters$sumjx, darters$J, k=53, nrep=nrep, nb=nb, theta=theta2[3], sigma=sigma[3], algo="polya-urn", mixing_type="2PD", basemeasure_hyperprio=TRUE, atau=2, btau=2, lb = 0, ub = 300, plot.it=TRUE, a_pi=20, b_pi=20, kap=var(darters$first), mu= mean_mu)
res.rg42<- removalsamplingRG(darters$x, darters$sumjx, darters$J, k=53, nrep=nrep, nb=nb, theta=theta2[4], sigma=sigma[4], algo="polya-urn", mixing_type="2PD", basemeasure_hyperprio=TRUE, atau = 2, btau=2, lb = 0, ub = 300, plot.it=TRUE, a_pi=20, b_pi=20, kap=var(darters$first), mu= mean_mu)
#--------

res.rg13<- removalsamplingRG(darters$x, darters$sumjx, darters$J, k=53, nrep=nrep, nb=nb, theta=theta3[1], sigma=sigma[1], algo="polya-urn", mixing_type="2PD", basemeasure_hyperprio=TRUE, atau =2, btau=2, lb = 0, ub = 300, plot.it=TRUE, a_pi=20, b_pi=20,  kap=var(darters$first), mu= mean_mu)
res.rg23<- removalsamplingRG(darters$x, darters$sumjx, darters$J, k=53, nrep=nrep, nb=nb, theta=theta3[2], sigma=sigma[2], algo="polya-urn", mixing_type="2PD", basemeasure_hyperprio=TRUE, atau =2, btau=2, lb = 0, ub = 300, plot.it=TRUE, a_pi=20, b_pi=20, kap=var(darters$first), mu= mean_mu)
res.rg33<- removalsamplingRG(darters$x, darters$sumjx, darters$J, k=53, nrep=nrep, nb=nb, theta=theta3[3], sigma=sigma[3], algo="polya-urn", mixing_type="2PD", basemeasure_hyperprio=TRUE, atau=2, btau=2, lb = 0, ub = 300, plot.it=TRUE, a_pi=20, b_pi=20, kap=var(darters$first), mu= mean_mu)
res.rg43<- removalsamplingRG(darters$x, darters$sumjx, darters$J, k=53, nrep=nrep, nb=nb, theta=theta3[4], sigma=sigma[4], algo="polya-urn", mixing_type="2PD", basemeasure_hyperprio=TRUE, atau = 2, btau=2, lb = 0, ub = 300, plot.it=TRUE, a_pi=20, b_pi=20, kap=var(darters$first), mu= mean_mu)
#-----

res.rg14<- removalsamplingRG(darters$x, darters$sumjx, darters$J, k=53, nrep=nrep, nb=nb, theta=theta4[1], sigma=sigma[1], algo="polya-urn", mixing_type="2PD", basemeasure_hyperprio=TRUE, atau =1, btau=2, lb = 0, ub = 300, plot.it=TRUE, a_pi=20, b_pi=20,  kap=var(darters$first), mu= mean_mu)
res.rg24<- removalsamplingRG(darters$x, darters$sumjx, darters$J, k=53, nrep=nrep, nb=nb, theta=theta4[2], sigma=sigma[2], algo="polya-urn", mixing_type="2PD", basemeasure_hyperprio=TRUE, atau =1, btau=2, lb = 0, ub = 300, plot.it=TRUE, a_pi=20, b_pi=20, kap=var(darters$first), mu= mean_mu)
res.rg34<- removalsamplingRG(darters$x, darters$sumjx, darters$J, k=53, nrep=nrep, nb=nb, theta=theta4[3], sigma=sigma[3], algo="polya-urn", mixing_type="2PD", basemeasure_hyperprio=TRUE, atau=1, btau=2, lb = 0, ub = 300, plot.it=TRUE, a_pi=20, b_pi=20, kap=var(darters$first), mu= mean_mu)
res.rg44<- removalsamplingRG(darters$x, darters$sumjx, darters$J, k=53, nrep=nrep, nb=nb, theta=theta4[4], sigma=sigma[4], algo="polya-urn", mixing_type="2PD", basemeasure_hyperprio=TRUE, atau = 1, btau=2, lb = 0, ub = 300, plot.it=TRUE, a_pi=20, b_pi=20, kap=var(darters$first), mu= mean_mu)


#--------------------
### FIT Poisson 
#--------------------

res.p11 <- removalsamplingP(darters$x, darters$sumjx, darters$J, k=53, nrep=nrep, nb=nb, theta=theta1[1], sigma=sigma[1], algo="polya-urn", mixing_type="2PD", basemeasure_hyperprio=TRUE, prior="normal",  a=3.6, b=1.78, lb = 0, ub = 300, plot.it=TRUE, a_pi=20, b_pi=20)
res.p21 <- removalsamplingP(darters$x, darters$sumjx, darters$J, k=53, nrep=nrep, nb=nb, theta=theta1[2], sigma=sigma[2], algo="polya-urn", mixing_type="2PD", basemeasure_hyperprio=TRUE, prior="normal",  a=3.6, b=1.78, lb = 0, ub = 300, plot.it=TRUE, a_pi=20, b_pi=20)
res.p31 <- removalsamplingP(darters$x, darters$sumjx, darters$J, k=53, nrep=nrep, nb=nb, theta=theta1[3], sigma=sigma[3], algo="polya-urn", mixing_type="2PD", basemeasure_hyperprio=TRUE, prior="normal",  a=3.6, b=1.78, lb = 0, ub = 300, plot.it=TRUE, a_pi=20, b_pi=20)
res.p41 <- removalsamplingP(darters$x, darters$sumjx, darters$J, k=53, nrep=nrep, nb=nb, theta=theta1[4], sigma=sigma[4], algo="polya-urn", mixing_type="2PD", basemeasure_hyperprio=TRUE, prior="normal",  a=3.6, b=1.78, lb = 0, ub = 300, plot.it=TRUE, a_pi=20, b_pi=20)


res.p12 <- removalsamplingP(darters$x, darters$sumjx, darters$J, k=53, nrep=nrep, nb=nb, theta=theta2[1], sigma=sigma[1], algo="polya-urn", mixing_type="2PD", basemeasure_hyperprio=TRUE, prior="normal",  a=3.6, b=1.78, lb = 0, ub = 300, plot.it=TRUE, a_pi=20, b_pi=20)
res.p22 <- removalsamplingP(darters$x, darters$sumjx, darters$J, k=53, nrep=nrep, nb=nb, theta=theta2[2], sigma=sigma[2], algo="polya-urn", mixing_type="2PD", basemeasure_hyperprio=TRUE, prior="normal",  a=3.6, b=1.78, lb = 0, ub = 300, plot.it=TRUE, a_pi=20, b_pi=20)
res.p32 <- removalsamplingP(darters$x, darters$sumjx, darters$J, k=53, nrep=nrep, nb=nb, theta=theta2[3], sigma=sigma[3], algo="polya-urn", mixing_type="2PD", basemeasure_hyperprio=TRUE, prior="normal",  a=3.6, b=1.78, lb = 0, ub = 300, plot.it=TRUE, a_pi=20, b_pi=20)
res.p42 <- removalsamplingP(darters$x, darters$sumjx, darters$J, k=53, nrep=nrep, nb=nb, theta=theta2[4], sigma=sigma[4], algo="polya-urn", mixing_type="2PD", basemeasure_hyperprio=TRUE, prior="normal",  a=3.6, b=1.78, lb = 0, ub = 300, plot.it=TRUE, a_pi=20, b_pi=20)

res.p13 <- removalsamplingP(darters$x, darters$sumjx, darters$J, k=53, nrep=nrep, nb=nb, theta=theta3[1], sigma=sigma[1], algo="polya-urn", mixing_type="2PD", basemeasure_hyperprio=TRUE, prior="normal",  a=3.6, b=1.78, lb = 0, ub = 300, plot.it=TRUE, a_pi=20, b_pi=20)
res.p23 <- removalsamplingP(darters$x, darters$sumjx, darters$J, k=53, nrep=nrep, nb=nb, theta=theta3[2], sigma=sigma[2], algo="polya-urn", mixing_type="2PD", basemeasure_hyperprio=TRUE, prior="normal",  a=3.6, b=1.78, lb = 0, ub = 300, plot.it=TRUE, a_pi=20, b_pi=20)
res.p33 <- removalsamplingP(darters$x, darters$sumjx, darters$J, k=53, nrep=nrep, nb=nb, theta=theta3[3], sigma=sigma[3], algo="polya-urn", mixing_type="2PD", basemeasure_hyperprio=TRUE, prior="normal",  a=3.6, b=1.78, lb = 0, ub = 300, plot.it=TRUE, a_pi=20, b_pi=20)
res.p43 <- removalsamplingP(darters$x, darters$sumjx, darters$J, k=53, nrep=nrep, nb=nb, theta=theta3[4], sigma=sigma[4], algo="polya-urn", mixing_type="2PD", basemeasure_hyperprio=TRUE, prior="normal",  a=3.6, b=1.78, lb = 0, ub = 300, plot.it=TRUE, a_pi=20, b_pi=20)


res.p14 <- removalsamplingP(darters$x, darters$sumjx, darters$J, k=53, nrep=nrep, nb=nb, theta=theta4[1], sigma=sigma[1], algo="polya-urn", mixing_type="2PD", basemeasure_hyperprio=TRUE, prior="normal",  a=3.6, b=1.78, lb = 0, ub = 300, plot.it=TRUE, a_pi=20, b_pi=20)
res.p24 <- removalsamplingP(darters$x, darters$sumjx, darters$J, k=53, nrep=nrep, nb=nb, theta=theta4[2], sigma=sigma[2], algo="polya-urn", mixing_type="2PD", basemeasure_hyperprio=TRUE, prior="normal",  a=3.6, b=1.78, lb = 0, ub = 300, plot.it=TRUE, a_pi=20, b_pi=20)
res.p34 <- removalsamplingP(darters$x, darters$sumjx, darters$J, k=53, nrep=nrep, nb=nb, theta=theta4[3], sigma=sigma[3], algo="polya-urn", mixing_type="2PD", basemeasure_hyperprio=TRUE, prior="normal",  a=3.6, b=1.78, lb = 0, ub = 300, plot.it=TRUE, a_pi=20, b_pi=20)
res.p44 <- removalsamplingP(darters$x, darters$sumjx, darters$J, k=53, nrep=nrep, nb=nb, theta=theta4[4], sigma=sigma[4], algo="polya-urn", mixing_type="2PD", basemeasure_hyperprio=TRUE, prior="normal",  a=3.6, b=1.78, lb = 0, ub = 300, plot.it=TRUE, a_pi=20, b_pi=20)


#-----------------------------------------------------------------------------------------------------------------------
# NUMBER OF CLUSTERS

g.p.1 <- c(res.p11$clustering$groups, 
           res.p21$clustering$groups,
           res.p31$clustering$groups,
           res.p41$clustering$groups)

g.p.2 <- c(res.p12$clustering$groups, 
           res.p22$clustering$groups,
           res.p32$clustering$groups,
           res.p42$clustering$groups)

g.p.3 <- c(res.p13$clustering$groups, 
           res.p23$clustering$groups,
           res.p33$clustering$groups,
           res.p43$clustering$groups)

g.p.4 <- c(res.p14$clustering$groups, 
           res.p24$clustering$groups,
           res.p34$clustering$groups,
           res.p44$clustering$groups)

#

g.rg.1 <- c(res.rg11$clustering$groups, 
           res.rg21$clustering$groups,
           res.rg31$clustering$groups,
           res.rg41$clustering$groups)

g.rg.2 <- c(res.rg12$clustering$groups, 
           res.rg22$clustering$groups,
           res.rg32$clustering$groups,
           res.rg42$clustering$groups)

g.rg.3 <- c(res.rg13$clustering$groups, 
           res.rg23$clustering$groups,
           res.rg33$clustering$groups,
           res.rg43$clustering$groups)

g.rg.4 <- c(res.rg14$clustering$groups, 
           res.rg24$clustering$groups,
           res.rg34$clustering$groups,
           res.rg44$clustering$groups)


par(mfrow=c(1,2), mar=c(5, 4.5, 4, 1.5) + 0.1)
plot(g.p.1~sigma, ty='b', col=2,  xaxt="n", ylab=expression(E*(K[n]*" | -")), xlim=c(0,.75), ylim=c(10,53),  xlab=expression(sigma),cex.axis=1.5,cex.lab=1.5)
axis(1, at=sigma,cex.axis=1.5)
points(g.p.2~sigma, ty='b', col=2, lty=2) 
points(g.p.3~sigma, ty='b', col=2, lty=3) 
points(g.p.4~sigma, ty='b', col=2, lty=4) 

plot(g.rg.1~sigma, ty='b', col=4, xaxt="n", ylab=expression(E*(K[n]*" | -")), ylim=c(10,53),  xlab=expression(sigma),cex.axis=1.5,cex.lab=1.5)
axis(1, at=sigma,cex.axis=1.5)
points(g.rg.2~sigma, ty='b', col=4, lty=2) 
points(g.rg.3~sigma, ty='b', col=4, lty=3) 
points(g.rg.4~sigma, ty='b', col=4, lty=4) 


