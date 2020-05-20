# Fitting California grassland demographic parameters using Bayesian models
### All priors informed by our knowledge of the system ####
# This is the main model presented in the text

### Required bug files:
# bayes_survival.bug, bayes_seed_production.bug, bayes_seedling_survival.bug

#################################
### Loading data and packages
#################################
require('rjags')
load('mordecai_data_sources.Rsave')

#############################
### Bayesian model settings
#############################

n.chains<-2
thin = 5*seq(1,1000)    # thinning the samples, taking 1 in 5

###################################
### Fit Bayesian models
###################################

###################################
### Adult perennial survival (s3)

## Fit a beta-binomial model
p.s3.g <- jags.model('bayes_survival.bug', data=list('Y'=po_s3$Y, 'n'=po_s3$n, 'N'=nrow(po_s3), 'paa'=8, 'pab'=1, 'pba' = 0.01, 'pbb' = 1), n.chains=n.chains, n.adapt=5000)
pom.s3.g<-coda.samples(p.s3.g,c('da', 'db'),5000)
plot(pom.s3.g)
s3.g = rbind(pom.s3.g[[1]][thin,],pom.s3.g[[2]][thin,])
bfun = function(df) (df[,1])/(df[,1] + df[,2])

# Plot the fitted model against the data
plot(density(bfun(s3.g)), xlim=range(po_s3$Y/po_s3$n), main="")
points(po_s3$Y/po_s3$n, rep(0,nrow(po_s3)), pch=16)
s3 = bfun(s3.g)

####################################
### Spring establishment (g1 and g2)

## Annuals
### Fitting a beta-binomial model
p.g1.g<-jags.model('bayes_survival.bug', data=list('Y'=post.g1$Y, 'n'=post.g1$n, 'N'=nrow(post.g1), 'paa'=1, 'pab'=1, 'pba' = 1, 'pbb' = 1), n.chains=n.chains, n.adapt=5000)
pom.g1.g<-coda.samples(p.g1.g,c('da', 'db'),5000)
plot(pom.g1.g)
g1.g = rbind(pom.g1.g[[1]][thin,],pom.g1.g[[2]][thin,])

# plotting the fitted model against the data
plot(density(bfun(g1.g)))
points(post.g1$Y/post.g1$n, c(0,0,0.05), pch=16, main="")

# save the parameter samples
g1 = bfun(g1.g)


## Perennial models
### Fitting a beta-binomial model
p.g2.g<-jags.model('bayes_survival.bug', data=list('Y'=post.g2$Y, 'n'=post.g2$n, 'N'=nrow(post.g2), 'paa'=1, 'pab'=1, 'pba' = 1, 'pbb' = 1), n.chains=n.chains, n.adapt=5000)
pom.g2.g<-coda.samples(p.g2.g,c('da', 'db'),5000)
plot(pom.g2.g)
g2.g = rbind(pom.g2.g[[1]][thin,],pom.g2.g[[2]][thin,])

# plotting the fitted model against the data
plot(density(bfun(g2.g)), main="")
points(post.g2$Y/post.g2$n, rep(0,nrow(post.g2)), pch=16)

# save the parameter samples
g2 = bfun(g2.g)

#########################
### Seed production

### Annual
# model fitting
p.com<-jags.model('bayes_seed_production.bug', data=list('Y'=log(po_com$Y),'db'=po_com$db, 'dA'=po_com$dA, 'N'=nrow(po_com),'a1' = 1, 'a2' = 9, 'c1' = 1, 'c2' = 9,'la' = 30,'lb' = 1, 't1' = 0.001, 't2' = 0.001), n.chains=n.chains,  n.adapt=5000)
pom.com<-coda.samples(p.com,c('lambda','ab','aA', 'sigma'),5000)
plot(pom.com)
seeds.com = rbind(pom.com[[1]][thin,-4],pom.com[[2]][thin,-4])
colnames(seeds.com) = c('alpha13', 'alpha11', 'lambda1')
seeds.com = data.frame(seeds.com)

# Plot the data and a sample of the model fits
par(mfrow=c(1,1))
plot(Y~db, data=po_com, col=2, pch=NA, xlab="density of annuals", ylab="per-capita seed production of annuals")
x=seq(1,max(po_com$db)+100)
n = sample(c(1:nrow(seeds.com)), size=100, replace=F)
for (i in 1:length(n)){
  y = with(seeds.com, lambda1[n[i]]/(1 + x*alpha11[n[i]] + mean(po_com$dA)*alpha13[n[i]]))
  lines(x,y, col="gray50")
}
y = with(seeds.com, mean(lambda1)/(1 + x*mean(alpha11) + mean(po_com$dA)*mean(alpha13)))
lines(x,y, lwd=3)
points(Y~db, data=po_com, col=2, pch=16)

plot(Y~dA, data=po_com, col=2, pch=NA, xlab="density of adult perennials", ylab="per-capita seed production of annuals")
x=seq(0,max(po_com$dA)+10)
for (i in 1:length(n)){
  y = with(seeds.com, lambda1[n[i]]/(1 + mean(po_com$db)*alpha11[n[i]] + x*alpha13[n[i]]))
  lines(x,y, col="gray50")
}
y = with(seeds.com, mean(lambda1)/(1 + mean(po_com$db)*mean(alpha11) + x*mean(alpha13)))
lines(x,y, lwd=3)
points(Y~dA, data=po_com, col=2, pch=16)

### Perennial
# model fitting
p.per<-jags.model('bayes_seed_production.bug', data=list('Y'=log(po_per$Y),'db'=po_per$db, 'dA'=po_per$dA, 'N'=nrow(po_per),'a1' = 1, 'a2' = 9, 'c1' = 1, 'c2' = 9, 'la' = 200,'lb' = 1, 't1' = 0.001, 't2' = 0.001), n.chains=n.chains,  n.adapt=5000)
pom.per<-coda.samples(p.per,c('lambda','ab','aA', 'sigma'),5000)
plot(pom.per)
seeds.per = rbind(pom.per[[1]][thin,-4],pom.per[[2]][thin,-4])
colnames(seeds.per) = c('alpha33', 'alpha31', 'lambda3')
seeds.per = data.frame(seeds.per)

# Plot the data against some model fits
par(mfrow=c(1,1))
plot(Y~db, data=po_per, col=2, pch=NA, xlab="density of annuals", ylab="per-capita seed production of perennials")
n = sample(c(1:nrow(seeds.per)), size=50, replace=F)
x=seq(1,max(po_per$db)+1)
for (i in 1:length(n)){
  y = with(seeds.per, lambda3[n[i]]/(1 + x*alpha31[n[i]] + mean(po_per$dA)*alpha33[n[i]]))
  lines(x,y, col="gray50")
}
y = with(seeds.per, mean(lambda3)/(1 + x*mean(alpha31) + mean(po_per$dA)*mean(alpha33)))
lines(x,y, lwd=3)
points(Y~db, data=po_per, col=2, pch=16)

plot(Y~dA, data=po_per, col=2, pch=NA, xlab="density of perennials", ylab="per-capita seed production of perennials")
x=seq(1,max(po_per$dA)+1)
for (i in 1:length(n)){
  y = with(seeds.per, lambda3[n[i]]/(1 + mean(po_per$db)*alpha31[n[i]] + x*alpha33[n[i]]))
  lines(x,y, col="gray50")
}
y = with(seeds.per, mean(lambda3)/(1 + mean(po_per$db)*mean(alpha31) + x*mean(alpha33)))
lines(x,y, lwd=3)
points(Y~dA, data=po_per, col=2, pch=16)

#########################################
## Over-summer survival of perennial seedlings
# model fitting
p.sl<-jags.model('bayes_seedling_survival.bug', data=list('Y'=po_sl$Y,'db'=po_sl$db,'de'=po_sl$de,'dA'=po_sl$dA,'n'=po_sl$n, 'N'=nrow(po_sl),'a1' = 1, 'a2' = 9, 'b1' = 1, 'b2' = 9, 'c1' = 1, 'c2' = 9, 'pa' = 1,'pb' = 1, 'ma' = 1, 'mb' = 9), n.chains=n.chains,  n.adapt=5000)
pom.sl<-coda.samples(p.sl,c('m','ab', 'ae', 'aA'),5000)
plot(pom.sl)
seeds.sl = rbind(pom.sl[[1]][thin,],pom.sl[[2]][thin,])
colnames(seeds.sl) = c('alpha23', 'alpha21', 'alpha22', 's2')

# Plot the data against some model fits
par(mfrow=c(1,1))
plot(Y/n~db, data=po_sl, col=2, pch=NA, xlab="density of annuals", ylab="survival probability of perennial seedlings")
n = sample(c(1:nrow(seeds.sl)), size=50, replace=F)
x=seq(0,max(po_sl$db)+10)
for (i in 1:length(n)){
     y=seeds.sl[n[i],4]/(1 + x*seeds.sl[n[i],2] + mean(po_sl$dA)*seeds.sl[n[i],1] + mean(po_sl$de)*seeds.sl[n[i],3])
     lines(x,y, col="gray50")
}
y=mean(seeds.sl[,4])/(1 + x*mean(seeds.sl[,2]) + mean(po_sl$dA)*mean(seeds.sl[,1]) + mean(po_sl$de)*mean(seeds.sl[,3]))
lines(x,y, lwd=3)
points(Y/n~db, data=po_sl, col=2, pch=16)

plot(Y/n~de, data=po_sl, col=2, pch=NA, xlab="density of perennial seedlings", ylab="survival probability of perennial seedlings")
n = sample(c(1:nrow(seeds.sl)), size=50, replace=F)
x=seq(0,max(po_sl$de)+10)
for (i in 1:length(n)){
     y=seeds.sl[n[i],4]/(1 + mean(po_sl$db)*seeds.sl[n[i],2] + mean(po_sl$dA)*seeds.sl[n[i],1] + x*seeds.sl[n[i],3])
     lines(x,y, col="gray50")
}
y=mean(seeds.sl[,4])/(1 + mean(po_sl$db)*mean(seeds.sl[,2]) + mean(po_sl$dA)*mean(seeds.sl[,1]) + x*mean(seeds.sl[,3]))
lines(x,y, lwd=3)
points(Y/n~de, data=po_sl, col=2, pch=16)

plot(Y/n~dA, data=po_sl, col=2, pch=NA, xlab="density of perennial adults", ylab="survival probability of perennial seedlings")
n = sample(c(1:nrow(seeds.sl)), size=50, replace=F)
x=seq(0,max(po_sl$dA)+10)
for (i in 1:length(n)){
     y=seeds.sl[n[i],4]/(1 + mean(po_sl$db)*seeds.sl[n[i],2] + x*seeds.sl[n[i],1] + mean(po_sl$de)*seeds.sl[n[i],3])
     lines(x,y, col="gray50")
}
y=mean(seeds.sl[,4])/(1 + mean(po_sl$db)*mean(seeds.sl[,2]) + x*mean(seeds.sl[,1]) + mean(po_sl$de)*mean(seeds.sl[,3]))
lines(x,y, lwd=3)
points(Y/n~dA, data=po_sl, col=2, pch=16)

## Perennial seedling competitive effects
# Since we have no information on alpha12 and alpha32, we assume that they are no greater than alpha11 and alpha31, respectively
# Set alpha12 = phi*alpha11 and alpha32 = phi*alpha31 where phi ranges from zero to one
phi = rbeta(2*length(thin), 1, 1)
alpha12 = phi*seeds.com$alpha11
alpha32 = phi*seeds.per$alpha31


####################################
### Save posterior MCMC samples
####################################
fits = data.frame(g1, g2, s3, seeds.com, seeds.per, seeds.sl, phi, alpha12, alpha32)
save(fits, file="samps.Rsave")

