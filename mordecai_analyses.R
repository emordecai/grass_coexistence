# Analyses of the parameterized population growth model

#####################################
### Load packages and model fits
#####################################
require('rjags')
load('samps.Rsave')

#####################################
### Population growth model functions
#####################################
# Effect of competition, as a function of densities x and competition coefficients a
comp = function(x1, x2, x3, a1, a2, a3){
  1/(1 + a1*x1 + a2*x2 + a3*x3)
}

# Population growth for the annual (N[1]) and the perennial seedlings (N[2]) and adults (N[3]), as a function of parameters (pars)
# pars contains the following parameters (subscript 1 is the annual, 2 and 3 are the perennial seedling and adult, respectively):
# g1, g2 = germination fraction
# lambda1, lambda3 = seed production in the absence of competition
# alpha11, alpha13, etc. = competition coefficients
# phi = factor scaling competitive effects of annuals to competitive effects of perennial seedlings
# s2, s3 = over-summer survival

growth = function(N, pars){
  with(pars, {
  M1 = N[1]*g1*lambda1*comp(g1*N[1], g2*N[2], N[3], alpha11, alpha11*phi, alpha13)
  M2 = N[3]*lambda3*comp(g1*N[1], g2*N[2], N[3], alpha31, alpha31*phi, alpha33)
  M3 = N[2]*g2*s2*comp(g1*N[1], g2*N[2], N[3], alpha21, alpha22, alpha23) + N[3]*s3
  c(M1, M2, M3)
  }
  )
}

# Growth rates when rare for the annual (N1) and the perennial (N2)
# First simulates to calculate resident densities, then calculates growth rates when rare as a function of the resident density of the other species
grwr = function(pars){
  # Simulate
  N1 = matrix(NA, nrow=t, ncol=3)
  N1[1,] = c(1, 0, 0)
  N2 = matrix(NA, nrow=t, ncol=3)
  N2[1,] = c(0, 0.5, 0.5)
  colnames(N1) = c("Na", "Np", "A")
  colnames(N2) = c("Na", "Np", "A")
  for (i in 1:(t-1)){
    N1[(i+1),] = growth(N1[i,], pars)
    N2[(i+1),] = growth(N2[i,], pars)
  }

  res.p = tail(N2,1)[2:3]
  res.a = tail(N1,1)[1]
  grwr.a = with(pars, g1*lambda1/(1 + g2*alpha11*phi*res.p[1] + alpha13*res.p[2]))
  grwr.p = with(pars, 0.5*(s3 + (4* lambda3/(1 + alpha31*res.a*g1) * g2*s2/(1 + alpha21*res.a*g1) + s3^2)^0.5))
  c(grwr.a, grwr.p)
}

#####################################
### Analyses
#####################################

# Calculate growth rates when rare across all parameter samples
t = 100
gr = matrix(NA, nrow(fits), 2)
colnames(gr) = c("annual", "perennial")
  
for (j in 1:nrow(fits)){
  pars = as.list(fits[j,])
  gr[j,] = grwr(pars)
}

int = HPDinterval(mcmc(gr))

# Plot the results
plot(gr[,1], gr[,2], pch=16, col=rgb(0,0,0,0.2), xlab="annual growth rate when rare", ylab="perennial growth rate when rare")
abline(h=1)
abline(v=1)
lines(int[1,], rep(mean(gr[,2]), 2), col=2, lwd=2)
lines(rep(mean(gr[,1]), 2), int[2,], col=2, lwd=2)

### Uncertainty analysis
# Uncertainty conbributed by each parameter relative to total model uncertainty
# This looks at uncertainty by grouping the parameters fitted together
# Parameter groupings: g1, g2, s3, annual seed production (lambda1, alpha11, alpha13), perennial seed production (lambda3, alpha31, alpha33), perennial seedling summer survival (s2, alpha21, alpha22, alpha23), perennial seedling competitive effects on annuals and perennial adults (phi)
ind = list(1, 2, 3, c(4:6), c(7:9), c(10:13), 14)
g.a = matrix(NA, nrow(fits), length(ind))
g.p = matrix(NA, nrow(fits), length(ind))
colnames(g.a) = colnames(g.p) = c("g1", "g2", "s3", "ann_seeds", "per_seeds", "seedl_surv", "seedl_effects")

for (k in 1:ncol(g.a)){
  for (j in 1:nrow(fits)){
    pars = as.list(colSums(fits)/nrow(fits))
    pars[ind[[k]]] <-fits[j,ind[[k]]]
    g = grwr(pars)
    g.a[j,k]<-g[1]
    g.p[j,k]<-g[2]
  }
}

h.a = HPDinterval(mcmc(g.a))
h.p = HPDinterval(mcmc(g.p))

# Calculate the percent of the full HPD covered by variation in each parameter
p.a = c()
p.p = c()
for (i in 1:nrow(h.a)){
  p.a[i] <- ((h.a[i,2] - h.a[i,1]))/(int[1,2] - int[1,1])
  p.p[i] <- ((h.p[i,2] - h.p[i,1]))/(int[2,2] - int[2,1])
}

names(p.a) = names(p.p) = colnames(g.a)

# Plot the results
par(mfrow=c(2,1))
barplot(p.a, las=3, names.arg = c("annual\nemergence", "perennial\nemergence", "adult\nsurvival", "annual\nseeds", "perennial\nseeds", "seedling\nsurvival", "seedling\neffects"), ylab = "fraction of total uncertainty", main="annual growth rate when rare")
mtext("A", 3, adj=-0.1, cex=1.3, line=2)
barplot(p.p, las=3, names.arg = c("annual\nemergence", "perennial\nemergence", "adult\nsurvival", "annual\nseeds", "perennial\nseeds", "seedling\nsurvival", "seedling\neffects"), ylab = "fraction of total uncertainty", main="perennial growth rate when rare")
mtext("B", 3, adj=-0.1, cex=1.3, line=2)
par(mfrow=c(1,1))

### Local sensitivity analysis
# Changing each focal parameter by +/- 5% while fixing the other parameters at their posterior means

gr.sens = matrix(NA, ncol(fits), 4)
colnames(gr.sens) = c("plus.a", "plus.p", "minus.a", "minus.p")
rownames(gr.sens) = colnames(fits)

for (j in 1:nrow(gr.sens)){
  pars = colSums(fits)/nrow(fits)
  pars[j]<-pars[j]*1.05
  pars = as.list(pars)
  gr.sens[j,c(1:2)] = grwr(pars)
  pars = colSums(fits)/nrow(fits)
  pars[j]<-pars[j]*0.95
  pars = as.list(pars)
  gr.sens[j,c(3:4)] = grwr(pars)
}
pars = lapply(fits, mean)
gr.base = grwr(pars)
g.sens = cbind((gr.sens[,1]-gr.base[1])/gr.base[1], (gr.sens[,2]-gr.base[2])/gr.base[2], (gr.sens[,3]-gr.base[1])/gr.base[1], (gr.sens[,4]-gr.base[2])/gr.base[2])
colnames(g.sens) = c("plus.a", "plus.p", "minus.a", "minus.p")

# reorder the columns for presentation
g.s = g.sens[c(1,2,3,6,5,4,9,8,7,13,11,12,10,14),]

# Plot the results
names = c(expression("g"[1]), expression("g"[2]), expression("s"[3]), expression(lambda[1]), expression(alpha[11]), expression(alpha[13]), expression(lambda[3]), expression(alpha[31]), expression(alpha[33]), expression("s"[2]), expression(alpha[21]), expression(alpha[22]), expression(alpha[23]), expression(phi))

par(mfrow=c(2,1))
barplot(g.s[,1], las=3, ylab = "effect of 5% increase", main="annual growth rate when rare", names.arg=names)
mtext("A", 3, adj=-0.2, cex=1.3, line=2)
abline(h=0)
barplot(g.s[,2], las=3, ylab = "effect of 5% increase", main="perennial growth rate when rare", names.arg=names)
abline(h=0)
mtext("B", 3, adj=-0.2, cex=1.3, line=2)
par(mfrow=c(1,1))

### Broader sensitivity analysis
# Varying each focal parameter across a feasible range, and calculating GRWR across all samples of the remaining parameters
# This will take a long time!
# Load our simulations to save time
load('sensitivity.Rsave')

# number of parameter values to use
n = 200

# # matrices to store the results
# gr.sens.a = matrix(NA, n, ncol(fits)-2)
# colnames(gr.sens.a) = colnames(fits)[c(1:(ncol(fits)-2))]
# gr.sens.p = matrix(NA, n, ncol(fits)-2)
# colnames(gr.sens.p) = colnames(fits)[c(1:(ncol(fits)-2))]

# gr.sens.a.u = matrix(NA, n, ncol(fits)-2)
# colnames(gr.sens.a.u) = colnames(fits)[c(1:(ncol(fits)-2))]
# gr.sens.p.u = matrix(NA, n, ncol(fits)-2)
# colnames(gr.sens.p.u) = colnames(fits)[c(1:(ncol(fits)-2))]

# gr.sens.a.l = matrix(NA, n, ncol(fits)-2)
# colnames(gr.sens.a.l) = colnames(fits)[c(1:(ncol(fits)-2))]
# gr.sens.p.l = matrix(NA, n, ncol(fits)-2)
# colnames(gr.sens.p.l) = colnames(fits)[c(1:(ncol(fits)-2))]

# range of values to use for parameters that are proportions
yp = seq(0.01, 0.99, length=n)

# range of values to multiply by the parameters that are not proportions
np = seq(0.01, 2, length=n)

# broader range of values to multiply by lambda3
op1 = seq(0.01, 15, length=n)

# broader range of values to multipy by alpha11
op2 = seq(0.01, 5, length=n)

# matrix of parameter values
pmod = cbind(yp, yp, yp, np, op2, np, np, np, op1, np, np, np, yp, yp)

# # list of parameters that are proportions
# prop.list = c(1,2,3,13,14)

# # vectors to store results
# mn.a = c()
# mn.p = c()
# hu.a = c()
# hu.p = c()
# hl.a = c()
# hl.p = c()

# # loop through values for the focal parameter
# for (i in 1:n){
	# # loop through focal parameters
     # for (j in 1:ncol(gr.sens.a)){
     	# # loop through samples of the remaining parameters
          # for (k in 1:nrow(fits)){
               # pars = fits[k,]
               # if (j %in% prop.list) pars[j]<-pmod[i,j] else pars[j]<-pmod[i,j]*pars[j]
               # pars = as.list(pars)
               # tmp = grwr(pars)
               # mn.a[k] = tmp[1]
               # mn.p[k] = tmp[2]
          # }
          # # record the mean and HPD intervals for each value of the focal parameter
          # gr.sens.a[i,j] = mean(mn.a)
          # gr.sens.p[i,j] = mean(mn.p)
          # hpd.a = HPDinterval(mcmc(mn.a))
          # hpd.p = HPDinterval(mcmc(mn.p))
          # gr.sens.a.u[i,j] = hpd.a[2]
          # gr.sens.p.u[i,j] = hpd.p[2]
          # gr.sens.a.l[i,j] = hpd.a[1]
          # gr.sens.p.l[i,j] = hpd.p[1]
     # }
# }

# Calculate the value of the focal parameter (if any) at which the qualitative outcome shifts for the annual and the perennial
tp = matrix(NA, ncol(gr.sens.a), 2)
colnames(tp) = c("annual", "perennial")
rownames(tp) = colnames(gr.sens.a)
for (i in 1:ncol(gr.sens.a)){
     tmpa = which(gr.sens.a[,i]<1)
     if (length(tmpa)>0 & gr.sens.a[1,i]<1) tp[i,1]<-pmod[tail(tmpa,1),i] else {
          if (length(tmpa)>0 & gr.sens.a[1,i]>1) tp[i,1]<-pmod[head(tmpa,1),i]
     }
     tmpp = which(gr.sens.p[,i]>1)
     if (length(tmpp)>0 & gr.sens.p[1,i]>1) tp[i,2]<-pmod[tail(tmpp,1),i] else {
          if (length(tmpp)>0 & gr.sens.p[1,i]<1) tp[i,2]<-pmod[head(tmpp,1),i]
     }
}

# display these values
tp

# Plot the relationships between parameter values and GRWR
# annual
for (i in 1:ncol(gr.sens.a)){
     matplot(pmod[,i], cbind(gr.sens.a[,i], gr.sens.a.u[,i], gr.sens.a.l[,i]), type="l", lty=c(1,2,2), col=1, main=colnames(gr.sens.a)[i], xlab="parameter value", ylab="annual growth rate when rare")
     abline(h=1, col="gray")
}

# perennial
for (i in 1:ncol(gr.sens.p)){
     matplot(pmod[,i], cbind(gr.sens.p[,i], gr.sens.p.u[,i], gr.sens.p.l[,i]), type="l", lty=c(1,2,2), col=1, main=colnames(gr.sens.a)[i], xlab="parameter value", ylab="perennial growth rate when rare")
     abline(h=1, col="gray")
}

