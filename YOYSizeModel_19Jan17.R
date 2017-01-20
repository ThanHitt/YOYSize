setwd("C:/Users/nhitt/Documents/PROJECTS/REGBKT/YOYsize/jagsMod")
require(jagsUI)
require(doBy)
require(car)
memory.limit(size=32692)

########################################################
# Load fish data
d <- read.csv("REGBKT_MASTER_16Dec16.csv", header=T, sep=",") # Individual body size data
d$Date <- as.POSIXct(d$Date, format="%m/%d/%Y")
d$doy <- as.numeric(strftime(d$Date, format = "%j"))
d <- subset(d, Species=="BKT")
d <- subset(d, doy < 305) # Remove samples collected after October
d <- subset(d, doy > 120) # Remove samples collected before May

######################################################
# Define adults and YOYs based on length and season
# Based on length-frequency histograms in Petty et al. (2005): "adults" defined as > 72 mm TL (65 mm SL) if sample collected before June 1 and as > 100 mm TL if sampled on or after June 1
# The purpose of this is to exclude small age-1 fish from counts of YOY when samples are collected in early spring
d$YOY <- NA
d$Adult <- NA
d$YOY[d$TL_mm <= 100] <- 1
d$YOY[d$TL_mm > 72 & d$doy < 152] <- 0
d$YOY[is.na(d$YOY)] <- 0
d$Adult[d$YOY==0] <- 1
d$Adult[d$YOY==1] <- 0

# Add YOY abundance
yoy <- summaryBy(YOY ~ VisitID, data=d, FUN=sum) # YOY abundance
yoy$YOY.sum <- (yoy$YOY.sum - mean(yoy$YOY.sum))/sd(yoy$YOY.sum) # z-scores by grand mean and sd
colnames(yoy) <- c("VisitID","yoyN")
d <- merge(d, yoy, by="VisitID", all=F) # limit dataset to samples where YOY were observed

# Add adult abundance
adu <- summaryBy(Adult ~ VisitID, data=d, FUN=sum) # adult abundance
d1 <- merge(d, adu, by="VisitID", all=T)
d1$Adult.sum[is.na(d1$Adult.sum)] <- 0 # indicate where no adults were observed
adu1 <- summaryBy(Adult.sum ~ VisitID, data=d1, FUN=mean)
colnames(adu1) <- c("VisitID","aduN")
adu1$aduN <- (adu1$aduN - mean(adu1$aduN))/sd(adu1$aduN) # z-scores by grand mean and sd
d <- merge(d, adu1, by="VisitID", all=T)

# Summarize YOY dataset
dy <- droplevels(subset(d, YOY==1))
nInds <- nrow(dy)
nSites <- length(unique(dy$SiteID))
nYears <- length(unique(dy$Year))
Years <- as.data.frame(seq(min(dy$Year),max(dy$Year),1))
colnames(Years) <- c("Year")

# Load covariates
sa <- read.csv("ClimateCovs.csv", header=T, sep=",") # Sample-level covariates
sa <- droplevels(subset(sa, SiteID %in% unique(dy$SiteID)))
si <- read.csv("SiteCovs.csv", header=T, sep=",") # Site-level covariates
si <- droplevels(subset(si, SiteID %in% unique(dy$SiteID)))
AllVisits <- as.data.frame(unique(sa$VisitID))
colnames(AllVisits) <- c("VisitID")
d <- subset(d, SiteID %in% unique(dy$SiteID))
d <- subset(d, PassNo==1)

# Index YOY BKT abundance in "yoyN1" by j=sites and m=years [j,m]
yN <- summaryBy(yoyN ~ SiteID + Year + VisitID, data=d, FUN=mean)
yNx <- merge(yN, AllVisits, by="VisitID", all=T)
yNx$yoyN.mean[is.na(yNx$yoyN.mean)] <- 0 # Assume mean-effects for unsampled sites/years
yoyN1 <- matrix(yNx$yoyN.mean, nrow=nSites, ncol=nYears, byrow=T)
rm(yN, yNx)

# Index adult BKT abundance in "aduN1" by j=sites and m=years [j,m]
aN <- summaryBy(aduN ~ SiteID + Year + VisitID, data=d, FUN=mean)
aNx <- merge(aN, AllVisits, by="VisitID", all=T)
aNx$aduN.mean[is.na(aNx$aduN.mean)] <- 0 # Assume mean-effects for unsampled sites/years
aduN1 <- matrix(aNx$aduN.mean, nrow=nSites, ncol=nYears, byrow=T)
rm(aN, aNx)

# Index day-of-year in "dy" by individuals [i]
doyx <- summaryBy(doy ~ VisitID, data=dy, FUN=mean)
doyx$doy.mean <- (doyx$doy.mean - mean(doyx$doy.mean))/sd(doyx$doy.mean) # z-scores by grand mean and sd
colnames(doyx) <- c("VisitID","doy")
dy <- merge(dy, doyx, by="VisitID", all=F)

# Index day-of-year in "doy1" by j=sites and m=years [j,m]
doyx <- summaryBy(doy.y ~ SiteID + Year + VisitID, data=dy, FUN=mean)
doyx$doy.y.mean <- (doyx$doy.y.mean-mean(doyx$doy.y.mean))/sd(doyx$doy.y.mean) # z-scores from grand mean and sd
doyxx <- merge(doyx, AllVisits, by="VisitID", all=T)
doyxx$doy.y.mean[is.na(doyxx$doy.y.mean)] <- 0 # Assume mean-effects for unsampled sites/years
doy1 <- matrix(doyxx$doy.y.mean, nrow=nSites, ncol=nYears, byrow=T)
rm(doyx, doyxx)

# Index fall precip in "dy" by individuals [i]
fpx <- summaryBy(Fall.ppt ~ VisitID, data=sa, FUN=mean)
fpx$Fall.ppt.mean <- (fpx$Fall.ppt.mean - mean(fpx$Fall.ppt.mean))/sd(fpx$Fall.ppt.mean) # z-scores by grand mean and sd
colnames(fpx) <- c("VisitID","fp")
dy <- merge(dy, fpx, by="VisitID", all=F)

# Index fall precip in "fp1" by j=sites and m=years [j,m]
fpx <- summaryBy(Fall.ppt ~ SiteID + Year, data=sa, FUN=mean)
fpx$zscore <- (fpx$Fall.ppt.mean-mean(fpx$Fall.ppt.mean))/sd(fpx$Fall.ppt.mean) 
fp1 <- matrix(fpx$zscore, nrow=nSites, ncol=nYears, byrow=T)
rm(fpx)

# Index fall temp in "dy" by individuals [i]
ftx <- summaryBy(Fall.tmean ~ VisitID, data=sa, FUN=mean)
ftx$Fall.tmean.mean <- (ftx$Fall.tmean.mean - mean(ftx$Fall.tmean.mean))/sd(ftx$Fall.tmean.mean) # z-scores by grand mean and sd
colnames(ftx) <- c("VisitID","ft")
dy <- merge(dy, ftx, by="VisitID", all=F)

# Index fall temp in "ft1" by j=sites and m=years [j,m]
ftx <- summaryBy(Fall.tmean ~ SiteID + Year, data=sa, FUN=mean)
ftx$zscore <- (ftx$Fall.tmean.mean-mean(ftx$Fall.tmean.mean))/sd(ftx$Fall.tmean.mean) # z-scores from grand mean and sd
ft1 <- matrix(ftx$zscore, nrow=nSites, ncol=nYears, byrow=T)
rm(ftx)

# Index winter precip in "dy" by individuals [i]
wpx <- summaryBy(Winter.ppt ~ VisitID, data=sa, FUN=mean)
wpx$Winter.ppt.mean <- (wpx$Winter.ppt.mean - mean(wpx$Winter.ppt.mean))/sd(wpx$Winter.ppt.mean) # z-scores by grand mean and sd
colnames(wpx) <- c("VisitID","wp")
dy <- merge(dy, wpx, by="VisitID", all=F)

# Index winter precip in "wp1" by j=sites and m=years [j,m]
wp1 <- matrix(nrow=nSites, ncol=nYears)
wpx <- summaryBy(Winter.ppt ~ SiteID + Year, data=sa, FUN=mean)
wpx$zscore <- (wpx$Winter.ppt.mean-mean(wpx$Winter.ppt.mean))/sd(wpx$Winter.ppt.mean) # z-scores from grand mean and sd
wp1 <- matrix(wpx$zscore, nrow=nSites, ncol=nYears, byrow=T)
rm(wpx)

# Index winter temp in "dy" by individuals [i]
wtx <- summaryBy(Winter.tmean ~ VisitID, data=sa, FUN=mean)
wtx$Winter.tmean.mean <- (wtx$Winter.tmean.mean - mean(wtx$Winter.tmean.mean))/sd(wtx$Winter.tmean.mean) # z-scores by grand mean and sd
colnames(wtx) <- c("VisitID","wt")
dy <- merge(dy, wtx, by="VisitID", all=F)

# Index winter temp in "wt1" by j=sites and m=years [j,m]
wtx <- summaryBy(Winter.tmean ~ SiteID + Year, data=sa, FUN=mean)
wtx$zscore <- (wtx$Winter.tmean.mean-mean(wtx$Winter.tmean.mean))/sd(wtx$Winter.tmean.mean) # z-scores from grand mean and sd
wt1 <- matrix(wtx$zscore, nrow=nSites, ncol=nYears, byrow=T)
rm(wtx)

# Index spring precip in "dy" by individuals [i]
spx <- summaryBy(Spring.ppt ~ VisitID, data=sa, FUN=mean)
spx$Spring.ppt.mean <- (spx$Spring.ppt.mean - mean(spx$Spring.ppt.mean))/sd(spx$Spring.ppt.mean) # z-scores by grand mean and sd
colnames(spx) <- c("VisitID","sp")
dy <- merge(dy, spx, by="VisitID", all=F)

# Index spring precip in "sp1" by j=sites and m=years [j,m]
spx <- summaryBy(Spring.ppt ~ SiteID + Year, data=sa, FUN=mean)
spx$zscore <- (spx$Spring.ppt.mean-mean(spx$Spring.ppt.mean))/sd(spx$Spring.ppt.mean) # z-scores from grand mean and sd
sp1 <- matrix(spx$zscore, nrow=nSites, ncol=nYears, byrow=T)
rm(spx)

# Index spring temp in "dy" by individuals [i]
stx <- summaryBy(Spring.tmean ~ VisitID, data=sa, FUN=mean)
stx$Spring.tmean.mean <- (stx$Spring.tmean.mean - mean(stx$Spring.tmean.mean))/sd(stx$Spring.tmean.mean) # z-scores by grand mean and sd
colnames(stx) <- c("VisitID","st")
dy <- merge(dy, stx, by="VisitID", all=F)

# Index spring temp in "st1" by j=sites and m=years [j,m]
stx <- summaryBy(Spring.tmean ~ SiteID + Year, data=sa, FUN=mean)
stx$zscore <- (stx$Spring.tmean.mean-mean(stx$Spring.tmean.mean))/sd(stx$Spring.tmean.mean) # z-scores from grand mean and sd
st1 <- matrix(stx$zscore, nrow=nSites, ncol=nYears, byrow=T)
rm(stx)

########################################################
# Write model
sink("model.txt")
cat("
    model {
    for(i in 1:n){
      y[i] ~ dnorm(y.hat[i], tau.y[sid[i],yid[i]])
      y.hat[i] <- alpha + site.ran.mu[sid[i]] + b1*doy[i] + b2*yoyN[i] + b3*aduN[i] + b4*fp[i] + 
      b5*ft[i] + b6*wp[i] + b7*wt[i] + b8*sp[i] + b9*st[i]
    
      res[i] <- y[i] - y.hat[i]
      pred[i] <- y.hat[i]    
      sq[i] <- pow(res[i], 2)
      y.new[i] ~ dnorm(y.hat[i], tau.y[sid[i],yid[i]])
      sq.new[i] <- pow(y.new[i] - pred[i], 2)
    }
    
    for(j in 1:J){
      # Random site effect on mean [j=sites]
      site.ran.mu[j] ~ dnorm(site.hat[j], tau.site.mu)
      site.hat[j] <- b10*elev[j] + b11*uba[j] + b12*grad[j]

      # Random site effect on SD [j=sites]
      site.ran.sigma[j] ~ dnorm(site.ran.sig.mu[j], tau.site.sigma)
      site.ran.sig.mu[j] <- b22*elev[j] + b23*uba[j] + b24*grad[j]

      # Random site SD [j=sites, m=years]
      for(m in 1:M){
        tau.y[j,m] <- pow(sigma[j,m], -2)
        log(sigma[j,m]) <- log.sigma[j,m]
        log.sigma[j,m] ~ dnorm(mu.sigma[j,m], inv.omega.sigma.squared)
        mu.sigma[j,m] <- Mu.Sigma + site.ran.sigma[j] + b13*doy1[j,m] + b14*yoyN1[j,m] + b15*aduN1[j,m] + b16*fp1[j,m] + 
                         b17*ft1[j,m] + b18*wp1[j,m] + b19*wt1[j,m] + b20*sp1[j,m] + b21*st1[j,m]
      }
    }

    # Priors
    b1 ~ dnorm(0, 0.001)
    b2 ~ dnorm(0, 0.001)
    b3 ~ dnorm(0, 0.001)
    b4 ~ dnorm(0, 0.001)
    b5 ~ dnorm(0, 0.001)
    b6 ~ dnorm(0, 0.001)
    b7 ~ dnorm(0, 0.001)
    b8 ~ dnorm(0, 0.001)
    b9 ~ dnorm(0, 0.001)
    b10 ~ dnorm(0, 0.001)
    b11 ~ dnorm(0, 0.001)
    b12 ~ dnorm(0, 0.001)
    b13 ~ dnorm(0, 0.001)
    b14 ~ dnorm(0, 0.001)
    b15 ~ dnorm(0, 0.001)
    b16 ~ dnorm(0, 0.001)
    b17 ~ dnorm(0, 0.001)
    b18 ~ dnorm(0, 0.001)
    b19 ~ dnorm(0, 0.001)
    b20 ~ dnorm(0, 0.001)
    b21 ~ dnorm(0, 0.001)
    b22 ~ dnorm(0, 0.001)
    b23 ~ dnorm(0, 0.001)
    b24 ~ dnorm(0, 0.001)
    alpha ~ dnorm(0, 0.001)
    tau.site.mu <- pow(sig.site.mu, -2)
    sig.site.mu ~ dunif(0, 20)
    Mu.Sigma ~ dnorm(0, 0.001)
    Med.Sigma <- exp(Mu.Sigma)
    inv.omega.sigma.squared <- pow(omega.sigma, -2)
    omega.sigma ~ dunif(0, 20)
    tau.site.sigma <- pow(sig.site.sigma, -2)
    sig.site.sigma ~ dunif(0, 20)
    
    # Derived parameters
    fit <- sum(sq[])
    fit.new <- sum(sq.new[])
    test <- step(fit.new - fit)
    bpvalue <- mean(test)
    } # end model
    ",fill = TRUE)
sink()

# Load data
dat <- list(y = dy$TL_mm,
            sid = as.numeric(as.factor(as.numeric(dy$SiteID))),
            yid = as.numeric(as.factor(dy$Year)),
            n = nInds,
            J = nSites,
            M = nYears,
            elev = as.numeric(scale(si$Elev_m)),
            uba = as.numeric(scale(si$UBA_ha)),
            grad = as.numeric(scale(asin(sqrt(si$Gradient_per/100)))),
            yoyN = dy$yoyN,
            yoyN1 = yoyN1,
            aduN = dy$aduN,
            aduN1 = aduN1,
            doy = dy$doy.y,
            doy1 = doy1,
            fp = dy$fp,
            fp1 = fp1,
            ft = dy$ft,
            ft1 = ft1,
            wp = dy$wp,
            wp1 = wp1,
            wt = dy$wt,
            wt1 = wt1,
            sp = dy$sp,
            sp1 = sp1,
            st = dy$st,
            st1 = st1)

# Initial values
inits <- function(){
  list(alpha = rnorm(1),
       Mu.Sigma = rnorm(1), 
       omega.sigma = runif(1),
       sig.site.mu = runif(1),
       log.sigma = matrix(runif(nSites*nYears), nrow=nSites, ncol=nYears),
       b1 = rnorm(1),
       b2 = rnorm(1),
       b3 = rnorm(1),
       b4 = rnorm(1),
       b5 = rnorm(1),
       b6 = rnorm(1),
       b7 = rnorm(1),
       b8 = rnorm(1),
       b9 = rnorm(1),
       b10 = rnorm(1),
       b11 = rnorm(1),
       b12 = rnorm(1),
       b13 = rnorm(1),
       b14 = rnorm(1),
       b15 = rnorm(1),
       b16 = rnorm(1),
       b17 = rnorm(1),
       b18 = rnorm(1),
       b19 = rnorm(1),
       b20 = rnorm(1),
       b21 = rnorm(1),
       b22 = rnorm(1),
       b23 = rnorm(1),
       b24 = rnorm(1))
}

# Set parameters to monitor
params <- c("b1","b2","b3","b4","b5","b6","b7","b8","b9","b10","b11","b12","b13","b14","b15","b16",
"b17","b18","b19","b20","b21","b22","b23","b24","alpha","sigma","site.ran.sigma","site.ran.mu","fit",
"fit.new","bpvalue")

# MCMC settings
ni <- 100000
nt <- 100
nb <- 1000
nc <- 3

########################################################
# Call JAGS from R
out <- jags(dat, inits, params, "model.txt", n.chains=nc, n.thin=nt, n.iter=ni, n.burnin=nb)
save(out, file="out_19Jan17_v3.R")

########################################################
pp.check(out, actual="fit", new="fit.new")
#plot(out$sims.list$fit, out$sims.list$fit.new, pch=16)
require(MASS)
chisq.test(out$sims.list$fit, out$sims.list$fit.new)
