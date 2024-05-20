# Heling Tong 4/10/2023
# functions for the main program

#######################GENERATION###################################

################################################################
# Generation of PK responses    
################################################################

## concentrations obtained from a cohort after receiving a specific dose
outcome.CON <- function(time, cohort.size, fixed.pkpara, var.comp.pkpara, error.var.pkpara, true.dose, h, d, k) {
  
  con.cohort=c() 
  
  for(a in 1:cohort.size) {    
    mean.b <- rep(0,length(fixed.pkpara))
    varcov.b <- diag(var.comp.pkpara, nrow=length(mean.b),ncol=length(mean.b))
    b <- rmvnorm(1, mean.b, varcov.b) # random effects for an individual 
    theta <- fixed.pkpara + b  # PK parameters for an individual
    
    x <- true.dose[h] # dose level
    
    Theta1 <- theta[1,1] # PK parameter
    Theta2 <- theta[1,2] # PK parameter
    
    f <- c() # concentration without residual
    for(g in 1:length(time)) {
      
      f[g] <- model.PK(x,Theta1,Theta2,time[g])
      
    }
    
    ## concentration for an individual
    mean.error <- rep(0,length(time))
    varcov.error <- diag(error.var.pkpara,nrow=length(mean.error),
                         ncol=length(mean.error))
    epsi <- rmvnorm(1, mean.error, varcov.error)
    y <- f+epsi   
    
    s <- 1
    for(u in k:(k+length(time)-1)){
      con.cohort[u] <- y[1,s]
      s <- s+1                              
    }
    k <- k+length(time)
  }
  
  subject <- c(rep(d,length(time)),rep(d+1,length(time)),rep(d+2,length(time)))
  
  dose <- rep(x, cohort.size*length(time))
  t <- rep(time, cohort.size)
  data.cohort.pk <- data.frame(subject=subject, dose=dose, t=t, conc=con.cohort)
  
  return(data.cohort.pk)
}

# data.cohort.pk <- outcome.CON(time, cohort.size, fixed.pkpara, var.comp.pkpara, true.dose, h, d)

################################################################
# Generation of DLT responses    
################################################################

# DLT data
outcome.DLT <- function(cohort.size, h, fixed.DLTpro) {
  
  p.j <- fixed.DLTpro[h] # probability of DLT observe
  DLT <- rbinom(cohort.size, rep(1, cohort.size), p.j)
  
  return(DLT)
  
}

# obsdata_DLT <- outcome.DLT(cohort.size, h, fixed.DLTpro)

################################################################
# Generation of PD responses    
################################################################

# effect data
outcome.EFF <- function(cohort.size, h, fixed.PDmean, var.comp.PDmean) {
  
  mu.cohort <- fixed.PDmean[h]   
  EFF <- rnorm(cohort.size, mu.cohort, var.comp.PDmean)
  
  return(EFF)
}

# obsdata_EFF <- outcome.EFF(cohort.size, h, fixed.PDmean, var.comp.PDmean)


#######################TOXICITY###################################

####################################################################
# dose-exposure model
####################################################################

# population pharmacokinetics (PK) model 
# One-compartment population with bolus input and first-order elimination
# x denotes dose
# v denotes volume
# cl denotes clearance
# t denotes time

model.PK <- function(x, v, cl, t){
  
  (x/v)*exp(-(cl/v)*t)
  
}

# AUC calculation
# t is the upper sampling time

model.AUC <- function(x, v, cl, t){ 
  
  auc <- (x/cl)*(1-exp(-(cl/v)*t))
  
  return(abs(auc))
  
}

# differential w.r.t. v

model.AUC.dv <- function(x, v, cl, t){ 
  
  dv <- ((x*t)/(v^2))*(-exp(-(cl/v)*t))
  
  return(dv)
  
}

# differential w.r.t. cl

model.AUC.dcl <- function(x, v, cl, t){  
  
  dcl <- (x/cl)*(exp(-(cl/v)*t))*((1/cl)+(t/v))-(x/cl^2)
  
  return(dcl)
   
}


# using Bayesian to obtain the estimation of PK parameters 

get.post.PK <- function(obsdata_CON, dose, t){
  
  PKprior <- c(set_prior("normal(0.5, 0.01)", class = "b", nlpar = "Theta1"),
               set_prior("normal(0.06, 0.001)", class = "b",nlpar = "Theta2"),
               set_prior("normal(0.0005, 0.2)", class="sigma"))
  
  PK_bayes_fit_formula <- bf(conc ~ (dose/Theta1)*exp(-(Theta2/Theta1)*t),
                             # Nonlinear variables
                             Theta1+Theta2~1,
                             # Nonlinear fit
                             nl = TRUE)
  
  PK_bayes_fit <- brm(PK_bayes_fit_formula,family=gaussian(), 
                      data = obsdata_CON,prior = PKprior,
                      warmup = 1000, iter = 3000, chains = 2)
  
  fit_res <- summary(PK_bayes_fit)
  
  theta1_mean <- fit_res[["fixed"]][["Estimate"]][1]
  theta2_mean <- fit_res[["fixed"]][["Estimate"]][2]
  theta1_sd <- fit_res[["fixed"]][["Est.Error"]][1]
  theta2_sd <- fit_res[["fixed"]][["Est.Error"]][2]
  sigma.est_mean <- fit_res[["spec_pars"]][["Estimate"]]
  
  # calculate the posterior estimator of parameters 
  theta1.fit <- (PK_bayes_fit[["fit"]]@sim[["samples"]][[1]][["b_Theta1_Intercept"]]+
                   PK_bayes_fit[["fit"]]@sim[["samples"]][[2]][["b_Theta1_Intercept"]])/2
  theta2.fit <- (PK_bayes_fit[["fit"]]@sim[["samples"]][[1]][["b_Theta2_Intercept"]]+
                   PK_bayes_fit[["fit"]]@sim[["samples"]][[2]][["b_Theta2_Intercept"]])/2
  sigma.fit <- (PK_bayes_fit[["fit"]]@sim[["samples"]][[1]][["sigma"]]+
                   PK_bayes_fit[["fit"]]@sim[["samples"]][[2]][["sigma"]])/2
  
  theta1 <- theta1.fit[seq(1, 3000, by=3)]
  theta2 <- theta2.fit[seq(1, 3000, by=3)]
  sigma.est <- sigma.fit[seq(1, 3000, by=3)]
  
  return(list(theta1_mean=theta1_mean, theta2_mean=theta2_mean, 
              theta1_sd=theta1_sd, theta2_sd=theta2_sd,
              sigma.est_mean=sigma.est_mean,
              theta1=theta1, theta2=theta2, sigma.est=sigma.est))
}

# calculate the posterior distribution of AUC 
get.post.AUC <- function(dose, theta1, theta2, t){
  
  AUC_fit <- c()
  
  for (i in 1:1000){
    
    AUC_fit <- rbind(AUC_fit, model.AUC(dose, theta1[i], theta2[i], t))
    
    
  }
  
  # estimation of DOSE-AUC function
  # density_est <- density(AUC_fit)
  # plot(density_est, main = paste("d=", dose[1]))
  
  return(AUC_fit)
  
}


###################################################################
# EXPOSURE-TOXICITY MODEL
# DOSE-DLT RATE MODEL
###################################################################

# logistic regression structural model
get.post.BLRM <- function(obsdata_DLT, auc.lim){
  
    # refDOse=auc_start denotes the reference exposure 
    # initialize the logistic regression structural model 
    model <- LogisticLogNormal(mean = c(-1.808, 0.032),
                               cov = matrix(c(1.17, 0.142, 0.142, 0.365), nrow = 2),
                               refDose = auc.lim)
                   
    # sample from the posterior distribution
    options <- McmcOptions(burnin=100,
                           step=2,
                           samples=2000)
    
    samples <- crmPack::mcmc(data = obsdata_DLT, model = model, options = options)
    
    # plot of probability of DLT (mean curve)
    # plot(samples, model, data) + ggtitle("Posterior")
    # emptydata <- Data(doseGrid = data@doseGrid, placebo = TRUE)
    # priorsamples <- mcmc(emptydata, model, options)
    # plot(priorsamples, model, emptydata) + ggtitle("Prior")
    
    # calculate the posterior probability of DLT 
    alpha1.fit <- samples@data[["alpha0"]]
    alpha2.fit <- samples@data[["alpha1"]]
    
    return(list(alpha1=alpha1.fit, alpha2=alpha2.fit))
    
}

# calculate the posterior distribution of DLT

get_invlogit <- function(x) {
  odds <- exp(x)
  probability <- odds / (1 + odds)
  return(probability)
}

get.post.DLT <-function(alpha1, alpha2, auc, auc.lim){
  
    DLT_rate <- c()
    
    for (i in 1:2000){
      
      DLT_rate <- cbind(DLT_rate, get_invlogit(alpha1[i] + alpha2[i]*log(auc/auc.lim)))
      
    }
    
    # estimation of AUC-DLT function
    # density_est <- density(DLT_rate)
    # plot(density_est, main = paste("AUC=", auc))
    
    # barplot of DLT_rate
    # prob <- table(DLT_rate) / length(DLT_rate)
    # barplot(prob, ylim=c(0,0.01))+title(paste("AUC=", auc))
    
    return(DLT_rate)
    
}


#######################EFFECT##################################### 

###################################################################
# DOSE-EFFECT MODEL
###################################################################

# Emax model using JAGS
model.Emax <- function() {
  
  for (i in 1:n.pat) {
    eff.data[i] ~ dnorm(muI[i],tau)  
    muI[i] <- eta + M*dose.data[i]^gamma/(beta^gamma+dose.data[i]^gamma)
  }
  
  
  eta ~ dgamma(1,10)
  beta ~ dgamma(4,8)
  M ~ dgamma(7.1,17.8)
  tau ~ dgamma(.1,.1)
  gamma ~ dgamma(1/9,1/18)
  
}

# get posterior distributions for the Emax model
get.post.Emax <- function(obsdata_EFF, dose, n.pat) {
  
  dat.list <- list(eff.data=obsdata_EFF, dose.data=dose, n.pat=n.pat)
  params <- c("eta","beta","M","gamma","tau")
  
  jags.inits <- function(){
    list("eta"=1,"beta"=0.5,"M"=10,"gamma"=2,"tau"=1/4)
  }
  
  mcmcmodel <- jags(data=dat.list, model.file=model.Emax,
                    inits=jags.inits,DIC=F,progress.bar="none",
                    param=c("eta","beta","M","gamma","tau"),
                    n.chains=1,n.iter=3000,n.burnin=1000,n.thin=1)
  
  param.mcmc <- as.mcmc(mcmcmodel)[[1]]
  param.mcmc <- cbind(param.mcmc[,c(2,1,4,3)],sqrt(1/param.mcmc[,5])) 
  
  return(param.mcmc)  
}


#######################SELECT##################################### 

###################################################################
# DOSE-SELECTION FUNCTION
###################################################################
 
# calculate the dose-selection function

dose.select <- function(pi.hat, mu.hat, W){
  
  dis <- W*mu.hat + (1-W)*pi.hat
  
  return(dis)
  
}

