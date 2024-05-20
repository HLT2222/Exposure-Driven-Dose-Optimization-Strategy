# This is a script for simulation study
# rm(list = ls())

full_model <- function (PFIM.directory,
                        Func.directory,
                        true.dose,
                        cohort.size,
                        cohort.num,
                        W,
                        stop.freq = 6,
                        dose.skip = 1,
                        up.down.num = 4,
                        target.toxicity = 0.25,
                        target.efficacy = 0.1,
                        tox_cut_off = 0.90,
                        eff_cut_off = 0.90,
                        lower.sampling = 0,
                        upper.sampling = 20,
                        fixed.pkpara = c(0.5, 0.06),
                        var.comp.pkpara = c(0.004, 0.00005),
                        error.var.pkpara = 0.000225,
                        fixed.DLTpro,
                        fixed.PDmean,
                        var.comp.PDmean,
                        start.dose = 1,
                        ini.pkpara = c(0.1, 0.005),
                        ini.varcomp = c(0.0007, 0.0000006),
                        ini.error.std = 0.002,
                        ntrial,
                        seed)

{
  if (!dir.exists(PFIM.directory)) {
    stop("PFIM Directory does not exist")
  }
  if (!dir.exists(Func.directory)) {
    stop("Function Directory does not exist")
  }
  if (stop.freq < 6) {
    warning("the value of stop.freq is too low to ensure good operating characteristics.")
  }
  
  set.seed(seed)
  
  source(file.path(Func.directory, "function.r"))
  
  # toxicity limits for up-and-down design
  pL <- target.toxicity / 3
  pM <- (2 * target.toxicity) / 3
  pU <- target.toxicity
  
  fixed.DLTpro <- fixed.DLTpro
  fixed.PDmean <- fixed.PDmean
  var.comp.PDmean <- var.comp.PDmean
  
  # True optimal dose level
  # calculate the true dose-selection values for all dose
  truesel.hat <- dose.select(fixed.DLTpro, fixed.PDmean, W)
  
  # identifying indexes of the doses for which DLT rate below toxicity rate upper bound
  true_index <- which(fixed.DLTpro < target.toxicity &
                      fixed.PDmean >= target.efficacy)
  
  if (length(true_index) == 0) {
    
    true_select_dose_pos <- 1
    
    true_optimal_dose <- 0
    
  } else {
     
    # extracting dose-selection values for these doses
    truesel_candi.hat <- dose.select(fixed.DLTpro[true_index], 
                                     fixed.PDmean[true_index], W)
    
    # the dose level with maximum dose-selection values
    true_select_dose_pos <- which(truesel.hat == max(truesel_candi.hat))
    
    true_optimal_dose <- true.dose[true_select_dose_pos]
    
  }
  
  # True AUC
  auc.true <- model.AUC(true.dose, fixed.pkpara[1], 
                        fixed.pkpara[2], upper.sampling)
  
  # reference AUC
  auc.lim <- auc.true[true_select_dose_pos]
  
  # auc.deviation (qnorm(0.75)=0.67)
  auc.dev <- auc.lim*0.67 

  #################################################################
  # Declaring variables to store the estimates over the simulations
  #################################################################
  
  # store the PK estimates
  pk1.hat <- c()
  pk2.hat <- c()
  varcomp1.hat <- c()
  varcomp2.hat <- c()
  stddev.error.hat <- c()
  
  # store the DLT estimates
  DLT1.hat <- c()
  DLT2.hat <- c()
  
  # store the PD estimates
  PD1.hat <- c()
  PD2.hat <- c()
  PD3.hat <- c()
  PD4.hat <- c()
  
  # store optimal doses
  od <- c()
  # optimal doses from complete analysis
  od.complete <- c()
  
  # counts trials with early stopping (achieving r)
  count.stops <- 0
  
  # up-and-down stage
  count.stops.updown <- 0 # counts trials stopped
  updown.dose <- c()  # dose at which trial stops
  updown.cohort <- c() # cohorts used
  
  # after up-and-down stage
  count.stops.over <- 0
  stoppedearly.dose <- c()  # dose at which trial stops
  stoppedearly.cohort <- c() # cohorts used
  
  # no. of cohorts used in a trial when it stops early
  stopped.at <- c()
  # optimal dose selected in such a trial
  stopped.dose <- c()
  
  # number of subjects allocated to different dose levels
  sub.data <-c()
  sub.num <- c()
  
  source(file.path(PFIM.directory, "model.r"))
  source(file.path(PFIM.directory, "PFIM.r"))
  source(file.path(Func.directory, "function.r"))
  
  #######################################################################
  # Main Program
  #######################################################################
  
  for (l in 1:ntrial) {
    # values to start with
    
    h <- start.dose  # starting dose level
    d <- 1  # to create an index for the subjects
    pk.esti <- ini.pkpara              # initial PK parameter values
    var.comp.esti <- ini.varcomp       # to obtain the optimal time
    error.esti <- ini.error.std        # points through the PFIM
    next.dose <- true.dose[h] # the dose for a trial to start with
    previ.level <- h # to be used in skip calculation
    
    stage <- 1 # stage of a trial
    
    alloc.dose <- c(true.dose[h]) # allocated doses to cohorts in a trial
    data.dlt <- data.frame() # store the DLT responses
    data.eff <- data.frame() # store the PD responses
    auc_used <- c() # store the fitted auc
    
    # loop for different cohorts
    
    for (i in 1:cohort.num) {
      
      stop.iter <- 0
      
      ## generating of PK responses
      k <- 1
      
      OptTime <- tryCatch(PFIM(model.file = "stdin.r", dose = true.dose[h],
                      beta = pk.esti, omega = diag(var.comp.esti),
                      sig.interA = error.esti), 
                 error = function(e)NULL)
      
      if (!is.null(OptTime)) {
        time <- round(OptTime[["opti.fed"]][["opti.tps"]][[1]], digits = 2)
        
        #extracting the optimal time points
      } else{
        stop.iter <- stop.iter + 1
        
      }
      
      # PK response for this cohort
      data.cohort.pk <- outcome.CON(time, cohort.size, fixed.pkpara, 
                                    var.comp.pkpara, error.var.pkpara,
                                    true.dose, h, d, k)
      
      # combine data
      if (d == 1) {
        data.pk <- data.cohort.pk
        
      } else{
        data.pk <- rbind(data.pk, data.cohort.pk)
        
      }
      
      obsdata_CON <- groupedData(formula = conc ~ t | subject, data = data.pk)
      
      ## generating dose-DLT outcomes
      data.cohort.dlt <- outcome.DLT(cohort.size, h, fixed.DLTpro)
      
      data.dlt <- rbind(data.dlt, data.cohort.dlt)
      colnames(data.dlt) <- c('subject1', 'subject2', 'subject3')
      
      ## generating PD outcomes
      data.cohort.eff <- outcome.EFF(cohort.size, h, 
                                     fixed.PDmean, var.comp.PDmean)
      
      data.eff <- rbind(data.eff, data.cohort.eff)
      colnames(data.eff) <- c('subject1', 'subject2', 'subject3')
      
      # subjects' index for the next trial
      d <- d + cohort.size
      
      #######################################################################
      # Dose selection based on up-and-down design
      #######################################################################
      
      if (stage < up.down.num) {
        if ((stage == (up.down.num - 1)) & (all(alloc.dose[1] == alloc.dose))) {
          h <- h + 1
          
        } else{
          prop.toxic <- ((sum(data.dlt)) / (cohort.size * stage))
          
          if (prop.toxic <= pL) {
            if (h == length(true.dose)) {
              h <- h
              
            } else{
              h <- h + 1
            }
          }
          else if ((prop.toxic > pL) & (prop.toxic < pM)) {
            h <- h
            
          }
          
          else if ((prop.toxic >= pM) & (prop.toxic < pU)) {
            if (h == 1) {
              h <- h
              
            } else{
              h <- h - 1
              
            }
            
          } else {
            count.stops.updown <- count.stops.updown + 1
            updown.dose <- c(updown.dose, true.dose[h])
            updown.cohort <- c(updown.cohort, stage)
            break # breaks the trial to the end
            
          }
          
        }
        
      }
      else {
        
        #####################################################################
        # Start of stage > up.down.num: model based approach
        #####################################################################
        
        while (TRUE) {
          # calculate the posterior estimator of PK model parameters
          pk_est <- get.post.PK(obsdata_CON = obsdata_CON,
                                dose = true.dose[h], t = time)
          
          if (pk_est[["theta1_mean"]]>0 & pk_est[["theta2_mean"]]>0) {
            # calculate the new AUC with posterior PK parameters
            auc_new <- model.AUC(x = true.dose, 
                                 v = pk_est[["theta1_mean"]],
                                 cl = pk_est[["theta2_mean"]],
                                 t = upper.sampling)
            
            break
          }
        }
        
        auc_used <- c(auc_used, auc_new[h])
        
        # calculate the posterior estimator of BLRM model parameters
        obsdata_DLT <- Data(x = rep(auc_used, each = cohort.size),
                            y = as.vector(as.matrix(t(data.dlt[-(1:3),]))),
                            cohort = rep(1:(i - 3), each = cohort.size),
                            doseGrid = auc_used,
                            ID = (1:(cohort.size * (i - 3))), placebo = FALSE)
        
        dlt_est <- get.post.BLRM(obsdata_DLT = obsdata_DLT, auc.lim = auc.lim)
        
        # fitting the posterior auc-DLT model
        # for each auc, it has a DLT distribution
        DLT_fit <- matrix(NA, nrow = length(auc_new),
                          ncol = length(dlt_est[["alpha1"]]))
        
        for (a.l in 1:length(auc_new)) {
          DLT_fit[a.l,] <- get.post.DLT(alpha1 = dlt_est[["alpha1"]],
                                        alpha2 = dlt_est[["alpha2"]],
                                        auc = auc_new[a.l],
                                        auc.lim = auc.lim)
          
        }
        
        # dose-DLT rate relationship
        pi.hat <- round(apply(DLT_fit, 1, mean), 2)
        
        # calculate the probability of pi higher than cut-off
        pi.hat_pro <- matrix(NA, nrow = 1, ncol = length(true.dose))
        for (d.l in 1:length(true.dose)) {
          # estimation of posterior function
          pi.hat_pro[, d.l] <- sum(DLT_fit[d.l , ] >= target.toxicity) / 
                               ncol(DLT_fit)
          
        }
        
        # fitting the posterior Emax model
        Emax.jags <- get.post.Emax(obsdata_EFF = as.vector(as.matrix(t(data.eff))),
                                   dose = rep(alloc.dose, each = cohort.size),
                                   n.pat = (cohort.size * i))
        emax_est <- apply(Emax.jags, 2, mean)
        mu.hat <- emax_est[1] + emax_est[3] * true.dose ^ emax_est[4] / 
                  (emax_est[2] ^ emax_est[4] + true.dose ^ emax_est[4])
        
        # calculate the probability of mu lower than cut-off
        mu.hat_pro <- matrix(NA, nrow = 1, ncol = length(true.dose))
        mu.fit <- c()
        
        for (m.l in 1:nrow(Emax.jags)) {
          mu.fit <- rbind(mu.fit, Emax.jags[m.l, 1] + 
                          Emax.jags[m.l, 3] * true.dose ^ Emax.jags[m.l, 4] /
                          (Emax.jags[m.l, 2] ^ Emax.jags[m.l, 4] + 
                             true.dose ^ Emax.jags[m.l, 4]))
          
        }
        
        for (d.l in 1:length(true.dose)) {
          # estimation of posterior function
          mu.hat_pro[, d.l] <- sum(mu.fit[, d.l] < target.efficacy) / 
                                   nrow(mu.fit)
          
        }
        
        # pk estimate for optimal time
        pk.esti <- c(pk_est[["theta1_mean"]], pk_est[["theta2_mean"]])
        var.comp.esti <- c((pk_est[["theta1_sd"]]) ^ 2, 
                           (pk_est[["theta2_sd"]]) ^ 2)
        error.esti <- c(pk_est[["sigma.est_mean"]])
        
        if (pi.hat_pro[1]>tox_cut_off | mu.hat_pro[length(true.dose)]>eff_cut_off){
          
          count.stops.over <- count.stops.over + 1
          h <- 0
          break
          
        }
        
        ###################################################################
        # Dose selection for the next cohort
        ###################################################################
        
        # calculate the dose-selection values for all dose
        sel.hat <- dose.select(pi.hat, mu.hat, W)
        
        # variance of AUC
        dv.esti <- model.AUC.dv(x=true.dose[h], v=pk_est[["theta1_mean"]], 
                                cl=pk_est[["theta2_mean"]], t=upper.sampling)
        dcl.esti <- model.AUC.dcl(x=true.dose[h], v=pk_est[["theta1_mean"]], 
                                  cl=pk_est[["theta2_mean"]], t=upper.sampling)
        c.var.esti <- (((dv.esti^2)*(pk_est[["theta1_sd"]])^2)+
                         ((dcl.esti^2)*(pk_est[["theta2_sd"]])^2))
        
        # relative AUC
        c.std.esti <- sqrt(c.var.esti)
        rel.c <- (auc_new - auc.lim)/c.std.esti  
        
        # identifying indexes of the doses for which both conditions are met
        index <- which(pi.hat < target.toxicity & mu.hat >= target.efficacy
                       & pi.hat_pro <= tox_cut_off & mu.hat_pro <= eff_cut_off 
                       & rel.c<=auc.dev)
        
        if (length(index) == 0) {
          select_dose_pos <- h
          
        } else {
          # extracting dose-selection values for these doses
          sel_candi.hat <- dose.select(pi.hat[index], mu.hat[index], W)
          
          # the dose level with maximum dose-selection values
          select_dose_pos <- which(sel.hat == max(sel_candi.hat))
          
        }
        
        
        ##################################################################
        # Restriction on dose skipping
        ##################################################################
        
        if (stage <= cohort.num) {
          # with <=, basically restriction on the last stage
          # leave it (just <) if we don't want restriction on the last stage
          if ((select_dose_pos - previ.level) > 0) {
            h <- previ.level + dose.skip
            
          } else if (((select_dose_pos - previ.level) < 0)) {
            h <- previ.level - dose.skip
            
          } else{
            h <- select_dose_pos
            
          }
        }
        
      }
      # end of stage>up.down.num
      
      # selected dose for the next cohort
      next.dose <- true.dose[h]
      
      # current dose level is going to be the previous dose level in the next stage
      previ.level <- h
      
      # counts stage for a trial
      stage <- stage + 1
      
      # determines the maximum frequent for the allocated doses to cohorts
      alloc.dose <- c(alloc.dose, true.dose[h])
      max.freq <- max(table(alloc.dose))
      
      if (max.freq == stop.freq & stage <= (cohort.num - 1)) {
        # break the current trial
        count.stops <- count.stops + 1
        stopped.at <- c(stopped.at, length(alloc.dose))
        stopped.dose <- c(stopped.dose, true.dose[h])
        break
        
      }
      
    }
    # end of the loop with index 'i'
    
    if (stage > up.down.num) {
      # store the model estimates
      pk1.hat <- c(pk1.hat, pk_est[["theta1_mean"]])
      pk2.hat <- c(pk2.hat, pk_est[["theta2_mean"]])
      
      varcomp1.hat <- c(varcomp1.hat, (pk_est[["theta1_sd"]]) ^ 2)
      varcomp2.hat <- c(varcomp2.hat, (pk_est[["theta2_sd"]]) ^ 2)
      stddev.error.hat <- c(stddev.error.hat, pk_est[["sigma.est_mean"]])
      
      DLT1.hat <- c(DLT1.hat, mean(dlt_est[["alpha1"]]))
      DLT2.hat <- c(DLT2.hat, mean(dlt_est[["alpha2"]]))
      
      PD1.hat <- c(PD1.hat, emax_est[1])
      PD2.hat <- c(PD2.hat, emax_est[2])
      PD3.hat <- c(PD3.hat, emax_est[3])
      PD4.hat <- c(PD4.hat, emax_est[4])
      
      # recommended dose for the next stage
      if (stage == cohort.num + 1) {
        # OD from complete analysis
        od.complete <- c(od.complete, true.dose[h])
        
        # removing the recommended dose at the last stage
        # from the list of the allocated doses
        alloc.sim <- capture.output(alloc.dose[-length(alloc.dose)])
        
        
      } else{
        alloc.sim <- capture.output(alloc.dose)
        
      }
      
      # records all the ODs: early and complete
      od <- c(od, true.dose[h])
      
      sub.data <- rbind(sub.data, count(data.pk, dose))
      
      
    }
    # end of stage>up.down.num
    
  }
  # end of the loop with index 'l'
  
  # calculate the percentage of subjects for allocated dose levels
  sub.data$n <- sub.data$n / 3
  sub.num <- sub.data %>%
    group_by(dose) %>%
    summarise(total_n = sum(n))
  
  out = list(
    
    PK.Para1 = pk1.hat, 
    PK.Para2 = pk2.hat,
    
    PK.Para.var1 = varcomp1.hat, 
    PK.Para.var2 = varcomp2.hat,
    PK.Para.std = stddev.error.hat,
    
    DLT.Para1 = DLT1.hat,
    DLT.Para2 = DLT2.hat, 
    
    PD.Para1 = PD1.hat, 
    PD.Para2 = PD2.hat,
    PD.Para3 = PD3.hat, 
    PD.Para4 = PD4.hat,
    
    num.stops.updown = count.stops.updown,
    num.stops.over = count.stops.over,
    num.stops = count.stops,
    true.opt.dose =  true_optimal_dose,
    sim.opt.dose = od,
    per.opt.dose = table(od) / (ntrial - count.stops.updown - count.stops.over),
    all.dose = sub.num$dose,
    all.num = sub.num$total_n,
    per.all.sub = sub.num$total_n / sum(sub.num$total_n) 
  )
  
}


