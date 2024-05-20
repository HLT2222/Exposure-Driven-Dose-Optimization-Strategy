library(mvtnorm)
library(nlme)
library(igraph)
library(arm)
library(brms)
library(ggplot2)
library(crmPack)
library(rjags)
library(R2jags)
library(MASS)
library(combinat)
library(dplyr)

rm(list = ls()) 

setwd(" ")
source('simple_simulation.r')
source('simulation.r')

scenarios_list <- list(scenario_1 = list(fixed.DLTpro <- c(0.01, 0.03, 0.05, 0.06, 0.08),
                                         fixed.PDmean <- c(0.03, 0.08, 0.13, 0.15, 0.17),
                                         var.comp.PDmean <- 0.1),
                       scenario_2 = list(fixed.DLTpro <- c(0.01, 0.02, 0.10, 0.17, 0.30),
                                         fixed.PDmean <- c(0.10, 0.21, 0.32, 0.45, 0.51),
                                         var.comp.PDmean <- 0.1),
                       scenario_3 = list(fixed.DLTpro <- c(0.03, 0.06, 0.10, 0.18, 0.32),
                                         fixed.PDmean <- c(0.08, 0.20, 0.45, 0.53, 0.55),
                                         var.comp.PDmean <- 0.1),
                       scenario_4 = list(fixed.DLTpro <- c(0.02, 0.02, 0.16, 0.34, 0.40),
                                         fixed.PDmean <- c(0.00, 0.03, 0.24, 0.45, 0.56),
                                         var.comp.PDmean <- 0.1),
                       scenario_5 = list(fixed.DLTpro <- c(0.05, 0.12, 0.18, 0.31, 0.45),
                                         fixed.PDmean <- c(0.12, 0.21, 0.35, 0.46, 0.49),
                                         var.comp.PDmean <- 0.1),
                       scenario_6 = list(fixed.DLTpro <- c(0.02, 0.15, 0.35, 0.40, 0.50),
                                         fixed.PDmean <- c(0.20, 0.44, 0.61, 0.62, 0.61),
                                         var.comp.PDmean <- 0.1),
                       scenario_7 = list(fixed.DLTpro <- c(0.04, 0.18, 0.32, 0.44, 0.50),
                                         fixed.PDmean <- c(0.00, 0.30, 0.45, 0.60, 0.55),
                                         var.comp.PDmean <- 0.1),
                       scenario_8 = list(fixed.DLTpro <- c(0.06, 0.31, 0.50, 0.60, 0.75),
                                         fixed.PDmean <- c(0.40, 0.43, 0.44, 0.45, 0.47),
                                         var.comp.PDmean <- 0.1),
                       scenario_9 = list(fixed.DLTpro <- c(0.40, 0.48, 0.57, 0.65, 0.74),
                                         fixed.PDmean <- c(0.05, 0.17, 0.31, 0.36, 0.40),
                                         var.comp.PDmean <- 0.1))


model_out <- full_model(PFIM.directory = "PFIM/PFIM4.0",
                             Func.directory = " ",
                             true.dose = seq(0.1, 0.9, 0.2),
                             cohort.size = 3,
                             cohort.num = 20,
                             W = 0.6,
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
                             fixed.DLTpro = scenarios_list[["scenario_1"]][[1]],
                             fixed.PDmean = scenarios_list[["scenario_1"]][[2]],
                             var.comp.PDmean = scenarios_list[["scenario_1"]][[3]],
                             start.dose = 1,
                             ini.pkpara = c(0.1, 0.005),
                             ini.varcomp = c(0.0007, 0.0000006),
                             ini.error.std = 0.002,
                             ntrial = 50,
                             seed = 12345)




