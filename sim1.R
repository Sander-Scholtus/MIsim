
rm(list = ls())
gc()

library(openxlsx)
library(sampling)
library(MASS)
library(foreach)
library(doParallel)

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

source('sim_functions.R')

#####

resMIa <- runSimulation(S = 100000L
                        , S_at_a_time = 10000L
                        , N = 1500L
                        , cfs = 2L
                        , m = 100L
                        , model.x = 'rnorm(n = N, mean = 0, sd = 1)'
                        , model.u = 'rnorm(n = N, mean = 0, sd = 1)'
                        , model.eps = 'rnorm(n = N, mean = 0, sd = 1)'
                        , model.y = '-0.5 + 0.5 * pop.x + pop.eps'
                        , model.logit.pi = '2.00 - 0.33 * pop.u + 0.5 * pop.y[ ,s]'
                        , model.logit.phi = '1 + 0.5 * pop.x + 0.5 * pop.u'
                        , sigma.est = 'Hajekplus'
                        , draw.u = 'hotdeck'
                        , informative = TRUE
                        , runKY = TRUE
                        , fpcKY = TRUE
                        , approxKY = TRUE
                        , runTradMI = TRUE
                        , writeExcel = TRUE
                        , nameExcel = 'sim1a_HD_20230305'
                        , seed = 20230305
)

resMIb <- runSimulation(S = 100000L
                        , S_at_a_time = 10000L
                        , N = 2500L
                        , cfs = 2L
                        , m = 100L
                        , model.x = 'rnorm(n = N, mean = 0, sd = 1)'
                        , model.u = 'rnorm(n = N, mean = 0, sd = 1)'
                        , model.eps = 'rnorm(n = N, mean = 0, sd = 1)'
                        , model.y = '-0.5 + 0.5 * pop.x + pop.eps'
                        , model.logit.pi = '0.25 - 0.33 * pop.u + 0.5 * pop.y[ ,s]'
                        , model.logit.phi = '1 + 0.5 * pop.x + 0.5 * pop.u'
                        , sigma.est = 'Hajekplus'
                        , draw.u = 'hotdeck'
                        , informative = TRUE
                        , runKY = TRUE
                        , fpcKY = TRUE
                        , approxKY = TRUE
                        , runTradMI = TRUE
                        , writeExcel = TRUE
                        , nameExcel = 'sim1b_HD_20230305'
                        , seed = 20230305
)

resMIc <- runSimulation(S = 100000L
                        , S_at_a_time = 10000L
                        , N = 5000L
                        , cfs = 2L
                        , m = 100L
                        , model.x = 'rnorm(n = N, mean = 0, sd = 1)'
                        , model.u = 'rnorm(n = N, mean = 0, sd = 1)'
                        , model.eps = 'rnorm(n = N, mean = 0, sd = 1)'
                        , model.y = '-0.5 + 0.5 * pop.x + pop.eps'
                        , model.logit.pi = '-0.96 - 0.33 * pop.u + 0.5 * pop.y[ ,s]'
                        , model.logit.phi = '1 + 0.5 * pop.x + 0.5 * pop.u'
                        , sigma.est = 'Hajekplus'
                        , draw.u = 'hotdeck'
                        , informative = TRUE
                        , runKY = TRUE
                        , fpcKY = TRUE
                        , approxKY = TRUE
                        , runTradMI = TRUE
                        , writeExcel = TRUE
                        , nameExcel = 'sim1c_HD_20230305'
                        , seed = 20230305
)

resMId <- runSimulation(S = 100000L
                        , S_at_a_time = 10000L
                        , N = 10000L
                        , cfs = 2L
                        , m = 100L
                        , model.x = 'rnorm(n = N, mean = 0, sd = 1)'
                        , model.u = 'rnorm(n = N, mean = 0, sd = 1)'
                        , model.eps = 'rnorm(n = N, mean = 0, sd = 1)'
                        , model.y = '-0.5 + 0.5 * pop.x + pop.eps'
                        , model.logit.pi = '-1.85 - 0.33 * pop.u + 0.5 * pop.y[ ,s]'
                        , model.logit.phi = '1 + 0.5 * pop.x + 0.5 * pop.u'
                        , sigma.est = 'Hajekplus'
                        , draw.u = 'hotdeck'
                        , informative = TRUE
                        , runKY = TRUE
                        , fpcKY = TRUE
                        , approxKY = TRUE
                        , runTradMI = TRUE
                        , writeExcel = TRUE
                        , nameExcel = 'sim1d_HD_20230305'
                        , seed = 20230305
)

resMIe <- runSimulation(S = 100000L
                        , S_at_a_time = 5000L
                        , N = 50000L
                        , cfs = 2L
                        , m = 100L
                        , model.x = 'rnorm(n = N, mean = 0, sd = 1)'
                        , model.u = 'rnorm(n = N, mean = 0, sd = 1)'
                        , model.eps = 'rnorm(n = N, mean = 0, sd = 1)'
                        , model.y = '-0.5 + 0.5 * pop.x + pop.eps'
                        , model.logit.pi = '-3.60 - 0.33 * pop.u + 0.5 * pop.y[ ,s]'
                        , model.logit.phi = '1 + 0.5 * pop.x + 0.5 * pop.u'
                        , sigma.est = 'Hajekplus'
                        , draw.u = 'hotdeck'
                        , informative = TRUE
                        , runKY = TRUE
                        , fpcKY = TRUE
                        , approxKY = TRUE
                        , runTradMI = TRUE
                        , writeExcel = TRUE
                        , nameExcel = 'sim1e_HD_20230305'
                        , seed = 20230305
)

# check: original scenario from Kim & Yang's paper
resMIf <- runSimulation(S = 100000L
                        , S_at_a_time = 5000L
                        , N = 50000L
                        , cfs = 2L
                        , m = 100L
                        , model.x = 'rnorm(n = N, mean = 0, sd = 1)'
                        , model.u = 'rnorm(n = N, mean = 0, sd = 1)'
                        , model.eps = 'rnorm(n = N, mean = 0, sd = 1)'
                        , model.y = '-0.5 + 0.5 * pop.x + pop.eps'
                        , model.logit.pi = '-3.66 - 0.33 * pop.u + 0.1 * pop.y[ ,s]'
                        , model.logit.phi = '1 + 0.5 * pop.x + 0.5 * pop.u'
                        , sigma.est = 'Hajekplus'
                        , draw.u = 'hotdeck'
                        , informative = TRUE
                        , runKY = TRUE
                        , fpcKY = TRUE
                        , approxKY = TRUE
                        , runTradMI = TRUE
                        , writeExcel = TRUE
                        , nameExcel = 'sim1f_HD_20230305'
                        , seed = 20230305
)



#####

resMIa.KY <- runSimulation(S = 10000L
                           , S_at_a_time = 10000L
                           , N = 1500L
                           , cfs = 2L
                           , m = 100L
                           , model.x = 'rnorm(n = N, mean = 0, sd = 1)'
                           , model.u = 'rnorm(n = N, mean = 0, sd = 1)'
                           , model.eps = 'rnorm(n = N, mean = 0, sd = 1)'
                           , model.y = '-0.5 + 0.5 * pop.x + pop.eps'
                           , model.logit.pi = '2.00 - 0.33 * pop.u + 0.5 * pop.y[ ,s]'
                           , model.logit.phi = '1 + 0.5 * pop.x + 0.5 * pop.u'
                           , informative = TRUE
                           , runKY = TRUE
                           , fpcKY = TRUE
                           , writeExcel = TRUE
                           , nameExcel = 'sim1a_KY_20230305'
                           , seed = 20230305
)

resMIb.KY <- runSimulation(S = 10000L
                           , S_at_a_time = 10000L
                           , N = 2500L
                           , cfs = 2L
                           , m = 100L
                           , model.x = 'rnorm(n = N, mean = 0, sd = 1)'
                           , model.u = 'rnorm(n = N, mean = 0, sd = 1)'
                           , model.eps = 'rnorm(n = N, mean = 0, sd = 1)'
                           , model.y = '-0.5 + 0.5 * pop.x + pop.eps'
                           , model.logit.pi = '0.25 - 0.33 * pop.u + 0.5 * pop.y[ ,s]'
                           , model.logit.phi = '1 + 0.5 * pop.x + 0.5 * pop.u'
                           , informative = TRUE
                           , runKY = TRUE
                           , fpcKY = TRUE
                           , writeExcel = TRUE
                           , nameExcel = 'sim1b_KY_20230305'
                           , seed = 20230305
)

resMIc.KY <- runSimulation(S = 10000L
                           , S_at_a_time = 10000L
                           , N = 5000L
                           , cfs = 2L
                           , m = 100L
                           , model.x = 'rnorm(n = N, mean = 0, sd = 1)'
                           , model.u = 'rnorm(n = N, mean = 0, sd = 1)'
                           , model.eps = 'rnorm(n = N, mean = 0, sd = 1)'
                           , model.y = '-0.5 + 0.5 * pop.x + pop.eps'
                           , model.logit.pi = '-0.96 - 0.33 * pop.u + 0.5 * pop.y[ ,s]'
                           , model.logit.phi = '1 + 0.5 * pop.x + 0.5 * pop.u'
                           , informative = TRUE
                           , runKY = TRUE
                           , fpcKY = TRUE
                           , writeExcel = TRUE
                           , nameExcel = 'sim1c_KY_20230305'
                           , seed = 20230305
)

resMId.KY <- runSimulation(S = 10000L
                           , S_at_a_time = 10000L
                           , N = 10000L
                           , cfs = 2L
                           , m = 100L
                           , model.x = 'rnorm(n = N, mean = 0, sd = 1)'
                           , model.u = 'rnorm(n = N, mean = 0, sd = 1)'
                           , model.eps = 'rnorm(n = N, mean = 0, sd = 1)'
                           , model.y = '-0.5 + 0.5 * pop.x + pop.eps'
                           , model.logit.pi = '-1.85 - 0.33 * pop.u + 0.5 * pop.y[ ,s]'
                           , model.logit.phi = '1 + 0.5 * pop.x + 0.5 * pop.u'
                           , informative = TRUE
                           , runKY = TRUE
                           , fpcKY = TRUE
                           , writeExcel = TRUE
                           , nameExcel = 'sim1d_KY_20230305'
                           , seed = 20230305
)

resMIe.KY <- runSimulation(S = 10000L
                           , S_at_a_time = 5000L
                           , N = 50000L
                           , cfs = 2L
                           , m = 100L
                           , model.x = 'rnorm(n = N, mean = 0, sd = 1)'
                           , model.u = 'rnorm(n = N, mean = 0, sd = 1)'
                           , model.eps = 'rnorm(n = N, mean = 0, sd = 1)'
                           , model.y = '-0.5 + 0.5 * pop.x + pop.eps'
                           , model.logit.pi = '-3.60 - 0.33 * pop.u + 0.5 * pop.y[ ,s]'
                           , model.logit.phi = '1 + 0.5 * pop.x + 0.5 * pop.u'
                           , informative = TRUE
                           , runKY = TRUE
                           , fpcKY = TRUE
                           , writeExcel = TRUE
                           , nameExcel = 'sim1e_KY_20230305'
                           , seed = 20230305
)

# check: original scenario from Kim & Yang's paper
resMIf.KY <- runSimulation(S = 10000L
                           , S_at_a_time = 5000L
                           , N = 50000L
                           , cfs = 2L
                           , m = 100L
                           , model.x = 'rnorm(n = N, mean = 0, sd = 1)'
                           , model.u = 'rnorm(n = N, mean = 0, sd = 1)'
                           , model.eps = 'rnorm(n = N, mean = 0, sd = 1)'
                           , model.y = '-0.5 + 0.5 * pop.x + pop.eps'
                           , model.logit.pi = '-3.66 - 0.33 * pop.u + 0.1 * pop.y[ ,s]'
                           , model.logit.phi = '1 + 0.5 * pop.x + 0.5 * pop.u'
                           , informative = TRUE
                           , runKY = TRUE
                           , fpcKY = TRUE
                           , writeExcel = TRUE
                           , nameExcel = 'sim1f_KY_20230305'
                           , seed = 20230305
)

