
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
                        , nameExcel = 'sim7a_HD_20230305'
                        , seed = 20230305
)

#####

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
                        , nameExcel = 'sim7b_HD_20230305'
                        , seed = 20230305
)

#####

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
                        , nameExcel = 'sim7c_HD_20230305'
                        , seed = 20230305
)

#####

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
                        , nameExcel = 'sim7d_HD_20230305'
                        , seed = 20230305
)

#####

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
                        , nameExcel = 'sim7e_HD_20230305'
                        , seed = 20230305
)

#####

resMIf <- runSimulation(S = 100000L
                        , S_at_a_time = 10000L
                        , N = 1500L
                        , cfs = 2L
                        , m = 100L
                        , model.x = 'rnorm(n = N, mean = 0, sd = 1)'
                        , model.u = 'rnorm(n = N, mean = 0, sd = 1)'
                        , model.eps = 'rnorm(n = N, mean = 0, sd = 1)'
                        , model.y = '-0.5 + 0.5 * pop.x + pop.eps'
                        , model.logit.pi = '1.76 - 0.33 * pop.u + 0.1 * pop.y[ ,s]'
                        , model.logit.phi = '1 + 0.5 * pop.x + 0.5 * pop.u'
                        , sigma.est = 'Hajekplus'
                        , draw.u = 'hotdeck'
                        , informative = TRUE
                        , runKY = TRUE
                        , fpcKY = TRUE
                        , approxKY = TRUE
                        , runTradMI = TRUE
                        , writeExcel = TRUE
                        , nameExcel = 'sim7f_HD_20230305'
                        , seed = 20230305
)

#####

resMIg <- runSimulation(S = 100000L
                        , S_at_a_time = 10000L
                        , N = 2500L
                        , cfs = 2L
                        , m = 100L
                        , model.x = 'rnorm(n = N, mean = 0, sd = 1)'
                        , model.u = 'rnorm(n = N, mean = 0, sd = 1)'
                        , model.eps = 'rnorm(n = N, mean = 0, sd = 1)'
                        , model.y = '-0.5 + 0.5 * pop.x + pop.eps'
                        , model.logit.pi = '0.07 - 0.33 * pop.u + 0.1 * pop.y[ ,s]'
                        , model.logit.phi = '1 + 0.5 * pop.x + 0.5 * pop.u'
                        , sigma.est = 'Hajekplus'
                        , draw.u = 'hotdeck'
                        , informative = TRUE
                        , runKY = TRUE
                        , fpcKY = TRUE
                        , approxKY = TRUE
                        , runTradMI = TRUE
                        , writeExcel = TRUE
                        , nameExcel = 'sim7g_HD_20230305'
                        , seed = 20230305
)

#####

resMIh <- runSimulation(S = 100000L
                        , S_at_a_time = 10000L
                        , N = 5000L
                        , cfs = 2L
                        , m = 100L
                        , model.x = 'rnorm(n = N, mean = 0, sd = 1)'
                        , model.u = 'rnorm(n = N, mean = 0, sd = 1)'
                        , model.eps = 'rnorm(n = N, mean = 0, sd = 1)'
                        , model.y = '-0.5 + 0.5 * pop.x + pop.eps'
                        , model.logit.pi = '-1.07 - 0.33 * pop.u + 0.1 * pop.y[ ,s]'
                        , model.logit.phi = '1 + 0.5 * pop.x + 0.5 * pop.u'
                        , sigma.est = 'Hajekplus'
                        , draw.u = 'hotdeck'
                        , informative = TRUE
                        , runKY = TRUE
                        , fpcKY = TRUE
                        , approxKY = TRUE
                        , runTradMI = TRUE
                        , writeExcel = TRUE
                        , nameExcel = 'sim7h_HD_20230305'
                        , seed = 20230305
)

#####

resMIi <- runSimulation(S = 100000L
                        , S_at_a_time = 10000L
                        , N = 10000L
                        , cfs = 2L
                        , m = 100L
                        , model.x = 'rnorm(n = N, mean = 0, sd = 1)'
                        , model.u = 'rnorm(n = N, mean = 0, sd = 1)'
                        , model.eps = 'rnorm(n = N, mean = 0, sd = 1)'
                        , model.y = '-0.5 + 0.5 * pop.x + pop.eps'
                        , model.logit.pi = '-1.93 - 0.33 * pop.u + 0.1 * pop.y[ ,s]'
                        , model.logit.phi = '1 + 0.5 * pop.x + 0.5 * pop.u'
                        , sigma.est = 'Hajekplus'
                        , draw.u = 'hotdeck'
                        , informative = TRUE
                        , runKY = TRUE
                        , fpcKY = TRUE
                        , approxKY = TRUE
                        , runTradMI = TRUE
                        , writeExcel = TRUE
                        , nameExcel = 'sim7i_HD_20230305'
                        , seed = 20230305
)

#####

resMIj <- runSimulation(S = 100000L
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
                        , nameExcel = 'sim7j_HD_20230305'
                        , seed = 20230305
)



##### Method of Kim & Yang

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
                        , sigma.est = 'Hajekplus'
                        , draw.u = 'hotdeck'
                        , informative = TRUE
                        , runKY = TRUE
                        , fpcKY = TRUE
                        , approxKY = FALSE
                        , writeExcel = TRUE
                        , nameExcel = 'sim7a_KY_20230305'
                        , seed = 20230305
)

#####

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
                        , sigma.est = 'Hajekplus'
                        , draw.u = 'hotdeck'
                        , informative = TRUE
                        , runKY = TRUE
                        , fpcKY = TRUE
                        , approxKY = FALSE
                        , writeExcel = TRUE
                        , nameExcel = 'sim7b_KY_20230305'
                        , seed = 20230305
)

#####

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
                        , sigma.est = 'Hajekplus'
                        , draw.u = 'hotdeck'
                        , informative = TRUE
                        , runKY = TRUE
                        , fpcKY = TRUE
                        , approxKY = FALSE
                        , writeExcel = TRUE
                        , nameExcel = 'sim7c_KY_20230305'
                        , seed = 20230305
)

#####

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
                        , sigma.est = 'Hajekplus'
                        , draw.u = 'hotdeck'
                        , informative = TRUE
                        , runKY = TRUE
                        , fpcKY = TRUE
                        , approxKY = FALSE
                        , writeExcel = TRUE
                        , nameExcel = 'sim7d_KY_20230305'
                        , seed = 20230305
)

#####

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
                        , sigma.est = 'Hajekplus'
                        , draw.u = 'hotdeck'
                        , informative = TRUE
                        , runKY = TRUE
                        , fpcKY = TRUE
                        , approxKY = FALSE
                        , writeExcel = TRUE
                        , nameExcel = 'sim7e_KY_20230305'
                        , seed = 20230305
)

#####

resMIf.KY <- runSimulation(S = 10000L
                        , S_at_a_time = 10000L
                        , N = 1500L
                        , cfs = 2L
                        , m = 100L
                        , model.x = 'rnorm(n = N, mean = 0, sd = 1)'
                        , model.u = 'rnorm(n = N, mean = 0, sd = 1)'
                        , model.eps = 'rnorm(n = N, mean = 0, sd = 1)'
                        , model.y = '-0.5 + 0.5 * pop.x + pop.eps'
                        , model.logit.pi = '1.76 - 0.33 * pop.u + 0.1 * pop.y[ ,s]'
                        , model.logit.phi = '1 + 0.5 * pop.x + 0.5 * pop.u'
                        , sigma.est = 'Hajekplus'
                        , draw.u = 'hotdeck'
                        , informative = TRUE
                        , runKY = TRUE
                        , fpcKY = TRUE
                        , approxKY = FALSE
                        , writeExcel = TRUE
                        , nameExcel = 'sim7f_KY_20230305'
                        , seed = 20230305
)

#####

resMIg.KY <- runSimulation(S = 10000L
                        , S_at_a_time = 10000L
                        , N = 2500L
                        , cfs = 2L
                        , m = 100L
                        , model.x = 'rnorm(n = N, mean = 0, sd = 1)'
                        , model.u = 'rnorm(n = N, mean = 0, sd = 1)'
                        , model.eps = 'rnorm(n = N, mean = 0, sd = 1)'
                        , model.y = '-0.5 + 0.5 * pop.x + pop.eps'
                        , model.logit.pi = '0.07 - 0.33 * pop.u + 0.1 * pop.y[ ,s]'
                        , model.logit.phi = '1 + 0.5 * pop.x + 0.5 * pop.u'
                        , sigma.est = 'Hajekplus'
                        , draw.u = 'hotdeck'
                        , informative = TRUE
                        , runKY = TRUE
                        , fpcKY = TRUE
                        , approxKY = FALSE
                        , writeExcel = TRUE
                        , nameExcel = 'sim7g_KY_20230305'
                        , seed = 20230305
)

#####

resMIh.KY <- runSimulation(S = 10000L
                        , S_at_a_time = 10000L
                        , N = 5000L
                        , cfs = 2L
                        , m = 100L
                        , model.x = 'rnorm(n = N, mean = 0, sd = 1)'
                        , model.u = 'rnorm(n = N, mean = 0, sd = 1)'
                        , model.eps = 'rnorm(n = N, mean = 0, sd = 1)'
                        , model.y = '-0.5 + 0.5 * pop.x + pop.eps'
                        , model.logit.pi = '-1.07 - 0.33 * pop.u + 0.1 * pop.y[ ,s]'
                        , model.logit.phi = '1 + 0.5 * pop.x + 0.5 * pop.u'
                        , sigma.est = 'Hajekplus'
                        , draw.u = 'hotdeck'
                        , informative = TRUE
                        , runKY = TRUE
                        , fpcKY = TRUE
                        , approxKY = FALSE
                        , writeExcel = TRUE
                        , nameExcel = 'sim7h_KY_20230305'
                        , seed = 20230305
)

#####

resMIi.KY <- runSimulation(S = 10000L
                        , S_at_a_time = 10000L
                        , N = 10000L
                        , cfs = 2L
                        , m = 100L
                        , model.x = 'rnorm(n = N, mean = 0, sd = 1)'
                        , model.u = 'rnorm(n = N, mean = 0, sd = 1)'
                        , model.eps = 'rnorm(n = N, mean = 0, sd = 1)'
                        , model.y = '-0.5 + 0.5 * pop.x + pop.eps'
                        , model.logit.pi = '-1.93 - 0.33 * pop.u + 0.1 * pop.y[ ,s]'
                        , model.logit.phi = '1 + 0.5 * pop.x + 0.5 * pop.u'
                        , sigma.est = 'Hajekplus'
                        , draw.u = 'hotdeck'
                        , informative = TRUE
                        , runKY = TRUE
                        , fpcKY = TRUE
                        , approxKY = FALSE
                        , writeExcel = TRUE
                        , nameExcel = 'sim7i_KY_20230305'
                        , seed = 20230305
)

#####

resMIj.KY <- runSimulation(S = 10000L
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
                        , approxKY = FALSE
                        , writeExcel = TRUE
                        , nameExcel = 'sim7j_KY_20230305'
                        , seed = 20230305
)
