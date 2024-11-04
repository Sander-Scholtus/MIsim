
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
                        , N = 1000L
                        , cfs = 2L
                        , m = 100L
                        , model.x = 'rnorm(n = N, mean = 0, sd = 1)'
                        , model.u = 'rnorm(n = N, mean = 0, sd = 1)'
                        , model.eps = 'rnorm(n = N, mean = 0, sd = 1)'
                        , model.y = '-0.5 + 0.5 * pop.x + pop.eps'
                        , model.pi = '0.35'
                        , model.logit.phi = '0 + 0.5 * pop.x + 0.5 * pop.u'
                        , sigma.est = 'Hajekplus'
                        , draw.u = 'hotdeck'
                        , informative = FALSE
                        , runKY = TRUE
                        , fpcKY = TRUE
                        , approxKY = TRUE
                        , runTradMI = FALSE
                        , writeExcel = TRUE
                        , nameExcel = 'sim2a_HD_20230305'
                        , seed = 20230305
)

resMIb <- runSimulation(S = 100000L
                        , S_at_a_time = 10000L
                        , N = 1000L
                        , cfs = 2L
                        , m = 100L
                        , model.x = 'rnorm(n = N, mean = 0, sd = 1)'
                        , model.u = 'rnorm(n = N, mean = 0, sd = 1)'
                        , model.eps = 'rnorm(n = N, mean = 0, sd = 1)'
                        , model.y = '-0.5 + 0.5 * pop.x + pop.eps'
                        , model.pi = '0.70'
                        , model.logit.phi = '-1 + 0.5 * pop.x + 0.5 * pop.u'
                        , sigma.est = 'Hajekplus'
                        , draw.u = 'hotdeck'
                        , informative = FALSE
                        , runKY = TRUE
                        , fpcKY = TRUE
                        , approxKY = TRUE
                        , runTradMI = FALSE
                        , writeExcel = TRUE
                        , nameExcel = 'sim2b_HD_20230305'
                        , seed = 20230305
)

resMIb_nofpc <- runSimulation(S = 100000L
                              , S_at_a_time = 10000L
                              , N = 1000L
                              , cfs = 2L
                              , m = 100L
                              , model.x = 'rnorm(n = N, mean = 0, sd = 1)'
                              , model.u = 'rnorm(n = N, mean = 0, sd = 1)'
                              , model.eps = 'rnorm(n = N, mean = 0, sd = 1)'
                              , model.y = '-0.5 + 0.5 * pop.x + pop.eps'
                              , model.pi = '0.70'
                              , model.logit.phi = '-1 + 0.5 * pop.x + 0.5 * pop.u'
                              , sigma.est = 'Hajekplus'
                              , draw.u = 'hotdeck'
                              , informative = FALSE
                              , runKY = TRUE
                              , fpcKY = FALSE
                              , approxKY = TRUE
                              , runTradMI = FALSE
                              , writeExcel = TRUE
                              , nameExcel = 'sim2b_HD_nofpc_20230305'
                              , seed = 20230305
)

# check: original scenario from Kim & Yang's paper
resMIc <- runSimulation(S = 100000L
                        , S_at_a_time = 5000L
                        , N = 50000L
                        , cfs = 2L
                        , m = 100L
                        , model.x = 'rnorm(n = N, mean = 0, sd = 1)'
                        , model.u = 'rnorm(n = N, mean = 0, sd = 1)'
                        , model.eps = 'rnorm(n = N, mean = 0, sd = 1)'
                        , model.y = '-0.5 + 0.5 * pop.x + pop.eps'
                        , model.logit.pi = '-4 - 0.5 * pop.x'
                        , model.logit.phi = '1 + 0.5 * pop.x + 0.5 * pop.u'
                        , sigma.est = 'Hajekplus'
                        , draw.u = 'hotdeck'
                        , informative = FALSE
                        , runKY = TRUE
                        , fpcKY = TRUE
                        , approxKY = TRUE
                        , runTradMI = FALSE
                        , writeExcel = TRUE
                        , nameExcel = 'sim2c_HD_20230305'
                        , seed = 20230305
)


#####

resMIa.KY <- runSimulation(S = 10000L
                           , S_at_a_time = 10000L
                           , N = 1000L
                           , cfs = 2L
                           , m = 100L
                           , model.x = 'rnorm(n = N, mean = 0, sd = 1)'
                           , model.u = 'rnorm(n = N, mean = 0, sd = 1)'
                           , model.eps = 'rnorm(n = N, mean = 0, sd = 1)'
                           , model.y = '-0.5 + 0.5 * pop.x + pop.eps'
                           , model.pi = '0.35'
                           , model.logit.phi = '0 + 0.5 * pop.x + 0.5 * pop.u'
                           , sigma.est = 'Hajekplus'
                           , draw.u = 'hotdeck'
                           , informative = FALSE
                           , runKY = TRUE
                           , fpcKY = TRUE
                           , approxKY = FALSE
                           , runTradMI = FALSE
                           , writeExcel = TRUE
                           , nameExcel = 'sim2a_KY_20230305'
                           , seed = 20230305
)

resMIb.KY <- runSimulation(S = 10000L
                           , S_at_a_time = 10000L
                           , N = 1000L
                           , cfs = 2L
                           , m = 100L
                           , model.x = 'rnorm(n = N, mean = 0, sd = 1)'
                           , model.u = 'rnorm(n = N, mean = 0, sd = 1)'
                           , model.eps = 'rnorm(n = N, mean = 0, sd = 1)'
                           , model.y = '-0.5 + 0.5 * pop.x + pop.eps'
                           , model.pi = '0.70'
                           , model.logit.phi = '-1 + 0.5 * pop.x + 0.5 * pop.u'
                           , sigma.est = 'Hajekplus'
                           , draw.u = 'hotdeck'
                           , informative = FALSE
                           , runKY = TRUE
                           , fpcKY = TRUE
                           , approxKY = FALSE
                           , runTradMI = FALSE
                           , writeExcel = TRUE
                           , nameExcel = 'sim2b_KY_20230305'
                           , seed = 20230305
)

resMIb_nofpc.KY <- runSimulation(S = 10000L
                                 , S_at_a_time = 10000L
                                 , N = 1000L
                                 , cfs = 2L
                                 , m = 100L
                                 , model.x = 'rnorm(n = N, mean = 0, sd = 1)'
                                 , model.u = 'rnorm(n = N, mean = 0, sd = 1)'
                                 , model.eps = 'rnorm(n = N, mean = 0, sd = 1)'
                                 , model.y = '-0.5 + 0.5 * pop.x + pop.eps'
                                 , model.pi = '0.70'
                                 , model.logit.phi = '-1 + 0.5 * pop.x + 0.5 * pop.u'
                                 , sigma.est = 'Hajekplus'
                                 , draw.u = 'hotdeck'
                                 , informative = FALSE
                                 , runKY = TRUE
                                 , fpcKY = FALSE
                                 , approxKY = FALSE
                                 , runTradMI = FALSE
                                 , writeExcel = TRUE
                                 , nameExcel = 'sim2b_KY_nofpc_20230305'
                                 , seed = 20230305
)

# check: original scenario from Kim & Yang's paper
resMIc.KY <- runSimulation(S = 10000L
                           , S_at_a_time = 5000L
                           , N = 50000L
                           , cfs = 2L
                           , m = 100L
                           , model.x = 'rnorm(n = N, mean = 0, sd = 1)'
                           , model.u = 'rnorm(n = N, mean = 0, sd = 1)'
                           , model.eps = 'rnorm(n = N, mean = 0, sd = 1)'
                           , model.y = '-0.5 + 0.5 * pop.x + pop.eps'
                           , model.logit.pi = '-4 - 0.5 * pop.x'
                           , model.logit.phi = '1 + 0.5 * pop.x + 0.5 * pop.u'
                           , sigma.est = 'Hajekplus'
                           , draw.u = 'hotdeck'
                           , informative = FALSE
                           , runKY = TRUE
                           , fpcKY = TRUE
                           , approxKY = FALSE
                           , runTradMI = FALSE
                           , writeExcel = TRUE
                           , nameExcel = 'sim2c_KY_20230305'
                           , seed = 20230305
)
