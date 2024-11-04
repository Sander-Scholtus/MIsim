
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
                        , N = 10000L
                        , cfs = 2L
                        , m = 100L
                        , model.x = 'runif(n = N, min = 2, max = 5)'
                        , model.eps = 'runif(n = N, min = -1.5, max = 1.5)'
                        , model.y = 'pop.x + pop.eps'
                        , model.pi = '200 * pop.x / (3.5 * N)'
                        , model.logit.phi = '8 - 2 * pop.x'
                        , sigma.est = 'Hajekplus'
                        , draw.u = 'hotdeck'
                        , informative = FALSE
                        , runKY = TRUE
                        , fpcKY = TRUE
                        , approxKY = TRUE
                        , runTradMI = TRUE
                        , writeExcel = TRUE
                        , nameExcel = 'sim4a_HD_20230305'
                        , seed = 20230305
)

resMIa1 <- runSimulation(S = 100000L
                         , S_at_a_time = 10000L
                         , N = 10000L
                         , cfs = 1L
                         , use.intercept = FALSE
                         , m = 100L
                         , model.x = 'runif(n = N, min = 2, max = 5)'
                         , model.eps = 'runif(n = N, min = -1.5, max = 1.5)'
                         , model.y = 'pop.x + pop.eps'
                         , model.pi = '200 * pop.x / (3.5 * N)'
                         , model.logit.phi = '8 - 2 * pop.x'
                         , sigma.est = 'Hajekplus'
                         , draw.u = 'hotdeck'
                         , informative = FALSE
                         , runKY = TRUE
                         , fpcKY = TRUE
                         , approxKY = TRUE
                         , runTradMI = TRUE
                         , writeExcel = TRUE
                         , nameExcel = 'sim4a_HD_no_icpt_20230305'
                         , seed = 20230305
)


resMIb <- runSimulation(S = 100000L
                        , S_at_a_time = 10000L
                        , N = 10000L
                        , cfs = 2L
                        , m = 100L
                        , model.x = 'runif(n = N, min = 2, max = 5)'
                        , model.eps = 'runif(n = N, min = -1.5, max = 1.5)'
                        , model.y = 'pop.x + pop.eps'
                        , model.pi = '200 * pop.y[ ,s] / (3.5 * N)'
                        , model.logit.phi = '8 - 2 * pop.x'
                        , sigma.est = 'Hajekplus'
                        , draw.u = 'hotdeck'
                        , informative = TRUE
                        , runKY = TRUE
                        , fpcKY = TRUE
                        , approxKY = TRUE
                        , runTradMI = TRUE
                        , writeExcel = TRUE
                        , nameExcel = 'sim4b_HD_20230305'
                        , seed = 20230305
)

resMIb1 <- runSimulation(S = 100000L
                         , S_at_a_time = 10000L
                         , N = 10000L
                         , cfs = 1L
                         , use.intercept = FALSE
                         , m = 100L
                         , model.x = 'runif(n = N, min = 2, max = 5)'
                         , model.eps = 'runif(n = N, min = -1.5, max = 1.5)'
                         , model.y = 'pop.x + pop.eps'
                         , model.pi = '200 * pop.y[ ,s] / (3.5 * N)'
                         , model.logit.phi = '8 - 2 * pop.x'
                         , sigma.est = 'Hajekplus'
                         , draw.u = 'hotdeck'
                         , informative = TRUE
                         , runKY = TRUE
                         , fpcKY = TRUE
                         , approxKY = TRUE
                         , runTradMI = TRUE
                         , writeExcel = TRUE
                         , nameExcel = 'sim4b_HD_no_icpt_20230305'
                         , seed = 20230305
)


##### Methode Kim & Yang

resMIa.KY <- runSimulation(S = 10000L
                           , S_at_a_time = 10000L
                           , N = 10000L
                           , cfs = 2L
                           , m = 100L
                           , model.x = 'runif(n = N, min = 2, max = 5)'
                           , model.eps = 'runif(n = N, min = -1.5, max = 1.5)'
                           , model.y = 'pop.x + pop.eps'
                           , model.pi = '200 * pop.x / (3.5 * N)'
                           , model.logit.phi = '8 - 2 * pop.x'
                           , sigma.est = 'Hajekplus'
                           , draw.u = 'hotdeck'
                           , informative = FALSE
                           , runKY = TRUE
                           , fpcKY = TRUE
                           , approxKY = FALSE
                           , runTradMI = TRUE
                           , writeExcel = TRUE
                           , nameExcel = 'sim4a_KY_20230305'
                           , seed = 20230305
)

resMIa1.KY <- runSimulation(S = 10000L
                            , S_at_a_time = 10000L
                            , N = 10000L
                            , cfs = 1L
                            , use.intercept = FALSE
                            , m = 100L
                            , model.x = 'runif(n = N, min = 2, max = 5)'
                            , model.eps = 'runif(n = N, min = -1.5, max = 1.5)'
                            , model.y = 'pop.x + pop.eps'
                            , model.pi = '200 * pop.x / (3.5 * N)'
                            , model.logit.phi = '8 - 2 * pop.x'
                            , sigma.est = 'Hajekplus'
                            , draw.u = 'hotdeck'
                            , informative = FALSE
                            , runKY = TRUE
                            , fpcKY = TRUE
                            , approxKY = FALSE
                            , runTradMI = TRUE
                            , writeExcel = TRUE
                            , nameExcel = 'sim4a_KY_no_icpt_20230305'
                            , seed = 20230305
)


resMIb.KY <- runSimulation(S = 10000L
                           , S_at_a_time = 10000L
                           , N = 10000L
                           , cfs = 2L
                           , m = 100L
                           , model.x = 'runif(n = N, min = 2, max = 5)'
                           , model.eps = 'runif(n = N, min = -1.5, max = 1.5)'
                           , model.y = 'pop.x + pop.eps'
                           , model.pi = '200 * pop.y[ ,s] / (3.5 * N)'
                           , model.logit.phi = '8 - 2 * pop.x'
                           , sigma.est = 'Hajekplus'
                           , draw.u = 'hotdeck'
                           , informative = TRUE
                           , runKY = TRUE
                           , fpcKY = TRUE
                           , approxKY = FALSE
                           , runTradMI = TRUE
                           , writeExcel = TRUE
                           , nameExcel = 'sim4b_KY_20230305'
                           , seed = 20230305
)

resMIb1.KY <- runSimulation(S = 10000L
                            , S_at_a_time = 10000L
                            , N = 10000L
                            , cfs = 1L
                            , use.intercept = FALSE
                            , m = 100L
                            , model.x = 'runif(n = N, min = 2, max = 5)'
                            , model.eps = 'runif(n = N, min = -1.5, max = 1.5)'
                            , model.y = 'pop.x + pop.eps'
                            , model.pi = '200 * pop.y[ ,s] / (3.5 * N)'
                            , model.logit.phi = '8 - 2 * pop.x'
                            , sigma.est = 'Hajekplus'
                            , draw.u = 'hotdeck'
                            , informative = TRUE
                            , runKY = TRUE
                            , fpcKY = TRUE
                            , approxKY = FALSE
                            , runTradMI = TRUE
                            , writeExcel = TRUE
                            , nameExcel = 'sim4b_KY_no_icpt_20230305'
                            , seed = 20230305
)

#####

# check: previously tested scenario with KY from Knottnerus' paper

resMIc.KY <- runSimulation(S = 10000L
                           , S_at_a_time = 10000L
                           , N = 10000L
                           , cfs = 2L
                           , m = 100L
                           , model.x = 'runif(n = N, min = 0.5, max = 3.5)'
                           , model.eps = 'rnorm(n = N, mean = 0, sd = 1)'
                           , model.y = 'pop.x + pop.eps'
                           , model.pi = '100 * pop.x / N'
                           , model.logit.phi = '5 - 2 * pop.x'
                           , sigma.est = 'Hajekplus'
                           , draw.u = 'hotdeck'
                           , informative = FALSE
                           , runKY = TRUE
                           , fpcKY = TRUE
                           , approxKY = FALSE
                           , runTradMI = TRUE
                           , writeExcel = TRUE
                           , nameExcel = 'sim4c_KY_20230305'
                           , seed = 20230305
)

resMIc1.KY <- runSimulation(S = 10000L
                            , S_at_a_time = 10000L
                            , N = 10000L
                            , cfs = 1L
                            , use.intercept = FALSE
                            , m = 100L
                            , model.x = 'runif(n = N, min = 0.5, max = 3.5)'
                            , model.eps = 'rnorm(n = N, mean = 0, sd = 1)'
                            , model.y = 'pop.x + pop.eps'
                            , model.pi = '100 * pop.x / N'
                            , model.logit.phi = '5 - 2 * pop.x'
                            , sigma.est = 'Hajekplus'
                            , draw.u = 'hotdeck'
                            , informative = FALSE
                            , runKY = TRUE
                            , fpcKY = TRUE
                            , approxKY = FALSE
                            , runTradMI = TRUE
                            , writeExcel = TRUE
                            , nameExcel = 'sim4c_KY_no_icpt_20230305'
                            , seed = 20230305
)
