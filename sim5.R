
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
                        , model.x = 'runif(n = N, min = 2.5, max = 5.5)'
                        , model.eps = '(exp(rnorm(n = N, mean = 3, sd = sqrt(0.22305))) - exp(3 + 0.22305/2)) / sqrt(126)'
                        , model.y = 'pop.x + pop.eps'
                        , model.pi = '200 * pop.x / (4 * N)'
                        , model.logit.phi = '9 - 2 * pop.x'
                        , draw.u = 'hotdeck'
                        , informative = FALSE
                        , runKY = TRUE
                        , fpcKY = TRUE
                        , approxKY = TRUE
                        , runTradMI = TRUE
                        , writeExcel = TRUE
                        , nameExcel = 'sim5a_HD_20230305'
                        , seed = 20230305
)

resMIa1 <- runSimulation(S = 100000L
                         , S_at_a_time = 10000L
                         , N = 10000L
                         , cfs = 1L
                         , use.intercept = FALSE
                         , m = 100L
                         , model.x = 'runif(n = N, min = 2.5, max = 5.5)'
                         , model.eps = '(exp(rnorm(n = N, mean = 3, sd = sqrt(0.22305))) - exp(3 + 0.22305/2)) / sqrt(126)'
                         , model.y = 'pop.x + pop.eps'
                         , model.pi = '200 * pop.x / (4 * N)'
                         , model.logit.phi = '9 - 2 * pop.x'
                         , sigma.est = 'Hajekplus'
                         , draw.u = 'hotdeck'
                         , informative = FALSE
                         , runKY = TRUE
                         , fpcKY = TRUE
                         , approxKY = TRUE
                         , runTradMI = TRUE
                         , writeExcel = TRUE
                         , nameExcel = 'sim5a_HD_no_icpt_20230305'
                         , seed = 20230305
)


resMIb <- runSimulation(S = 100000L
                        , S_at_a_time = 10000L
                        , N = 10000L
                        , cfs = 2L
                        , m = 100L
                        , model.x = 'runif(n = N, min = 2.5, max = 5.5)'
                        , model.eps = '(exp(rnorm(n = N, mean = 3, sd = sqrt(0.22305))) - exp(3 + 0.22305/2)) / sqrt(126)'
                        , model.y = 'pop.x + pop.eps'
                        , model.pi = '200 * pop.y[ ,s] / (4 * N)'
                        , model.logit.phi = '9 - 2 * pop.x'
                        , sigma.est = 'Hajekplus'
                        , draw.u = 'hotdeck'
                        , informative = TRUE
                        , runKY = TRUE
                        , fpcKY = TRUE
                        , approxKY = TRUE
                        , runTradMI = TRUE
                        , writeExcel = TRUE
                        , nameExcel = 'sim5b_HD_20230305'
                        , seed = 20230305
)

resMIb1 <- runSimulation(S = 100000L
                         , S_at_a_time = 10000L
                         , N = 10000L
                         , cfs = 1L
                         , use.intercept = FALSE
                         , m = 100L
                         , model.x = 'runif(n = N, min = 2.5, max = 5.5)'
                         , model.eps = '(exp(rnorm(n = N, mean = 3, sd = sqrt(0.22305))) - exp(3 + 0.22305/2)) / sqrt(126)'
                         , model.y = 'pop.x + pop.eps'
                         , model.pi = '200 * pop.y[ ,s] / (4 * N)'
                         , model.logit.phi = '9 - 2 * pop.x'
                         , sigma.est = 'Hajekplus'
                         , draw.u = 'hotdeck'
                         , informative = TRUE
                         , runKY = TRUE
                         , fpcKY = TRUE
                         , approxKY = TRUE
                         , runTradMI = TRUE
                         , writeExcel = TRUE
                         , nameExcel = 'sim5b_HD_no_icpt_20230305'
                         , seed = 20230305
)

##### Methode Kim & Yang

resMIa.KY <- runSimulation(S = 10000L
                           , S_at_a_time = 10000L
                           , N = 10000L
                           , cfs = 2L
                           , m = 100L
                           , model.x = 'runif(n = N, min = 2.5, max = 5.5)'
                           , model.eps = '(exp(rnorm(n = N, mean = 3, sd = sqrt(0.22305))) - exp(3 + 0.22305/2)) / sqrt(126)'
                           , model.y = 'pop.x + pop.eps'
                           , model.pi = '200 * pop.x / (4 * N)'
                           , model.logit.phi = '9 - 2 * pop.x'
                           , sigma.est = 'Hajekplus'
                           , informative = FALSE
                           , runKY = TRUE
                           , fpcKY = TRUE
                           , writeExcel = TRUE
                           , nameExcel = 'sim5a_KY_20230305'
                           , seed = 20230305
)

resMIa1.KY <- runSimulation(S = 10000L
                            , S_at_a_time = 10000L
                            , N = 10000L
                            , cfs = 1L
                            , use.intercept = FALSE
                            , m = 100L
                            , model.x = 'runif(n = N, min = 2.5, max = 5.5)'
                            , model.eps = '(exp(rnorm(n = N, mean = 3, sd = sqrt(0.22305))) - exp(3 + 0.22305/2)) / sqrt(126)'
                            , model.y = 'pop.x + pop.eps'
                            , model.pi = '200 * pop.x / (4 * N)'
                            , model.logit.phi = '9 - 2 * pop.x'
                            , sigma.est = 'Hajekplus'
                            , informative = FALSE
                            , runKY = TRUE
                            , fpcKY = TRUE
                            , writeExcel = TRUE
                            , nameExcel = 'sim5a_KY_no_icpt_20230305'
                            , seed = 20230305
)

resMIb.KY <- runSimulation(S = 10000L
                           , S_at_a_time = 10000L
                           , N = 10000L
                           , cfs = 2L
                           , m = 100L
                           , model.x = 'runif(n = N, min = 2.5, max = 5.5)'
                           , model.eps = '(exp(rnorm(n = N, mean = 3, sd = sqrt(0.22305))) - exp(3 + 0.22305/2)) / sqrt(126)'
                           , model.y = 'pop.x + pop.eps'
                           , model.pi = '200 * pop.y[ ,s] / (4 * N)'
                           , model.logit.phi = '9 - 2 * pop.x'
                           , sigma.est = 'Hajekplus'
                           , informative = TRUE
                           , runKY = TRUE
                           , fpcKY = TRUE
                           , writeExcel = TRUE
                           , nameExcel = 'sim5b_KY_20230305'
                           , seed = 20230305
)

resMIb1.KY <- runSimulation(S = 10000L
                            , S_at_a_time = 10000L
                            , N = 10000L
                            , cfs = 1L
                            , use.intercept = FALSE
                            , m = 100L
                            , model.x = 'runif(n = N, min = 2.5, max = 5.5)'
                            , model.eps = '(exp(rnorm(n = N, mean = 3, sd = sqrt(0.22305))) - exp(3 + 0.22305/2)) / sqrt(126)'
                            , model.y = 'pop.x + pop.eps'
                            , model.pi = '200 * pop.y[ ,s] / (4 * N)'
                            , model.logit.phi = '9 - 2 * pop.x'
                            , sigma.est = 'Hajekplus'
                            , informative = TRUE
                            , runKY = TRUE
                            , fpcKY = TRUE
                            , writeExcel = TRUE
                            , nameExcel = 'sim5b_KY_no_icpt_20230305'
                            , seed = 20230305
)
