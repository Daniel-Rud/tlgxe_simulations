
# Libraries #########################################
# this line is for BART machine 
options(java.parameters = "-Xmx5g")
library(ggplot2)
library(magrittr)
library(future)
library(future.apply)
library(bigmemory)
library(ggpubr)
library(mvtnorm)
library(lmtest)
library(progressr)
library(speedglm)
library(nleqslv)
library(bartMachine)
library(xgboost)
library(dplyr) # for pivot longer
library(tidyr)
source("/Users/danielrud/Desktop/USC/Targeted Learning/tlgxe/tmle_gxe_function.R")
library(ipw) # for inverse prob trt wts
#####################################################



# bernoulli random observations 
rbern = function(p) return(rbinom(length(p), 1, p))

# logit function
logit = function(x) return(log(x / (1 - x)))

# function to perform G computation 
gcomp = function(outcome_data, family = "binomial",
                 formula = NULL, 
                 pred_ids = NULL)
{
  
  if(is.null(formula))
  {
    formula = Y ~ . 
  }
  
  Q0_mod = glm(formula = formula, data = outcome_data, family = family)
  
  gcomp_copy_1 = outcome_data
  gcomp_copy_0 = outcome_data
  
  gcomp_copy_1$A = 1
  gcomp_copy_0$A = 0
  
  if(is.null(pred_ids))
  {
    pred_ids = 1:nrow(outcome_data)
  }
  
  suppressWarnings(Q0AW <- predict(Q0_mod, newdata = outcome_data[pred_ids, ], type = "response"))
  suppressWarnings(Q1W <- predict(Q0_mod, newdata = gcomp_copy_1[pred_ids, ], type = "response"))
  suppressWarnings(Q0W <- predict(Q0_mod, newdata = gcomp_copy_0[pred_ids, ], type = "response"))
  
  return(list(Q0AW  = Q0AW ,Q1W = Q1W, 
              Q0W = Q0W ))
}

# function to perform G computation given a prefit glm on the exposure names `A`
gcomp_glm_mod = function(outcome_data, Q0_mod)
{
  gcomp_copy_1 = outcome_data
  gcomp_copy_0 = outcome_data
  
  gcomp_copy_1$A = 1
  gcomp_copy_0$A = 0
  
  Q0AW = predict(Q0_mod, newdata = outcome_data, type = "response")
  Q1W = predict(Q0_mod, newdata = gcomp_copy_1, type = "response")
  Q0W = predict(Q0_mod, newdata = gcomp_copy_0, type = "response")
  
  return(list(Q0AW  = Q0AW ,Q1W = Q1W, 
              Q0W = Q0W ))
}

# output ACE effects from g computation results 
produce_ACE_from_gcomp = function(g_comp, ACE_type)
{
  EY1 = mean(g_comp$Q1W)
  EY0 = mean(g_comp$Q0W)
  ACE = 0 
  if(ACE_type == "ATE" )
  {
    ACE = EY1 - EY0
  }else
  {
    ACE = (EY1 / (1-EY1)) / (EY0 / (1-EY0))
  }
  return(ACE)
}

# used to find A effects for generating data according to prespecified ATE levels
ATE_func_uniroot = function(x, S_0, S_1, ATE, family)
{
  f = NULL
  if(family == "binomial")
  {
    f = mean( plogis(x + S_1) - plogis(S_0)) - ATE 
  }else
  {
    f = mean(x + S_1 - S_0) - ATE
  }
  
  return(f)
}

# used to find A effects for generating data according to prespecified MOR levels
MOR_func_uniroot = function(x, S_0, S_1, MOR, family)
{
  m1 = NULL 
  m0 = NULL
  
  if(family == "binomial")
  {
    m1 = mean(plogis(x + S_1))
    m0 = mean(plogis(S_0))
  }else
  {
    m1 = mean(x + S_1)
    m0 = mean(S_0)
  }
  
  f = (m1 / (1-m1)) / (m0 / (1-m0)) - MOR
  
  return(f)
  
}


# generate the intercept effect for propensity model to obey exposure prevalence
# and the A effect at each level according to an ACE 
gen_A_prop_int_effect =  function(n = 1000000, n_snp = 10, 
                                  family = "binomial", 
                                  snps_propensity = c(2:5), # which SNPs other than the first SNP to consider in propensity 
                                  snp_propensity_effects = c(.5,.1,-.2,-.1), 
                                  exposure_prev = .5,  # this will control the intercept 
                                  SNP_MAF = rep(.5, n_snp),
                                  SNP1_intercept_effs = c(-1,-1,-1), 
                                  SNP_adj_effs_0 = rnorm(n_snp - 1, sd = .1), 
                                  SNP_adj_effs_1 = rnorm(n_snp - 1, sd = .1),
                                  SNP_adj_effs_2 = rnorm(n_snp - 1, sd = .1),
                                  ACE_type = c("ATE", "MOR"), 
                                  ACE = c(.2,.2,.2),
                                  rho = -1, 
                                  confounder_data_propensity_formula, 
                                  confounder_propensity_effects, 
                                  confounder_data_outcome_formula, 
                                  confounder_outcome_effects)
{
  
  ### Define correlation structure ############################################
  corr_mat = 0 
  # if rho is set to -1, generate independent SNPs
  if(rho == -1)
  {
    corr_mat = diag(n_snp)
  }else # if rho != - 1, generate according to proximity based decay.  
  {
    corr_mat = gen_corr_mat(n_snp, corr_factor = rho) 
  }
  ############################################################################
  
  ### Generate SNP data ######################################################
  snp_data = rmvnorm(n = n, mean = rep(0,n_snp), sigma = corr_mat)
  snp_data = sapply(1:n_snp,
                    FUN = function(x)
                    {
                      
                      col = snp_data[,x]
                      # SNPs discretized according to binomial probabilities from MAF
                      quantiles = quantile(col, probs = c((1-SNP_MAF[x])^2,
                                                          (1-SNP_MAF[x])^2 + 2*(SNP_MAF[x])*(1-SNP_MAF[x]))) # p^2 + 2pq
                      SNP = numeric(length(col))
                      
                      SNP_0_ind  = which(col <= quantiles[1])
                      SNP_1_ind = which(col <= quantiles[2] & col > quantiles[1])
                      SNP_2_ind = which(col > quantiles[2])
                      
                      SNP = col
                      SNP[SNP_0_ind] = 0
                      SNP[SNP_1_ind] = 1
                      SNP[SNP_2_ind] = 2
                      return(SNP)
                    })
  
  ###########################################################################
  # GENERATE CONFOUNDING DATA WITH NONLINEAR CONFOUNDING
  ###########################################################################
  age = rnorm(n, mean = 50, sd = 5)
  sex = sample(c(0,1), size = n, replace = T)
  cohort = sample(1:5, size = n, replace = T) %>% factor 
  ancestry_1 = rnorm(n, mean = 0, sd = 1)
  ancestry_2 = rnorm(n, mean = 0, sd = 1)
  
  confounder_data = data.frame(age = age,
                               sex = sex, 
                               cohort = cohort, 
                               ancestry_1 = ancestry_1, 
                               ancestry_2 = ancestry_2, 
                               X1 = snp_data[,1])
  
  confounder_data_propensity = model.matrix(confounder_data_propensity_formula, 
                                            data = confounder_data)[,-1]

  
  
  ### Generate exposure data#################################################
  
  exposure_data = numeric(n)
  propensity_intercept = 0
  
  if(!is.na(snps_propensity[1])) # if we supply 
  {
    # find intercept for propensity model so that we retrieve exposure prevalence
    propensity_intercept = uniroot(f = function(int, p, XB) { mean(plogis(XB+int)) - p}, 
                                   interval = c(-10, 10), 
                                   p = exposure_prev, 
                                   XB = tcrossprod(snp_data[, snps_propensity], matrix(snp_propensity_effects, nrow = 1)) + 
                                     tcrossprod(confounder_data_propensity, matrix(confounder_propensity_effects, nrow = 1)))$root
    
    exposure_data = (tcrossprod(snp_data[, snps_propensity], matrix(snp_propensity_effects, nrow = 1)) + 
                       tcrossprod(confounder_data_propensity, matrix(confounder_propensity_effects, nrow = 1)) +  propensity_intercept) %>% plogis %>% rbern
  }else
  {
    # find intercept for propensity model so that we retrieve exposure prevalence
    propensity_intercept = uniroot(f = function(int, p, XB) { mean(plogis(XB+int)) - p}, 
                                   interval = c(-10, 10), 
                                   p = exposure_prev, 
                                   XB = tcrossprod(confounder_data_propensity, matrix(confounder_propensity_effects, nrow = 1)))$root
    
    exposure_data = (tcrossprod(confounder_data_propensity, matrix(confounder_propensity_effects, nrow = 1)) +  
                       propensity_intercept) %>% plogis %>% rbern
  }
  
  # generate outcome model model matrix 
  confounder_data_outcome = model.matrix(confounder_data_outcome_formula, 
                                         data = cbind(A = exposure_data, confounder_data))[,-1]
  confounder_data_outcome_0 = model.matrix(confounder_data_outcome_formula, 
                                         data = cbind(A = 0, confounder_data))[,-1]
  confounder_data_outcome_1 = model.matrix(confounder_data_outcome_formula, 
                                           data = cbind(A = 1, confounder_data))[,-1]
  
  ###########################################################################
  
  total_data = data.frame(A = exposure_data, snp_data, confounder_data_outcome)
  total_data_0 = data.frame(A = 0, snp_data, confounder_data_outcome_0)
  total_data_1 = data.frame(A = 1, snp_data, confounder_data_outcome_1)
  
  data_subsets = list("X1_0" = total_data[which(total_data[,2] == 0), ],
                      "X1_1" = total_data[which(total_data[,2] == 1), ], 
                      "X1_2" = total_data[which(total_data[,2] == 2), ])
  
  data_subsets_0 = list("X1_0" = total_data_0[which(total_data_0[,2] == 0), ],
                      "X1_1" = total_data_0[which(total_data_0[,2] == 1), ], 
                      "X1_2" = total_data_0[which(total_data_0[,2] == 2), ])
  
  data_subsets_1 = list("X1_0" = total_data_1[which(total_data_1[,2] == 0), ],
                      "X1_1" = total_data_1[which(total_data_1[,2] == 1), ], 
                      "X1_2" = total_data_1[which(total_data_1[,2] == 2), ])
  
  SNP_adj_effects = rbind(SNP_adj_effs_0, SNP_adj_effs_1, SNP_adj_effs_2)
  
  Y = numeric(length = n)
  
  ind = 1 
  A_effects = numeric(3)
  for(i in 1:3)
  {
      current_data_0 = data_subsets_0[[i]] %>% as.matrix
      current_data_1 = data_subsets_1[[i]] %>% as.matrix
      
      # remember -- these wont always be the same if there are A*confounder terms 
      X_other_0 = current_data_0[, -c(1:2)]
      X_other_1 = current_data_1[, -c(1:2)]
      
      S_0 =  X_other_0 %*% matrix(c(SNP_adj_effects[i,], confounder_outcome_effects), ncol = 1) + 
        SNP1_intercept_effs[i]
      S_1 = X_other_1 %*% matrix(c(SNP_adj_effects[i,], confounder_outcome_effects), ncol = 1) + 
        SNP1_intercept_effs[i]
      
      
      A_effect = ifelse(ACE_type == "ATE",
                        uniroot(ATE_func_uniroot, S_0 = S_0, S_1 = S_1, ATE = ACE[i], family = family, 
                                interval = c(-10,10), tol = 1E-18, 
                                extendInt = "yes")[[1]],
                        uniroot(MOR_func_uniroot, S_0 = S_0, S_1 = S_1, MOR = ACE[i], family = family, 
                                interval = c(-10,10), tol = 1E-18, 
                                extendInt = "yes")[[1]])
      
      
      A_effects[i] = A_effect
      
      if(family == "binomial")
      {
        Y1 = (A_effect + S_1) %>% plogis %>% rbern
        Y0 = (S_0) %>% plogis %>% rbern
        A = data_subsets[[i]][,1]
        Y[ind:(ind +nrow(current_data_0) - 1)] =  A*Y1 + (1-A)*Y0 
      }else
      {
        Y1 = (A_effect + S_1) 
        Y0 = (S_0) 
        A = data_subsets[[i]][,1]
        Y[ind:(ind +nrow(current_data_0) - 1)] =  A*Y1 + (1-A)*Y0 
      }
      
      ind = ind +nrow(current_data_0)
  }
  
  
  # test with gcomp, output results of gcomp 
  
  # make a dataframe with the data with confounders in original form -- 
  # need to be careful merging with Y!!!!! 
  
  data = data.frame(A = exposure_data, snp_data, 
                    confounder_data[,-which(colnames(confounder_data) == "X1")])
  
  data = list(data[which(data$X1 == 0), ], 
              data[which(data$X1 == 1), ], 
              data[which(data$X1 == 2), ])
  
  data = data.frame(Y = Y, do.call(rbind, data) )

  colnames(data)[3:(3+n_snp - 1)] = paste0("X",1:n_snp )
  
  data_subset_0 = data[which(data$X1 == 0), -3]
  data_subset_1 = data[which(data$X1 == 1), -3]
  data_subset_2 = data[which(data$X1 == 2), -3]
  
  formula = update(confounder_data_outcome_formula, 
                formula(paste0("Y ~ A + ", paste0("X", 2:n_snp, collapse = "+"), "+ .")))
  
  
  gcomp_0 = gcomp(outcome_data = data_subset_0, family = family, 
                  formula = formula)
  gcomp_1 = gcomp(outcome_data = data_subset_1, family = family,
                  formula = formula)
  gcomp_2 = gcomp(outcome_data = data_subset_2, family = family, 
                  formula = formula)
  
  est_ATE_0 = mean(gcomp_0$Q1W - gcomp_0$Q0W)
  est_ATE_1 = mean(gcomp_1$Q1W - gcomp_1$Q0W)
  est_ATE_2 = mean(gcomp_2$Q1W - gcomp_2$Q0W)
  
  est_MOR_0 = (mean(gcomp_0$Q1W) / (1 - mean(gcomp_0$Q1W))) /
    (mean(gcomp_0$Q0W) / (1 - mean(gcomp_0$Q0W)))
  est_MOR_1 = (mean(gcomp_1$Q1W) / (1 - mean(gcomp_1$Q1W))) /
    (mean(gcomp_1$Q0W) / (1 - mean(gcomp_1$Q0W)))
  est_MOR_2 = (mean(gcomp_2$Q1W) / (1 - mean(gcomp_2$Q1W))) /
    (mean(gcomp_2$Q0W) / (1 - mean(gcomp_2$Q0W)))
  
  
  effects = c(ATE_0 = est_ATE_0, ATE_1 = est_ATE_1, ATE_2 = est_ATE_2,
              MOR_0 = est_MOR_0, MOR_1 = est_MOR_1, MOR_2 = est_MOR_2)
  
  
  return(list(propensity_intercept = propensity_intercept, 
              A_effects = A_effects, 
              sim_ACEs = effects, 
              disease_prev = mean(Y) ))
  
}

# generate data given the A effects found from `gen_A_prop_int_effect`
gen_data_ACE = function(n = 1000, n_snp = 10, 
                        family = "binomial",
                        snps_propensity = c(2:5), # which SNPs other than the first SNP to consider in propensity 
                        snp_propensity_effects = c(.5,.1,-.2,-.1), 
                        propensity_intercept = 0, 
                        SNP_MAF = rep(.5, n_snp),
                        SNP1_intercept_effs = c(-1,-1,-1), 
                        SNP_adj_effs_0 = rnorm(n_snp - 1, sd = .1), 
                        SNP_adj_effs_1 = rnorm(n_snp - 1, sd = .1),
                        SNP_adj_effs_2 = rnorm(n_snp - 1, sd = .1),
                        ACE_type = c("ATE", "MOR"), 
                        A_effects = rep(0, 3),
                        rho = -1, 
                        confounder_data_propensity_formula, 
                        confounder_propensity_effects, 
                        confounder_data_outcome_formula, 
                        confounder_outcome_effects)
{
  
  ### Define correlation structure ############################################
  corr_mat = 0 
  # if rho is set to -1, generate independent SNPs
  if(rho == -1)
  {
    corr_mat = diag(n_snp)
  }else # if rho != - 1, generate according to proximity based decay.  
  {
    corr_mat = gen_corr_mat(n_snp, corr_factor = rho) 
  }
  ############################################################################
  
  ### Generate SNP data ######################################################
  snp_data = rmvnorm(n = n, mean = rep(0,n_snp), sigma = corr_mat)
  snp_data = sapply(1:n_snp,
                    FUN = function(x)
                    {
                      
                      col = snp_data[,x]
                      # SNPs discretized according to binomial probabilities from MAF
                      quantiles = quantile(col, probs = c((1-SNP_MAF[x])^2,
                                                          (1-SNP_MAF[x])^2 + 2*(SNP_MAF[x])*(1-SNP_MAF[x]))) # p^2 + 2pq
                      SNP = numeric(length(col))
                      
                      
                      
                      SNP_0_ind  = which(col <= quantiles[1])
                      SNP_1_ind = which(col <= quantiles[2] & col > quantiles[1])
                      SNP_2_ind = which(col > quantiles[2])
                      
                      SNP = col
                      SNP[SNP_0_ind] = 0
                      SNP[SNP_1_ind] = 1
                      SNP[SNP_2_ind] = 2
                      return(SNP)
                    })
  
  ###########################################################################
  
  ###########################################################################
  # GENERATE CONFOUNDING DATA WITH NONLINEAR CONFOUNDING
  ###########################################################################
  age = rnorm(n, mean = 50, sd = 5)
  sex = sample(c(0,1), size = n, replace = T)
  cohort = sample(1:5, size = n, replace = T) %>% factor 
  ancestry_1 = rnorm(n, mean = 0, sd = 1)
  ancestry_2 = rnorm(n, mean = 0, sd = 1)
  
  confounder_data_orig = data.frame(age = age,
                                    sex = sex, 
                                    cohort = cohort, 
                                    ancestry_1 = ancestry_1, 
                                    ancestry_2 = ancestry_2)
  
  confounder_data_propensity = model.matrix(confounder_data_propensity_formula, 
                                            data = cbind(confounder_data_orig, 
                                                         'X1' = snp_data[,1]))[,-1]
  
  
  ### Generate exposure data#################################################
  
  exposure_data = numeric(n_snp)
  
  if(!is.na(snps_propensity[1])) # if we supply 
  {
    exposure_data = (tcrossprod(snp_data[, snps_propensity], matrix(snp_propensity_effects, nrow = 1)) + 
                       tcrossprod(confounder_data_propensity, matrix(confounder_propensity_effects, nrow = 1)) +  propensity_intercept) %>% 
      plogis %>% rbern
  }else
  {
    exposure_data = (tcrossprod(confounder_data_propensity, matrix(confounder_propensity_effects, nrow = 1)) +  
                       propensity_intercept) %>% plogis %>% rbern
  }
  
  confounder_data_outcome = model.matrix(confounder_data_outcome_formula, 
                                         data = cbind(A = exposure_data , confounder_data_orig))[,-1]
  confounder_data_outcome_0 = model.matrix(confounder_data_outcome_formula, 
                                           data = cbind(A = 0, confounder_data_orig))[,-1]
  confounder_data_outcome_1 = model.matrix(confounder_data_outcome_formula, 
                                           data = cbind(A = 1, confounder_data_orig))[,-1]
  
  ###########################################################################
  
  total_data = data.frame(A = exposure_data, snp_data, confounder_data_outcome)
  total_data_0 = data.frame(A = 0, snp_data, confounder_data_outcome_0)
  total_data_1 = data.frame(A = 1, snp_data, confounder_data_outcome_1)
  
  # this is the original confounder data 
  total_data_orig = data.frame(A = exposure_data, snp_data, confounder_data_orig)
  
  data_subsets = list("X1_0" = total_data[which(total_data[,2] == 0), ],
                      "X1_1" = total_data[which(total_data[,2] == 1), ], 
                      "X1_2" = total_data[which(total_data[,2] == 2), ])
  
  data_subsets_0 = list("X1_0" = total_data_0[which(total_data_0[,2] == 0), ],
                        "X1_1" = total_data_0[which(total_data_0[,2] == 1), ], 
                        "X1_2" = total_data_0[which(total_data_0[,2] == 2), ])
  
  data_subsets_1 = list("X1_0" = total_data_1[which(total_data_1[,2] == 0), ],
                        "X1_1" = total_data_1[which(total_data_1[,2] == 1), ], 
                        "X1_2" = total_data_1[which(total_data_1[,2] == 2), ])
  
  data_subsets_orig = list("X1_0" = total_data_orig[which(total_data_orig[,2] == 0), ],
                           "X1_1" = total_data_orig[which(total_data_orig[,2] == 1), ], 
                           "X1_2" = total_data_orig[which(total_data_orig[,2] == 2), ])
  
  SNP_adj_effects = rbind(SNP_adj_effs_0, SNP_adj_effs_1, SNP_adj_effs_2)
  
  Y = numeric(length = n)
  
  ind = 1
  for(i in 1:3)
  {
    current_data_0 = data_subsets_0[[i]] %>% as.matrix
    current_data_1 = data_subsets_1[[i]] %>% as.matrix
    
    # remember -- these wont always be the same if there are A*confounder terms 
    X_other_0 = current_data_0[, -c(1:2)]
    X_other_1 = current_data_1[, -c(1:2)]
    
    S_0 =  X_other_0 %*% matrix(c(SNP_adj_effects[i,], confounder_outcome_effects), ncol = 1) + 
      SNP1_intercept_effs[i]
    S_1 = X_other_1 %*% matrix(c(SNP_adj_effects[i,], confounder_outcome_effects), ncol = 1) + 
      SNP1_intercept_effs[i]
    
    A_effect = A_effects[i]
    
    if(family == "binomial")
    {
      Y1 = (A_effect + S_1) %>% plogis %>% rbern
      Y0 = (S_0) %>% plogis %>% rbern
      A = data_subsets[[i]][,1]
      Y[ind:(ind +nrow(current_data_0) - 1)] =  A*Y1 + (1-A)*Y0 
    }else
    {
      Y1 = (A_effect + S_1) 
      Y0 = (S_0) 
      A = data_subsets[[i]][,1]
      Y[ind:(ind +nrow(current_data_0) - 1)] =  A*Y1 + (1-A)*Y0 + 
        rnorm(length(Y0),mean = 0, sd = 2.5)
    }
    
    ind = ind +nrow(current_data_0)
  }
  
  data = data.frame(A = exposure_data, snp_data, confounder_data_orig)
  
  data = list(data[which(data$X1 == 0), ], 
              data[which(data$X1 == 1), ], 
              data[which(data$X1 == 2), ])
  
  data = data.frame(Y = Y, do.call(rbind, data) )
  
  colnames(data)[3:(3+n_snp - 1)] = paste0("X",1:n_snp)
  
  # 
  # 
  # # Code for testing with gcomp 
  # data_subset_0 = data[which(data$X1 == 0), -3]
  # data_subset_1 = data[which(data$X1 == 1), -3]
  # data_subset_2 = data[which(data$X1 == 2), -3]
  # 
  # formula = update(confounder_data_outcome_formula, 
  #                  formula(paste0("Y ~ A + ", paste0("X", 2:n_snp, collapse = "+"), "+ .")))
  # 
  # 
  # gcomp_0 = gcomp(outcome_data = data_subset_0, family = family, 
  #                 formula = formula)
  # gcomp_1 = gcomp(outcome_data = data_subset_1, family = family,
  #                 formula = formula)
  # gcomp_2 = gcomp(outcome_data = data_subset_2, family = family, 
  #                 formula = formula)
  # 
  # est_ATE_0 = mean(gcomp_0$Q1W - gcomp_0$Q0W)
  # est_ATE_1 = mean(gcomp_1$Q1W - gcomp_1$Q0W)
  # est_ATE_2 = mean(gcomp_2$Q1W - gcomp_2$Q0W)
  # 
  # est_MOR_0 = (mean(gcomp_0$Q1W) / (1 - mean(gcomp_0$Q1W))) /
  #   (mean(gcomp_0$Q0W) / (1 - mean(gcomp_0$Q0W)))
  # est_MOR_1 = (mean(gcomp_1$Q1W) / (1 - mean(gcomp_1$Q1W))) /
  #   (mean(gcomp_1$Q0W) / (1 - mean(gcomp_1$Q0W)))
  # est_MOR_2 = (mean(gcomp_2$Q1W) / (1 - mean(gcomp_2$Q1W))) /
  #   (mean(gcomp_2$Q0W) / (1 - mean(gcomp_2$Q0W)))
  # 
  # 
  # effects = c(ATE_0 = est_ATE_0, ATE_1 = est_ATE_1, ATE_2 = est_ATE_2,
  #             MOR_0 = est_MOR_0, MOR_1 = est_MOR_1, MOR_2 = est_MOR_2)
  # effects
  
  return(data)
}


tmle_comparison_1snp = function(num_sims = 1000, 
                                n_smp = 5000, 
                                n_snp = 10,
                                family = "binomial",
                                snps_propensity = c(2:5), 
                                snp_propensity_effects = c(.5,.1,-.2,-.1),
                                exposure_prev = .5, 
                                SNP_MAF = rep(.5, n_snp),
                                SNP1_intercept_effs = c(-1,-1,-1), 
                                SNP_adj_effs_0 = rnorm(n_snp - 1, sd = .1), 
                                SNP_adj_effs_1 = rnorm(n_snp - 1, sd = .1),
                                SNP_adj_effs_2 = rnorm(n_snp - 1, sd = .1),
                                ACE_type = c("ATE", "MOR"),
                                ACE = c(.2,.2,.2),
                                rho = -1, 
                                W_exposure_null = F, 
                                W_outcome_null = F, 
                                TMLE_args_list = list(
                                  outcome_method = c("glmnet_int", "glmnet", "gesso", "logistf"), 
                                  npv_thresh = (5/sqrt(length(Y)))/log(length(Y)), 
                                  near_positivity_method = c("trim", "rebound"), 
                                  nfolds_cv_Q_init = 10, 
                                  nfolds_cv_glmnet_propensity = 3, 
                                  nfolds_cv_glmnet_outcome = 3,
                                  alpha_outcome = .5, 
                                  alpha_propensity = .5, 
                                  clever_cov_propensity_wt = T
                                ), 
                                confounder_data_propensity_formula, 
                                confounder_propensity_effects, 
                                confounder_data_outcome_formula, 
                                confounder_outcome_effects, 
                                propensity_SL.library = c("SL.gam", 
                                                          "SL.glmnet", 
                                                          "SL.ranger", 
                                                          "SL.xgboost"), 
                                propensity_formula = NULL, 
                                future.seed = T)
# if W_exposure_null / W_outcome_null is true, fit propoensity / outcome model with no 
# confounders
{
  
  p = progressor(steps = num_sims)
  
  # J is used in W exposure and W outcome 
  j = 1 
  
  # need to set seed for generating A_prop_int_effects
  set.seed(future.seed)
  # generate A effects and propensity intercept for data generation 
  A_prop_int_effects = gen_A_prop_int_effect(n = 1000000, 
                                             n_snp = n_snp,
                                             family= family,
                                             snps_propensity = snps_propensity, 
                                             snp_propensity_effects = snp_propensity_effects, 
                                             exposure_prev = exposure_prev, 
                                             SNP_MAF = SNP_MAF, 
                                             SNP1_intercept_effs = SNP1_intercept_effs, 
                                             SNP_adj_effs_0 =  SNP_adj_effs_0, 
                                             SNP_adj_effs_1 = SNP_adj_effs_1, 
                                             SNP_adj_effs_2  = SNP_adj_effs_2 , 
                                             ACE_type = ACE_type, 
                                             ACE = ACE,
                                             rho = rho, 
                                             confounder_data_propensity_formula = confounder_data_propensity_formula, 
                                             confounder_propensity_effects = confounder_propensity_effects, 
                                             confounder_data_outcome_formula = confounder_data_outcome_formula, 
                                             confounder_outcome_effects = confounder_outcome_effects)
  {  # this bracket is used to run the parallelization commands and the future command all together
  plan(multisession)
  EM = future_lapply(1:num_sims, FUN = function(i) 
  {
    data = gen_data_ACE(n = n_smp, 
                        n_snp = n_snp, 
                        family = family, 
                        snps_propensity = snps_propensity, 
                        snp_propensity_effects = snp_propensity_effects, 
                        propensity_intercept = A_prop_int_effects$propensity_intercept,
                        SNP_MAF = SNP_MAF, 
                        SNP1_intercept_effs = SNP1_intercept_effs, 
                        SNP_adj_effs_0 =  SNP_adj_effs_0, 
                        SNP_adj_effs_1 = SNP_adj_effs_1, 
                        SNP_adj_effs_2  = SNP_adj_effs_2 , 
                        ACE_type = ACE_type, 
                        A_effects =  A_prop_int_effects$A_effects,
                        rho = rho, 
                        confounder_data_propensity_formula = confounder_data_propensity_formula, 
                        confounder_propensity_effects = confounder_propensity_effects, 
                        confounder_data_outcome_formula = confounder_data_outcome_formula, 
                        confounder_outcome_effects = confounder_outcome_effects)
    
    ############################################################################
    # Diagnostic plots for simulated data ######################################
    # table(data$Y)
    # plot(density(predict(glm(A ~ . - Y, data = data, family = "binomial"),
    #                      type = "response")))
    # summary(predict(glm(A ~ . - Y, data = data, family = "binomial"), 
    #                 type = "response"))
    ############################################################################
    
    # tests for power 
      # oracle model ################################################################
        # Create the string for the new terms X1, X2, ..., Xn
        new_terms <- paste0("X", 1:n_snp, collapse = " + ")
        # Combine the existing formula, interaction term, and new terms into a complete formula string
        formula_string <- paste("Y ~ A*X1 + ", new_terms, "+ .")
        # Convert the string to a formula and update the existing formula
        oracle_formula <- update(confounder_data_outcome_formula, as.formula(formula_string))
        oracle_mod = glm(oracle_formula, data = data, family = family)
        
      # oracle model w/o SNPS #######################################################
        
        oracle_formula_no_snp <- update(confounder_data_outcome_formula, Y ~ A*X1 + .)
        oracle_mod_no_snp = glm(oracle_formula_no_snp, data = data, family = family)
      
      # gxe 1df ################################################################
        glm_formula_std = formula(paste0("Y ~ A*X", j, "+  age +sex + cohort + ancestry_1 + ancestry_2"))
        standard_GxE_mod = glm(glm_formula_std, data = data, family = family,)
        
        glm_1df_power = ifelse(summary(standard_GxE_mod)$coefficients["A:X1",4] <= 0.05, 1, 0)
        
        
      # gxe 2df ################################################################
        glm_formula_no_int = formula(paste0("Y ~ A + factor(X", j,") +  age +sex + cohort + ancestry_1 + ancestry_2" ))
        glm_formula_int = formula(paste0("Y ~  A* factor(X", j,") +  age +sex + cohort + ancestry_1 + ancestry_2"))
        
        glm_no_int = glm(glm_formula_no_int, data = data, family = family)
        glm_int = glm(glm_formula_int, data = data, family = family)
        
        lr_test = lrtest(glm_no_int, glm_int)
        
        
      # tmle  ##################################################################
        # to turn cohort into categorical for tmle
        tmle_data = model.matrix( ~ . , data)[,-1] %>% data.frame()
        
        tmle_mod = tmle_gxe(Y = tmle_data$Y,
                            A = tmle_data$A, 
                            G = tmle_data[,3:12], 
                            W = tmle_data[, -c(1:2, 3:12)], 
                            family = family,
                            case_control_design = F, 
                            disease_prevalence = NULL,
                            weights = NULL, 
                            propensity_SL.library = propensity_SL.library, 
                            propensity_SL.cvControl = list(V = 10L, 
                                                           stratifyCV = T, 
                                                           shuffle = TRUE, 
                                                           validRows = NULL), 
                            include_G_propensity = F,
                            include_W_outcome = T, 
                            TMLE_args_list = TMLE_args_list,
                            SNP_results = 1,
                            parallel = F,
                            ncores = NULL, 
                            progress = F, 
                            verbose = F)[,1, drop = F]
    
    # bias tests 
      # oracle -- gcomp with correct outcome formula 
        # gcomp for 3 levels of X1 
        
        gcomp_oracle_0 = gcomp(outcome_data = data, 
                               family = family,
                               formula = oracle_formula, 
                               pred_ids = which(data$X1 == 0))
        E01  = mean(gcomp_oracle_0$Q1W); E00 = mean(gcomp_oracle_0$Q0W)
        
        gcomp_oracle_1 = gcomp(outcome_data = data, 
                               family = family,
                               formula = oracle_formula, 
                               pred_ids = which(data$X1 == 1))
        E11  = mean(gcomp_oracle_1$Q1W); E10 = mean(gcomp_oracle_1$Q0W)
        
        gcomp_oracle_2 = gcomp(outcome_data = data, 
                               family = family,
                               formula = oracle_formula, 
                               pred_ids = which(data$X1 == 2))
        E21  = mean(gcomp_oracle_2$Q1W); E20 = mean(gcomp_oracle_2$Q0W)
        
        oracle_ATE = c(E01 - E00, E11 - E10, E21 - E20)
        
        oracle_MOR = c( (E01 / (1-E01)) / (E00 / (1 - E00)), 
                        (E11 / (1 - E11)) / (E10 / (1-E10)), 
                        (E21 / (1- E21)) / (E20 / (1 - E20)))
        
     # oracle -- gcomp with correct outcome formula, no other snps included 
        # gcomp for 3 levels of X1 
        
        gcomp_oracle_no_snp_0 = gcomp(outcome_data = data, 
                               family = family,
                               formula = oracle_formula_no_snp, 
                               pred_ids = which(data$X1 == 0))
        E01  = mean(gcomp_oracle_no_snp_0$Q1W); E00 = mean(gcomp_oracle_no_snp_0$Q0W)
        
        gcomp_oracle_no_snp_1 = gcomp(outcome_data = data, 
                               family = family,
                               formula = oracle_formula_no_snp, 
                               pred_ids = which(data$X1 == 1))
        E11  = mean(gcomp_oracle_no_snp_1$Q1W); E10 = mean(gcomp_oracle_no_snp_1$Q0W)
        
        gcomp_oracle_no_snp_2 = gcomp(outcome_data = data, 
                               family = family,
                               formula = oracle_formula_no_snp, 
                               pred_ids = which(data$X1 == 2))
        E21  = mean(gcomp_oracle_no_snp_2$Q1W); E20 = mean(gcomp_oracle_no_snp_2$Q0W)
        
        oracle_no_snp_ATE = c(E01 - E00, E11 - E10, E21 - E20)
        
        oracle_no_snp_MOR = c( (E01 / (1-E01)) / (E00 / (1 - E00)), 
                        (E11 / (1 - E11)) / (E10 / (1-E10)), 
                        (E21 / (1- E21)) / (E20 / (1 - E20)))
        
        # oracle -- gcomp with correct outcome formula at each level of X1 = 0,1,2
        # gcomp for 3 levels of X1 
        
        gcomp_oracle_3lvl_0 = gcomp(outcome_data = data[which(data$X1 == 0), ], 
                                      family = family,
                                      formula = oracle_formula_no_snp)
        E01  = mean(gcomp_oracle_3lvl_0$Q1W); E00 = mean(gcomp_oracle_3lvl_0$Q0W)
        
        gcomp_oracle_3lvl_1 = gcomp(outcome_data = data[which(data$X1 == 1), ], 
                                      family = family,
                                      formula = oracle_formula_no_snp)
        E11  = mean(gcomp_oracle_3lvl_1$Q1W); E10 = mean(gcomp_oracle_3lvl_1$Q0W)
        
        gcomp_oracle_3lvl_2 = gcomp(outcome_data = data[which(data$X1 == 2), ], 
                                      family = family,
                                      formula = oracle_formula_no_snp)
        
        E21  = mean(gcomp_oracle_3lvl_2$Q1W); E20 = mean(gcomp_oracle_3lvl_2$Q0W)
        
        oracle_3lvl_ATE = c(E01 - E00, E11 - E10, E21 - E20)
        
        oracle_3lvl_MOR = c( (E01 / (1-E01)) / (E00 / (1 - E00)), 
                               (E11 / (1 - E11)) / (E10 / (1-E10)), 
                               (E21 / (1- E21)) / (E20 / (1 - E20)))
        
      # gcomp -- compute ACE for SNP X1 = 0, 1, 2 w/o other SNPs, include confounders 
        
        gcomp_0 = gcomp(outcome_data = data[which(data$X1 == 0), ], 
                            family = family,
                            formula = glm_formula_std)
        
        E01  = mean(gcomp_0$Q1W); E00 = mean(gcomp_0$Q0W)
        
        gcomp_1 = gcomp(outcome_data = data[which(data$X1 == 1), ], 
                            family = family,
                            formula = glm_formula_std)
        E11  = mean(gcomp_1$Q1W); E10 = mean(gcomp_1$Q0W)
        
        gcomp_2 = gcomp(outcome_data = data[which(data$X1 == 2), ], 
                            family = family,
                            formula = glm_formula_std)
        E21  = mean(gcomp_2$Q1W); E20 = mean(gcomp_2$Q0W)
        
        gcomp_ATE = c(E01 - E00, E11 - E10, E21 - E20)
        
        gcomp_MOR = c( (E01 / (1-E01)) / (E00 / (1 - E00)), 
                     (E11 / (1 - E11)) / (E10 / (1-E10)), 
                     (E21 / (1- E21)) / (E20 / (1 - E20))
        )  
      # iptw -- compute ACE for SNP X1 = 0, 1, 2 w/o other SNPs, include confounders 
        
        # X1 = 0 
        ipw_weights_0 = ipwpoint(exposure = A, 
                                 family = "binomial", link = "logit", 
                               denominator = ~  age +sex + cohort + ancestry_1 + ancestry_2, 
                               data = data[which(data$X1 == 0), ])$ipw.weights
        
        
        suppressWarnings(iptw_0 <- glm(Y ~ A, data = data[which(data$X1 == 0), ], 
                     family = family, weights = ipw_weights_0))
        
        E01  = (iptw_0$coefficients %>% sum) ; E00 = iptw_0$coefficients[1]
        
        # X1 = 1
        ipw_weights_1 = ipwpoint(exposure = A, 
                                 family = "binomial", link = "logit", 
                                 denominator = ~  age +sex + cohort + ancestry_1 + ancestry_2, 
                                 data = data[which(data$X1 == 1), ])$ipw.weights
        
        
        suppressWarnings(iptw_1 <- glm(Y ~ A, data = data[which(data$X1 == 1), ], 
                                       family = family, weights = ipw_weights_1))
        
        E11  = (iptw_1$coefficients %>% sum) ; E10 = iptw_1$coefficients[1]
        
        # X1 = 2
        ipw_weights_2 = ipwpoint(exposure = A, 
                                 family = "binomial", link = "logit", 
                                 denominator = ~  age +sex + cohort + ancestry_1 + ancestry_2, 
                                 data = data[which(data$X1 == 2), ])$ipw.weights
        
        
        suppressWarnings(iptw_2 <- glm(Y ~ A, data = data[which(data$X1 == 2), ], 
                                       family = family, weights = ipw_weights_2))
        
        E21  = (iptw_2$coefficients %>% sum) ; E20 = iptw_2$coefficients[1]

        
        if(family == "binomial")
        {
          E01 = plogis(E01); E00 = plogis(E00)
          E11 = plogis(E11); E10 = plogis(E10)
          E21 = plogis(E21); E20 = plogis(E20)
        }
        
        iptw_ATE = c(E01 - E00, E11 - E10, E21 - E20) %>% as.numeric
        
        iptw_MOR = c( (E01 / (1-E01)) / (E00 / (1 - E00)), 
                       (E11 / (1 - E11)) / (E10 / (1-E10)), 
                       (E21 / (1- E21)) / (E20 / (1 - E20))
        )  %>% as.numeric

    power = c(tmle_ATE = ifelse(tmle_mod["ATE_EM_pvalue",1] < 0.05, 1, 0) %>% as.numeric,
              tmle_ATE_linear = ifelse(tmle_mod["ATE_EM_lin_pvalue",1] < 0.05, 1, 0)%>% as.numeric,
              tmle_MOR = ifelse(tmle_mod["MOR_EM_pvalue",1] < 0.05, 1, 0)%>% as.numeric,
              tmle_MOR_mult = ifelse(tmle_mod["MOR_EM_mult_pvalue",1] < 0.05, 1, 0)%>% as.numeric,
              oracle_mod = ifelse(summary(oracle_mod)$coefficients["A:X1", 4] < 0.05, 1, 0)%>% as.numeric,
              oracle_no_snp = ifelse(summary(oracle_mod_no_snp)$coefficients["A:X1", 4] < 0.05, 1, 0)%>% as.numeric,
              gxe_1df = glm_1df_power%>% as.numeric,
              gxe_2df = ifelse(lr_test$`Pr(>Chisq)`[2] < 0.05, 1, 0)%>% as.numeric)
    
    ATE_res = rbind(TMLE = c(tmle_mod[1],tmle_mod[2], tmle_mod[3]), 
                    TMLE_linear = c(tmle_mod["ATE_EM_lin_baseline_est",1],
                                    tmle_mod["ATE_EM_lin_baseline_est",1] + tmle_mod["ATE_EM_lin_est",1] , 
                                    tmle_mod["ATE_EM_lin_baseline_est",1] + 2*tmle_mod["ATE_EM_lin_est",1]),
                    oracle_ATE = oracle_ATE,
                    oracle_no_snp_ATE = oracle_no_snp_ATE, 
                    oracle_3lvl_ATE = oracle_3lvl_ATE,
                    gcomp_ATE = gcomp_ATE, 
                    iptw_ATE = iptw_ATE
                    )
    
    MOR_res = rbind(TMLE = c(tmle_mod[7],tmle_mod[8], tmle_mod[9]), 
                    TMLE_mult = c(tmle_mod["MOR_EM_mult_baseline_est",1], 
                                  tmle_mod["MOR_EM_mult_baseline_est",1] * tmle_mod["MOR_EM_mult_est",1], 
                                  tmle_mod["MOR_EM_mult_baseline_est",1] * tmle_mod["MOR_EM_mult_est",1]^2), 
                    oracle_MOR = oracle_MOR,
                    oracle_no_snp_MOR = oracle_no_snp_MOR,
                    oracle_3lvl_MOR = oracle_3lvl_MOR,
                    gcomp_MOR = gcomp_MOR,
                    iptw_MOR = iptw_MOR)
    
    ACE_res = list(signif = power, ATE_res = ATE_res, 
                   MOR_res = MOR_res)
    
    p()
    return(ACE_res)
  }, future.seed = future.seed) %>% t
  
  plan(sequential)
  }
  
  # PROCESS Power 
  
  power = sapply(EM, "[[", 1) %>% t %>% colMeans(na.rm = T) # na.RM is for MOR_mult that sometimes 
  # is NA
  
  # Process bias and MSE estimates 
  
  true_ACEs = A_prop_int_effects$sim_ACEs
  
  if(ACE_type == "ATE")
  {
    true_ACEs[1:3] = ACE
  }else # if controlling MOR
  {
    true_ACEs[4:6] = ACE
  }
  
  return_list = NULL
  
  if(family == "binomial")
  {
    # methods_ATE = EM[[1]]$ATE_res %>% rownames
    # ATE_estimates = data.frame(do.call(rbind, lapply(EM,"[[", 2)) %>% 
    #                              apply(MARGIN = 1, FUN = function(x) return(x - true_ACEs[1:3])) %>% t, 
    #                            method = rep(methods_ATE, times = length(EM)) %>% factor)
    # 
    # colnames(ATE_estimates)[1:3] = paste0("ATE_", 0:2)
    # 
    # ATE_long_estimates = reshape(ATE_estimates, varying = 1:3, direction = "long", sep = "_",
    #                              timevar = "EM_level")
    # ATE_long_estimates$EM_level = ATE_long_estimates$EM_level %>% factor
    
    methods_MOR = EM[[1]]$MOR_res %>% rownames
    MOR_estimates = data.frame(do.call(rbind, lapply(EM,"[[", 3)) %>% 
                                 apply(MARGIN = 1, FUN = function(x) return(x - true_ACEs[4:6])) %>% t, 
                               method = rep(methods_MOR, times = length(EM)) %>% factor)
    
    colnames(MOR_estimates)[1:3] = paste0("MOR_", 0:2)
    
    MOR_long_estimates = reshape(MOR_estimates, varying = 1:3, direction = "long", sep = "_",
                                 timevar = "EM_level")
    MOR_long_estimates$EM_level = MOR_long_estimates$EM_level %>% factor
  
    
    power_plot_data = data.frame(power = power, method = names(power), SNP = "X1")
    
    power_plot = ggplot(data = power_plot_data, mapping = aes(x = SNP, y = power, color = method)) + 
      geom_jitter(size = 2, width = .1) + geom_hline(yintercept = 0.05, color = "red", linetype = "dashed") + 
      geom_hline(yintercept = 0.04, color = "black", linetype = "dashed") + 
      geom_hline(yintercept = 0.06, color = "black", linetype = "dashed")
    
    # ATE_plot = ggplot(data = ATE_long_estimates, mapping = aes(x = EM_level, y = ATE, fill = method)) + 
    #   geom_boxplot() + xlab("Effect Modification Level of SNP X1 = 0,1,2") + 
    #   ylab("Bias of ATE estimate") + geom_hline(yintercept = 0, color = "red", linetype = "dashed")
    
    
    MOR_plot = ggplot(data = MOR_long_estimates, mapping = aes(x = EM_level, y = MOR, fill = method)) + 
      geom_boxplot() + xlab("Effect Modification Level of SNP X1 = 0,1,2") + 
      ylab("Bias of MOR estimate") + geom_hline(yintercept = 0, color = "red", linetype = "dashed")
    
    EM_P = ggarrange(power_plot, MOR_plot, ncol = 2)
    
    # make bias plot 
    bias_mean_data = tapply(MOR_long_estimates$MOR, list(MOR_long_estimates$EM_level, 
                                                    MOR_long_estimates$method), mean)%>% t %>% data.frame 
    
    bias_se_data = tapply(MOR_long_estimates$MOR, list(MOR_long_estimates$EM_level, 
                                                       MOR_long_estimates$method), sd)%>% t %>% data.frame 
    # bias_data$method = rownames(bias_data) %>% factor
    # 
    # colnames(bias_data) = c("SNP_0", "SNP_1", "SNP_2", "method")
    # 
    # bias_data = bias_data %>% 
    #   pivot_longer(cols = starts_with("SNP_"), names_to = "SNP_lvl", values_to = "Bias") %>% data.frame
    # 
    # bias_data$SNP_lvl = substr(bias_data$SNP_lvl, 5, 5) %>% as.numeric
    # 
    # bias_plot = ggplot(data = bias_data, mapping = aes(x = SNP_lvl, y = Bias, color = method)) + 
    #   geom_line() + geom_point() + xlab("SNP Level") + scale_x_continuous(breaks = c(0,1,2)) + 
    #   ylab("| Bias |") + scale_color_brewer(palette = "Set1")
    
    return_list = list(power = power, EM_P = EM_P, 
                       true_ACEs = true_ACEs, 
                       A_effects = A_prop_int_effects$A_effects, 
                       bias_mean_data = bias_mean_data, 
                       bias_se_data = bias_se_data)
  }else
  {
    methods_ATE = EM[[1]]$ATE_res %>% rownames
    ATE_estimates = data.frame(do.call(rbind, lapply(EM,"[[", 2)) %>% 
                                 apply(MARGIN = 1, FUN = function(x) return(x - true_ACEs[1:3])) %>% t, 
                               method = rep(methods_ATE, times = length(EM)) %>% factor)
    
    colnames(ATE_estimates)[1:3] = paste0("ATE_", 0:2)
    
    ATE_long_estimates = reshape(ATE_estimates, varying = 1:3, direction = "long", sep = "_",
                                 timevar = "EM_level")
    ATE_long_estimates$EM_level = ATE_long_estimates$EM_level %>% factor
    
    
    power = power[-which(names(power) %in% c("tmle_MOR", "tmle_MOR_mult"))]
    power_plot_data = data.frame(power = power, method = names(power), SNP = "X1")
    
    power_plot = ggplot(data = power_plot_data, mapping = aes(x = SNP, y = power, color = method)) + 
      geom_jitter(size = 2, width = .1) + geom_hline(yintercept = 0.05, color = "red", linetype = "dashed") + 
      geom_hline(yintercept = 0.04, color = "black", linetype = "dashed") + 
      geom_hline(yintercept = 0.06, color = "black", linetype = "dashed")
    
    ATE_plot = ggplot(data = ATE_long_estimates, mapping = aes(x = EM_level, y = ATE, fill = method)) + 
      geom_boxplot() + xlab("Effect Modification Level of SNP X1 = 0,1,2") + 
      ylab("Bias of ATE estimate") + geom_hline(yintercept = 0, color = "red", linetype = "dashed")
    
    bias_mean_data = tapply(ATE_long_estimates$ATE, list(ATE_long_estimates$EM_level, 
                                                         ATE_long_estimates$method), mean)%>% t %>% data.frame 
    
    bias_se_data = tapply(ATE_long_estimates$ATE, list(ATE_long_estimates$EM_level, 
                                                       ATE_long_estimates$method), sd)%>% t %>% data.frame 
    
    
    EM_P = ggarrange(power_plot, ATE_plot, ncol = 2)
    
    return_list = list(power = power, EM_P = EM_P, 
                       true_ACEs = true_ACEs, 
                       A_effects = A_prop_int_effects$A_effects, 
                       bias_mean_data = bias_mean_data, 
                       bias_se_data = bias_se_data)
  }
  
  
  return(return_list)
}

output_results = function(sim)
{
  cat("Gcomp results for ATE and MOR -- USING LINEAR ADJUSTMENTS FOR CONFOUNDING\n")
  print(sim$true_ACEs)
  cat("\nA effects: ", paste(sim$A_effects %>% round(4), collapse = ", "))
  cat("\n\n")
  print(sim$power)
  print(sim$EM_P)
}


# function for comparing power results
process_sim_results = function(sim_results)
{
  
  rownames = paste0("EM_lvl_", 1:length(sim_results))
  ACE_effects = sapply(sim_results, "[[", 3) %>% t
  A_effects = sapply(sim_results, "[[", 4) %>% t
  
  rownames(ACE_effects) = rownames
  rownames(A_effects) = rownames
  colnames(A_effects) = c("X1_0","X1_1", "X1_2")
  
  power = sapply(sim_results, "[[", 1) %>% t 
  rownames(power) = rownames
  
  plot_data = data.frame(power = as.vector(power), 
                         method = rep(colnames(power), each = nrow(power)), 
                         ACE_lvl = rep(1:nrow(power), times = ncol(power)))
  
  plot = ggplot(data = plot_data, mapping = aes(x = ACE_lvl, y = power, color = method)) + 
    geom_point() + geom_line() + 
    geom_hline(yintercept = 0.05, color = "red", linetype = "dashed") + 
    geom_hline(yintercept = 0.04, color = "black", linetype = "dashed") + 
    geom_hline(yintercept = 0.06, color = "black", linetype = "dashed") + 
    scale_color_brewer(palette = "Set1")
  
  return(list(ACE_effects = ACE_effects, 
              A_effects = A_effects, 
              plot = plot,
              power_data = power))
}

process_power_plots = function(sim)
{
  power_data = sapply(sim, "[[", 1) %>% t
  
  power_data = cbind(sim_num = 1:nrow(power_data), power_data)
  
  power_data = power_data %>% data.frame %>% pivot_longer(cols =2:ncol(power_data), 
                                                          names_to = "method", 
                                                          values_to = "power") %>% data.frame
  
  power_plot = ggplot(data = power_data, mapping = aes(x = sim_num, y = power, color = method)) + 
    geom_point() + geom_line() + xlab("ACE sim level") + ylab("Power") +
    scale_x_continuous(breaks = unique(power_data$sim_num)) + 
    scale_color_brewer(palette = "Set1") + geom_hline(yintercept = 0.05, color = "red", 
                                                      linetype = "dashed")
  
  return(power_plot)
  
}



