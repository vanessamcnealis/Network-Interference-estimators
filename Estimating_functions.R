
# Functions that collect the neighbourhood covariate
# and treatment vectors for a given node --------------------------------------

x.alters <- function(i, design.matrix){
  x <- A[i,]
  X.N <- design.matrix[x==1,]
  return(rbind(design.matrix[i,],X.N))
}

z.alters <- function(i, Z){
  x <- A[i,]
  Z.N <- Z[x==1]
  return(c(Z[i],Z.N))
}

avg.po <- function(pos, geex.object){
  avg <- geex.object@estimates[pos]
  vcov <- geex.object@vcov[pos, pos]
  return(c(avg, vcov))
}

contrast <- function(pos1, pos2, geex.object){
  avg1 <- geex.object@estimates[pos1]
  avg2 <- geex.object@estimates[pos2]
  diff <- avg1 - avg2
  
  vcov <- geex.object@vcov[c(pos1,pos2), c(pos1, pos2)]
  var.pos1 <- vcov[1,1]
  var.pos2 <- vcov[2,2]
  cov.1.2 <- vcov[1,2]
  
  var.contrast <- var.pos1 + var.pos2 - 2*cov.1.2
  
  return(c(diff, var.contrast))
}

# Conditional joint propensity score --------------------------------------

g <- function(b, i, data, exposure.model, Z){
  
  gamma.hat <- as.numeric(fixef(exposure.model))
  
  var.b <- as.numeric(VarCorr(exposure.model))
  
  design.matrix <- model.matrix(exposure.model)
  
  # i-th row of the adjacency matrix
  x <- A[i,]
  X.N <- rbind(design.matrix[i,], design.matrix[x==1,])
  
  # Predict neighbourhood level propensity scores
  p.N <- inv.logit(as.numeric(X.N %*% as.matrix(gamma.hat)) + rep(b, nrow(X.N)))
  
  # Neighborhood treatments
  Z.N <- c(Z[i], Z[x==1])
  
  return(prod(p.N^Z.N*(1-p.N)^(1-Z.N)))
}


# Function which computes the marginal joint propensity score -------------

int.GH <- function(i, exposure.model, Z){
  
  var.b <- as.numeric(VarCorr(exposure.model))
  
  w <- gauss.quad(10, kind="hermite")
  
  gx <- sapply(w$nodes*sqrt(2*var.b), FUN = g, i = i, 
               exposure.model = exposure.model,
               data=data, Z=Z)
  
  return(sum(w$weights*gx/sqrt(pi)))
}

# Initial values for IPW estimator of the avg potential outcome under treatment
# coverage alpha ----------------------------------------

initial_values_ipw <- function(data, propensity_formula, alpha){
  
  mod.ps <- glmer(propensity_formula,
                  data = data,
                  family = binomial(link = 'logit'))
  
  gamma.hat <- as.numeric(fixef(mod.ps))
  var.b <- as.numeric(VarCorr(mod.ps))
  design.matrix <- model.matrix(mod.ps)
  
  Z <- data$Z
  
  f <- sapply(data$id, FUN = int.GH, exposure.model = mod.ps, Z=Z)
  
  temp <- data.frame(component.id = as.numeric(data$component.id), Y=as.numeric(data$Y), 
                     Z = as.numeric(data$Z),
                     degree = as.numeric(data$degree), k_treated=as.numeric(data$k_treated))
  
  num1 <- (alpha[1]^temp$k_treated)*((1-alpha[1])^(temp$degree-temp$k_treated))
  temp$ipw0.alpha1 <- (temp$Y*as.numeric(temp$Z==0)*num1)/f
  temp$ipw1.alpha1 <- (temp$Y*as.numeric(temp$Z==1)*num1)/f
  temp$ipw.alpha1 <- (1-alpha[1])*temp$ipw0.alpha1 + alpha[1]*temp$ipw1.alpha1
  num2 <- (alpha[2]^temp$k_treated)*((1-alpha[2])^(temp$degree-temp$k_treated))
  temp$ipw0.alpha2 <- (temp$Y*as.numeric(temp$Z==0)*num2)/f
  temp$ipw1.alpha2 <- (temp$Y*as.numeric(temp$Z==1)*num2)/f
  temp$ipw.alpha2 <- (1-alpha[2])*temp$ipw0.alpha2 + alpha[2]*temp$ipw1.alpha2
  num3 <- (alpha[3]^temp$k_treated)*((1-alpha[3])^(temp$degree-temp$k_treated))
  temp$ipw0.alpha3 <- (temp$Y*as.numeric(temp$Z==0)*num3)/f
  temp$ipw1.alpha3 <- (temp$Y*as.numeric(temp$Z==1)*num3)/f
  temp$ipw.alpha3 <- (1-alpha[3])*temp$ipw0.alpha3 + alpha[3]*temp$ipw1.alpha3
  
  Var <- c("ipw0.alpha1", "ipw1.alpha1", "ipw.alpha1", 
           "ipw0.alpha2", "ipw1.alpha2", "ipw.alpha2",
           "ipw0.alpha3", "ipw1.alpha3", "ipw.alpha3")
  
  component_summaries_ipw <- temp %>% group_by(component.id) %>%
    summarise_at(.vars = vars(all_of(Var)),
                 .funs = list(Component = ~mean(.) ) )
  
  ipw0.alpha1 <- mean(component_summaries_ipw$ipw0.alpha1_Component)
  ipw1.alpha1 <- mean(component_summaries_ipw$ipw1.alpha1_Component)
  ipw.alpha1 <- mean(component_summaries_ipw$ipw.alpha1_Component)
  ipw0.alpha2 <- mean(component_summaries_ipw$ipw0.alpha2_Component)
  ipw1.alpha2 <- mean(component_summaries_ipw$ipw1.alpha2_Component)
  ipw.alpha2 <- mean(component_summaries_ipw$ipw.alpha2_Component)
  ipw0.alpha3 <- mean(component_summaries_ipw$ipw0.alpha3_Component)
  ipw1.alpha3 <- mean(component_summaries_ipw$ipw1.alpha3_Component)
  ipw.alpha3 <- mean(component_summaries_ipw$ipw.alpha3_Component)
  
  return(c(gamma.hat, var.b, ipw0.alpha1, ipw1.alpha1, ipw.alpha1,
           ipw0.alpha2, ipw1.alpha2, ipw.alpha2,
           ipw0.alpha3, ipw1.alpha3, ipw.alpha3))
  
}


# Estimating functions for the IPW estimator ------------------------------

ipw_estfun <- function(data, models, alpha, N, m){ 
  id <- data$id
  Z <- data$Z
  Y <- data$Y
  X.N <- data$X.N
  Z.N <- data$Z.N
  k_treated <- data$k_treated
  degree <- data$degree
  Xe <- grab_design_matrix(data = data,
                           rhs_formula = grab_fixed_formula(models$e))
  gamma_pos <- 1:(ncol(Xe))
  varb_pos <- length(gamma_pos)+1
  # Grab estimating functions
  e_scores <- grab_psiFUN(models$e, data)
  
  num1 <- (alpha[1]^k_treated)*((1-alpha[1])^(degree-k_treated))
  num2 <- (alpha[2]^k_treated)*((1-alpha[2])^(degree-k_treated))
  num3 <- (alpha[3]^k_treated)*((1-alpha[3])^(degree-k_treated))
  
  g <- function(X.N, Z.N, b, coef){
    
    # Predict neighborhood level propensity scores
    p.n <- inv.logit(X.N %*% coef + b)
    
    return(prod(p.n^Z.N*(1-p.n)^(1-Z.N)))
  }
  
  int.GH <- function(X.N, Z.N, coef, var.b){
    w <- gauss.quad(10, kind="hermite")
    
    gx <- sapply(w$nodes*sqrt(2*var.b), FUN = g, X.N = X.N,
                 Z.N = Z.N, coef=coef)
    return(sum(w$weights*gx/sqrt(pi)))
  }
  
  
  function(theta) {
    coef <- theta[gamma_pos]
    var.b <- (theta[varb_pos])^2
    e <- map2_dbl(X.N, Z.N, int.GH, coef = coef, var.b = var.b)
    ipw0.alpha1 <- mean((Y*as.numeric(Z==0)*num1)/e)
    ipw1.alpha1 <- mean((Y*as.numeric(Z==1)*num1)/e)
    ipw.alpha1 <- alpha[1]*ipw1.alpha1 + (1-alpha[1])*ipw0.alpha1
    ipw0.alpha2 <- mean((Y*as.numeric(Z==0)*num2)/e)
    ipw1.alpha2 <- mean((Y*as.numeric(Z==1)*num2)/e)
    ipw.alpha2 <- alpha[2]*ipw1.alpha2 + (1-alpha[2])*ipw0.alpha2
    ipw0.alpha3 <- mean((Y*as.numeric(Z==0)*num3)/e)
    ipw1.alpha3 <- mean((Y*as.numeric(Z==1)*num3)/e)
    ipw.alpha3 <- alpha[3]*ipw1.alpha3 + (1-alpha[3])*ipw0.alpha3
    
    return(c(e_scores(theta[c(gamma_pos, varb_pos)]),
             ipw0.alpha1 - theta[length(theta)-8],
             ipw1.alpha1 - theta[length(theta)-7],
             ipw.alpha1 - theta[length(theta)-6],
             ipw0.alpha2 - theta[length(theta)-5],
             ipw1.alpha2 - theta[length(theta)-4],
             ipw.alpha2 - theta[length(theta)-3],
             ipw0.alpha3 - theta[length(theta)-2],
             ipw1.alpha3 - theta[length(theta)-1],
             ipw.alpha3 - theta[length(theta)]))
  }
}

estimate_ipw <- function(data, propensity_formula, init, alpha, m){
  N <- nrow(data)
  e_model <- glmer(propensity_formula,
                   data = data,
                   family = binomial(link = 'logit'))
  models <- list(e = e_model)
  m_estimate(estFUN = ipw_estfun, data = data, units="component.id",
             root_control = setup_root_control(start = init),
             outer_args = list(models = models,
                               alpha = alpha, N=N, m=m))
}

ipw_estimation <- function(data, propensity_formula, alpha){
  
  mod.ps <- glmer(propensity_formula,
                  data = data,
                  family = binomial(link = 'logit'))
  
  # Collect the neighbourhood covariate and treatment vectors for each node
  design.matrix <- model.matrix(mod.ps)
  
  #create a copy of the data set to alter it
  
  Z <- data$Z
  
  data$X.N  <- sapply(data$id, FUN = x.alters, 
                     design.matrix=design.matrix)
  data$Z.N  <- sapply(data$id, FUN = z.alters, Z=Z)
  
  # Compute initial values for the M-estimation routine (in this case,
  # they correspond to the IPW point estimates of the average potential
  # outcomes)
  
  # If q is the number of columns of the design matrix, then the first
  # q elements of init_ipw correspond to the MLE of gamma, the parameter
  # vector of the logistic regression model, the q+1 element to the random
  # intercept variance, and the last 6 elements correspond to the ipw
  # estimates of the average potential outcomes under treatment and no
  # treatment for alpha1, alpha2 and alpha3
  
  init_ipw <- initial_values_ipw(data=data, propensity_formula, alpha)
  
  # Call the M-estimation routine with the above initial values
  # This might take some time
  ipw.estimates <- estimate_ipw(data = data, 
                                propensity_formula = propensity_formula, 
                                alpha = alpha,
                                init = init_ipw,
                                m=m)
  
  # Retrieve the quantities of interest and calculate causal contrasts
  
  nparm <- ncol(design.matrix)+1
  pos.0.alpha1 <- nparm+1
  pos.1.alpha1 <- nparm+2
  pos.alpha1 <- nparm+3
  pos.0.alpha2 <- nparm+4
  pos.1.alpha2 <- nparm+5
  pos.alpha2 <- nparm+6
  pos.0.alpha3 <- nparm+7
  pos.1.alpha3 <- nparm+8
  pos.alpha3 <- nparm + 9
  
  # Table with estimates of average potential outcomes
  
  temp <- avg.po(pos = pos.0.alpha1, ipw.estimates)
  row1 <- data.frame(estimand = "Average potential outcome z=0", alpha = alpha[1],
                     estimate = temp[1], variance = temp[2], 
                     CI_lb = temp[1] - qnorm(0.975)*sqrt(temp[2]),
                     CI_ub = temp[1] + qnorm(0.975)*sqrt(temp[2]))
  temp <- avg.po(pos = pos.1.alpha1, ipw.estimates)
  row2 <- data.frame(estimand = "Average potential outcome z=1", alpha = alpha[1],
                     estimate = temp[1], variance = temp[2], 
                     CI_lb = temp[1] - qnorm(0.975)*sqrt(temp[2]),
                     CI_ub = temp[1] + qnorm(0.975)*sqrt(temp[2]))
  temp <- avg.po(pos = pos.alpha1, ipw.estimates)
  row3 <- data.frame(estimand = "Marginal average potential outcome", alpha = alpha[1],
                     estimate = temp[1], variance = temp[2], 
                     CI_lb = temp[1] - qnorm(0.975)*sqrt(temp[2]),
                     CI_ub = temp[1] + qnorm(0.975)*sqrt(temp[2]))
  
  temp <- avg.po(pos = pos.0.alpha2, ipw.estimates)
  row4 <- data.frame(estimand = "Average potential outcome z=0", alpha = alpha[2],
                     estimate = temp[1], variance = temp[2], 
                     CI_lb = temp[1] - qnorm(0.975)*sqrt(temp[2]),
                     CI_ub = temp[1] + qnorm(0.975)*sqrt(temp[2]))
  temp <- avg.po(pos = pos.1.alpha2, ipw.estimates)
  row5<- data.frame(estimand = "Average potential outcome z=1", alpha = alpha[2],
                     estimate = temp[1], variance = temp[2], 
                     CI_lb = temp[1] - qnorm(0.975)*sqrt(temp[2]),
                     CI_ub = temp[1] + qnorm(0.975)*sqrt(temp[2]))
  temp <- avg.po(pos = pos.alpha2, ipw.estimates)
  row6 <- data.frame(estimand = "Marginal average potential outcome", alpha = alpha[2],
                     estimate = temp[1], variance = temp[2], 
                     CI_lb = temp[1] - qnorm(0.975)*sqrt(temp[2]),
                     CI_ub = temp[1] + qnorm(0.975)*sqrt(temp[2]))
  
  temp <- avg.po(pos = pos.0.alpha3, ipw.estimates)
  row7 <- data.frame(estimand = "Average potential outcome z=0", alpha = alpha[3],
                     estimate = temp[1], variance = temp[2], 
                     CI_lb = temp[1] - qnorm(0.975)*sqrt(temp[2]),
                     CI_ub = temp[1] + qnorm(0.975)*sqrt(temp[2]))
  temp <- avg.po(pos = pos.1.alpha3, ipw.estimates)
  row8 <- data.frame(estimand = "Average potential outcome z=1", alpha = alpha[3],
                     estimate = temp[1], variance = temp[2], 
                     CI_lb = temp[1] - qnorm(0.975)*sqrt(temp[2]),
                     CI_ub = temp[1] + qnorm(0.975)*sqrt(temp[2]))
  temp <- avg.po(pos = pos.alpha3, ipw.estimates)
  row9 <- data.frame(estimand = "Marginal average potential outcome", alpha = alpha[3],
                     estimate = temp[1], variance = temp[2], 
                     CI_lb = temp[1] - qnorm(0.975)*sqrt(temp[2]),
                     CI_ub = temp[1] + qnorm(0.975)*sqrt(temp[2]))
  
  table1 <- rbind(row1, row2, row3, row4, row5, row6, row7, row8, row9)
  
  # Table with estimates of relevant contrasts
  # Direct effect
  
  temp <- contrast(pos.1.alpha1, pos.0.alpha1, ipw.estimates)
  row1 <- data.frame(estimand = "Direct effect", alpha0 = alpha[1], alpha1 = alpha[1],
                     estimate = temp[1], variance = temp[2], 
                     CI_lb = temp[1] - qnorm(0.975)*sqrt(temp[2]),
                     CI_ub = temp[1] + qnorm(0.975)*sqrt(temp[2]))
  temp <- contrast(pos.1.alpha2, pos.0.alpha2, ipw.estimates)
  row2 <- data.frame(estimand = "Direct effect", alpha0 = alpha[2], alpha1 = alpha[2],
                     estimate = temp[1], variance = temp[2], 
                     CI_lb = temp[1] - qnorm(0.975)*sqrt(temp[2]),
                     CI_ub = temp[1] + qnorm(0.975)*sqrt(temp[2]))
  temp <- contrast(pos.1.alpha3, pos.0.alpha3, ipw.estimates)
  row3 <- data.frame(estimand = "Direct effect", alpha0 = alpha[3], alpha1 = alpha[3],
                     estimate = temp[1], variance = temp[2], 
                     CI_lb = temp[1] - qnorm(0.975)*sqrt(temp[2]),
                     CI_ub = temp[1] + qnorm(0.975)*sqrt(temp[2]))
  
  temp <- contrast(pos.0.alpha2, pos.0.alpha1, ipw.estimates)
  row4 <- data.frame(estimand = "Indirect effect", alpha0 = alpha[1], alpha1 = alpha[2],
                     estimate = temp[1], variance = temp[2], 
                     CI_lb = temp[1] - qnorm(0.975)*sqrt(temp[2]),
                     CI_ub = temp[1] + qnorm(0.975)*sqrt(temp[2]))
  temp <- contrast(pos.0.alpha3, pos.0.alpha1, ipw.estimates)
  row5 <- data.frame(estimand = "Indirect effect", alpha0 = alpha[1], alpha1 = alpha[3],
                     estimate = temp[1], variance = temp[2], 
                     CI_lb = temp[1] - qnorm(0.975)*sqrt(temp[2]),
                     CI_ub = temp[1] + qnorm(0.975)*sqrt(temp[2]))
  temp <- contrast(pos.0.alpha3, pos.0.alpha2, ipw.estimates)
  row6 <- data.frame(estimand = "Indirect effect", alpha0 = alpha[2], alpha1 = alpha[3],
                     estimate = temp[1], variance = temp[2], 
                     CI_lb = temp[1] - qnorm(0.975)*sqrt(temp[2]),
                     CI_ub = temp[1] + qnorm(0.975)*sqrt(temp[2]))
  
  temp <- contrast(pos.1.alpha2, pos.0.alpha1, ipw.estimates)
  row7 <- data.frame(estimand = "Total effect", alpha0 = alpha[1], alpha1 = alpha[2],
                     estimate = temp[1], variance = temp[2], 
                     CI_lb = temp[1] - qnorm(0.975)*sqrt(temp[2]),
                     CI_ub = temp[1] + qnorm(0.975)*sqrt(temp[2]))
  temp <- contrast(pos.1.alpha3, pos.0.alpha1, ipw.estimates)
  row8 <- data.frame(estimand = "Total effect", alpha0 = alpha[1], alpha1 = alpha[3],
                     estimate = temp[1], variance = temp[2], 
                     CI_lb = temp[1] - qnorm(0.975)*sqrt(temp[2]),
                     CI_ub = temp[1] + qnorm(0.975)*sqrt(temp[2]))
  temp <- contrast(pos.1.alpha3, pos.0.alpha2, ipw.estimates)
  row9 <- data.frame(estimand = "Total effect", alpha0 = alpha[2], alpha1 = alpha[3],
                     estimate = temp[1], variance = temp[2], 
                     CI_lb = temp[1] - qnorm(0.975)*sqrt(temp[2]),
                     CI_ub = temp[1] + qnorm(0.975)*sqrt(temp[2]))
  
  temp <- contrast(pos.alpha2, pos.alpha1, ipw.estimates)
  row10 <- data.frame(estimand = "Overall effect", alpha0 = alpha[1], alpha1 = alpha[2],
                     estimate = temp[1], variance = temp[2], 
                     CI_lb = temp[1] - qnorm(0.975)*sqrt(temp[2]),
                     CI_ub = temp[1] + qnorm(0.975)*sqrt(temp[2]))
  temp <- contrast(pos.alpha3, pos.alpha1, ipw.estimates)
  row11 <- data.frame(estimand = "Overall effect", alpha0 = alpha[1], alpha1 = alpha[3],
                     estimate = temp[1], variance = temp[2], 
                     CI_lb = temp[1] - qnorm(0.975)*sqrt(temp[2]),
                     CI_ub = temp[1] + qnorm(0.975)*sqrt(temp[2]))
  temp <- contrast(pos.alpha3, pos.alpha2, ipw.estimates)
  row12 <- data.frame(estimand = "Overall effect", alpha0 = alpha[2], alpha1 = alpha[3],
                     estimate = temp[1], variance = temp[2], 
                     CI_lb = temp[1] - qnorm(0.975)*sqrt(temp[2]),
                     CI_ub = temp[1] + qnorm(0.975)*sqrt(temp[2]))
  
  table2 <- rbind(row1, row2, row3, row4, row5, row6, row7, row8, row9,
                  row10, row11, row12)
  
  return(list(table1, table2))
  
}

# Estimating functions for the REG estimator ------------------------------

reg_estfun <- function(data, models, alpha, counterfactual_prop_treated, m, N){ 
  
  Xm <- grab_design_matrix(data = data,
                           rhs_formula = grab_fixed_formula(models$m))
  beta_pos <- 1:(ncol(Xm))
  
  Xm_df <- data.frame(Xm)
  Xm_df$id <- data$id
  Xm_df <- left_join(Xm_df, counterfactual_prop_treated, by = "id")
  
  Xm1 <- Xm_df %>% mutate(prop_treated = counterfactual_prop_treated,
                          Z=1) %>%
    mutate(Z.prop_treated = Z*prop_treated) %>% 
    dplyr::select(-c(id, counterfactual_prop_treated, k_treated, degree)) %>% 
    as.matrix()
  colnames(Xm1) <- colnames(Xm)
  
  Xm0 <- Xm_df %>% mutate(prop_treated = counterfactual_prop_treated,
                          Z=0) %>%
    mutate(Z.prop_treated = Z*prop_treated) %>% 
    dplyr::select(-c(id, counterfactual_prop_treated, k_treated, degree)) %>% 
    as.matrix()
  colnames(Xm0) <- colnames(Xm)
  
  # Grab estimating functions
  m_scores <- grab_psiFUN(models$m, data)
  
  Xm_df <- Xm_df %>% mutate(weights1 = dbinom(x = k_treated, size = degree, 
                                              prob = alpha[1]),
                            weights2 = dbinom(x = k_treated, size = degree, 
                                              prob = alpha[2]),
                            weights3 = dbinom(x = k_treated, size = degree, 
                                              prob = alpha[3]))
  function(theta) {
    beta <- theta[beta_pos]
    pred.0 <- Xm0 %*% beta
    pred.1 <- Xm1 %*% beta
    Xm_df <- Xm_df %>% add_column(pred.0, pred.1)
    temp <- Xm_df %>% group_by(id) %>% summarise(reg.0.alpha1 = sum(weights1*pred.0),
                                                 reg.1.alpha1 = sum(weights1*pred.1),
                                                 reg.alpha1 = alpha[1]*reg.1.alpha1 + (1-alpha[1])*reg.0.alpha1,
                                                 reg.0.alpha2 = sum(weights2*pred.0),
                                                 reg.1.alpha2 = sum(weights2*pred.1),
                                                 reg.alpha2 = alpha[2]*reg.1.alpha2 + (1-alpha[2])*reg.0.alpha2,
                                                 reg.0.alpha3 = sum(weights3*pred.0),
                                                 reg.1.alpha3 = sum(weights3*pred.1),
                                                 reg.alpha3 = alpha[3]*reg.1.alpha3 + (1-alpha[3])*reg.0.alpha3) 
    reg0.alpha1 <- mean(temp$reg.0.alpha1)
    reg1.alpha1 <- mean(temp$reg.1.alpha1)
    reg.alpha1 <- mean(temp$reg.alpha1)
    reg0.alpha2 <- mean(temp$reg.0.alpha2)
    reg1.alpha2 <- mean(temp$reg.1.alpha2)
    reg.alpha2 <- mean(temp$reg.alpha2)
    reg0.alpha3 <- mean(temp$reg.0.alpha3)
    reg1.alpha3 <- mean(temp$reg.1.alpha3)
    reg.alpha3 <- mean(temp$reg.alpha3)
    
    return(c(m_scores(theta[beta_pos]),
             reg0.alpha1 - theta[length(theta)-8],
             reg1.alpha1 - theta[length(theta)-7],
             reg.alpha1 - theta[length(theta)-6],
             reg0.alpha2 - theta[length(theta)-5],
             reg1.alpha2 - theta[length(theta)-4],
             reg.alpha2 - theta[length(theta)-3],
             reg0.alpha3 - theta[length(theta)-2],
             reg1.alpha3 - theta[length(theta)-1],
             reg.alpha3 - theta[length(theta)]))
  }
}


estimate_reg <- function(data, outcome_formula, init, alpha, 
                         counterfactual_prop_treated, m){
  m_model <- glm(outcome_formula,
                 data = data, family = "gaussian")
  models <- list(m = m_model)
  N <- nrow(data)
  m_estimate(estFUN = reg_estfun, data = data, units="component.id",
             root_control = setup_root_control(start = init),
             outer_args = list(models = models,
                               alpha = alpha,
                               counterfactual_prop_treated = counterfactual_prop_treated,
                               m=m,
                               N=N))
}

# Initial values for REG estimator of the avg potential outcome under treatment
# coverage alpha ----------------------------------------

initial_values_reg <- function(data, outcome_formula, counterfactual, alpha){
  
  mod <- lm(outcome_formula,
              data = data)
  beta.hat <- as.numeric(coef(mod))
  
  Xm <- model.matrix(object = outcome_formula, data = data)
  
  Xm_df <- data.frame(Xm)
  Xm_df$id <- data$id
  Xm_df$component.id <-  data$component.id
  Xm_df <- left_join(Xm_df, counterfactual, by = "id")
  
  Xm1 <- Xm_df %>% mutate(prop_treated = counterfactual_prop_treated,
                          Z=1) %>%
    mutate(Z.prop_treated = Z*prop_treated) %>% 
    dplyr::select(-c(id, component.id, counterfactual_prop_treated, k_treated, degree)) %>% 
    as.matrix()
  colnames(Xm1) <- colnames(Xm)
  
  Xm0 <- Xm_df %>% mutate(prop_treated = counterfactual_prop_treated,
                          Z=0) %>%
    mutate(Z.prop_treated = Z*prop_treated) %>% 
    dplyr::select(-c(id, component.id, counterfactual_prop_treated, k_treated, degree)) %>% 
    as.matrix()
  colnames(Xm0) <- colnames(Xm)
  
  Xm_df <- Xm_df %>% mutate(weights1 = dbinom(x = k_treated, size = degree, 
                                             prob = alpha[1]),
                            weights2 = dbinom(x = k_treated, size = degree, 
                                              prob = alpha[2]),
                            weights3 = dbinom(x = k_treated, size = degree, 
                                              prob = alpha[3]))
  
  pred.0 <- Xm0 %*% beta.hat
  pred.1 <- Xm1 %*% beta.hat
  Xm_df <- Xm_df %>% add_column(pred.0, pred.1)
  temp <- Xm_df %>% group_by(id) %>% dplyr::summarise(reg.0.alpha1 = sum(weights1*pred.0),
                                                      reg.1.alpha1 = sum(weights1*pred.1),
                                                      reg.alpha1 = alpha[1]*reg.1.alpha1 + (1-alpha[1])*reg.0.alpha1,
                                                      reg.0.alpha2 = sum(weights2*pred.0),
                                                      reg.1.alpha2 = sum(weights2*pred.1),
                                                      reg.alpha2 = alpha[2]*reg.1.alpha2 + (1-alpha[2])*reg.0.alpha2,
                                                      reg.0.alpha3 = sum(weights3*pred.0),
                                                      reg.1.alpha3 = sum(weights3*pred.1),
                                                      reg.alpha3 = alpha[3]*reg.1.alpha3 + (1-alpha[3])*reg.0.alpha3) 
  
  temp <- temp %>% left_join(., data %>% dplyr::select(id, component.id), by = "id")
  
  Var <- c("reg.0.alpha1", "reg.1.alpha1", "reg.alpha1", 
           "reg.0.alpha2", "reg.1.alpha2", "reg.alpha2",
           "reg.0.alpha3", "reg.1.alpha3", "reg.alpha3")
  
  component_summaries_reg <- temp %>% group_by(component.id) %>%
    summarise_at(.vars = vars(all_of(Var)),
                 .funs = list(Component = ~mean(.) ) )
  
  reg0.alpha1 <- mean(component_summaries_reg$reg.0.alpha1_Component)
  reg1.alpha1 <- mean(component_summaries_reg$reg.1.alpha1_Component)
  reg.alpha1 <- mean(component_summaries_reg$reg.alpha1_Component)
  reg0.alpha2 <- mean(component_summaries_reg$reg.0.alpha2_Component)
  reg1.alpha2 <- mean(component_summaries_reg$reg.1.alpha2_Component)
  reg.alpha2 <- mean(component_summaries_reg$reg.alpha2_Component)
  reg0.alpha3 <- mean(component_summaries_reg$reg.0.alpha3_Component)
  reg1.alpha3 <- mean(component_summaries_reg$reg.1.alpha3_Component)
  reg.alpha3 <- mean(component_summaries_reg$reg.alpha3_Component)
  
  return(c(beta.hat, reg0.alpha1, reg1.alpha1, reg.alpha1,
           reg0.alpha2, reg1.alpha2, reg.alpha2,
           reg0.alpha3, reg1.alpha3, reg.alpha3))
  
}

reg_estimation <- function(data, outcome_formula, alpha){
  
  Xm <- model.matrix(object = outcome_formula, data = data)
  
  counterfactual <- data %>% uncount(degree+1) %>% 
    dplyr::select(-k_treated) %>%
    mutate(prop_treated = 0) %>% group_by(id) %>% 
    mutate_at(vars(prop_treated), ~. + 0:(first(degree)))  %>%
    rename(k_treated = prop_treated) %>% 
    mutate(prop_treated = ifelse(degree==0, 0, k_treated/degree)) %>%
    dplyr::select(id, prop_treated, k_treated, degree) %>% 
    rename(counterfactual_prop_treated = prop_treated)

  
  # Compute initial values for the M-estimation routine (in this case,
  # they correspond to the REG point estimates of the average potential
  # outcomes)
  
  # If p is the number of columns of the design matrix, then the first
  # p elements of init_reg correspond to the MLE of beta, the parameter
  # vector of the outcome regression model, and the last 9 elements correspond 
  # to the reg estimates of the (marginal) average potential outcomes 
  # under treatment and no treatment for alpha1, alpha2 and alpha3
  
  init_reg <- initial_values_reg(data=data, outcome_formula = outcome_formula,
                                 counterfactual = counterfactual, alpha=alpha)
  
  # Call the M-estimation routine with the above initial values
  # This might take some time
  reg.estimates <- estimate_reg(data = data, 
                                outcome_formula = outcome_formula, 
                                alpha = alpha,
                                init = init_reg,
                                counterfactual_prop_treated=counterfactual,
                                m=m)
  
  # Retrieve the quantities of interest and calculate causal contrasts
  
  nparm <- ncol(Xm)
  pos.0.alpha1 <- nparm+1
  pos.1.alpha1 <- nparm+2
  pos.alpha1 <- nparm+3
  pos.0.alpha2 <- nparm+4
  pos.1.alpha2 <- nparm+5
  pos.alpha2 <- nparm+6
  pos.0.alpha3 <- nparm+7
  pos.1.alpha3 <- nparm+8
  pos.alpha3 <- nparm + 9
  
  # Table with estimates of average potential outcomes
  
  temp <- avg.po(pos = pos.0.alpha1, reg.estimates)
  row1 <- data.frame(estimand = "Average potential outcome z=0", alpha = alpha[1],
                     estimate = temp[1], variance = temp[2], 
                     CI_lb = temp[1] - qnorm(0.975)*sqrt(temp[2]),
                     CI_ub = temp[1] + qnorm(0.975)*sqrt(temp[2]))
  temp <- avg.po(pos = pos.1.alpha1, reg.estimates)
  row2 <- data.frame(estimand = "Average potential outcome z=1", alpha = alpha[1],
                     estimate = temp[1], variance = temp[2], 
                     CI_lb = temp[1] - qnorm(0.975)*sqrt(temp[2]),
                     CI_ub = temp[1] + qnorm(0.975)*sqrt(temp[2]))
  temp <- avg.po(pos = pos.alpha1, reg.estimates)
  row3 <- data.frame(estimand = "Marginal average potential outcome", alpha = alpha[1],
                     estimate = temp[1], variance = temp[2], 
                     CI_lb = temp[1] - qnorm(0.975)*sqrt(temp[2]),
                     CI_ub = temp[1] + qnorm(0.975)*sqrt(temp[2]))
  
  temp <- avg.po(pos = pos.0.alpha2, reg.estimates)
  row4 <- data.frame(estimand = "Average potential outcome z=0", alpha = alpha[2],
                     estimate = temp[1], variance = temp[2], 
                     CI_lb = temp[1] - qnorm(0.975)*sqrt(temp[2]),
                     CI_ub = temp[1] + qnorm(0.975)*sqrt(temp[2]))
  temp <- avg.po(pos = pos.1.alpha2, reg.estimates)
  row5<- data.frame(estimand = "Average potential outcome z=1", alpha = alpha[2],
                    estimate = temp[1], variance = temp[2], 
                    CI_lb = temp[1] - qnorm(0.975)*sqrt(temp[2]),
                    CI_ub = temp[1] + qnorm(0.975)*sqrt(temp[2]))
  temp <- avg.po(pos = pos.alpha2, reg.estimates)
  row6 <- data.frame(estimand = "Marginal average potential outcome", alpha = alpha[2],
                     estimate = temp[1], variance = temp[2], 
                     CI_lb = temp[1] - qnorm(0.975)*sqrt(temp[2]),
                     CI_ub = temp[1] + qnorm(0.975)*sqrt(temp[2]))
  
  temp <- avg.po(pos = pos.0.alpha3, reg.estimates)
  row7 <- data.frame(estimand = "Average potential outcome z=0", alpha = alpha[3],
                     estimate = temp[1], variance = temp[2], 
                     CI_lb = temp[1] - qnorm(0.975)*sqrt(temp[2]),
                     CI_ub = temp[1] + qnorm(0.975)*sqrt(temp[2]))
  temp <- avg.po(pos = pos.1.alpha3, reg.estimates)
  row8 <- data.frame(estimand = "Average potential outcome z=1", alpha = alpha[3],
                     estimate = temp[1], variance = temp[2], 
                     CI_lb = temp[1] - qnorm(0.975)*sqrt(temp[2]),
                     CI_ub = temp[1] + qnorm(0.975)*sqrt(temp[2]))
  temp <- avg.po(pos = pos.alpha3, reg.estimates)
  row9 <- data.frame(estimand = "Marginal average potential outcome", alpha = alpha[3],
                     estimate = temp[1], variance = temp[2], 
                     CI_lb = temp[1] - qnorm(0.975)*sqrt(temp[2]),
                     CI_ub = temp[1] + qnorm(0.975)*sqrt(temp[2]))
  
  table1 <- rbind(row1, row2, row3, row4, row5, row6, row7, row8, row9)
  
  # Table with estimates of relevant contrasts
  # Direct effect
  
  temp <- contrast(pos.1.alpha1, pos.0.alpha1, reg.estimates)
  row1 <- data.frame(estimand = "Direct effect", alpha0 = alpha[1], alpha1 = alpha[1],
                     estimate = temp[1], variance = temp[2], 
                     CI_lb = temp[1] - qnorm(0.975)*sqrt(temp[2]),
                     CI_ub = temp[1] + qnorm(0.975)*sqrt(temp[2]))
  temp <- contrast(pos.1.alpha2, pos.0.alpha2, reg.estimates)
  row2 <- data.frame(estimand = "Direct effect", alpha0 = alpha[2], alpha1 = alpha[2],
                     estimate = temp[1], variance = temp[2], 
                     CI_lb = temp[1] - qnorm(0.975)*sqrt(temp[2]),
                     CI_ub = temp[1] + qnorm(0.975)*sqrt(temp[2]))
  temp <- contrast(pos.1.alpha3, pos.0.alpha3, reg.estimates)
  row3 <- data.frame(estimand = "Direct effect", alpha0 = alpha[3], alpha1 = alpha[3],
                     estimate = temp[1], variance = temp[2], 
                     CI_lb = temp[1] - qnorm(0.975)*sqrt(temp[2]),
                     CI_ub = temp[1] + qnorm(0.975)*sqrt(temp[2]))
  
  temp <- contrast(pos.0.alpha2, pos.0.alpha1, reg.estimates)
  row4 <- data.frame(estimand = "Indirect effect", alpha0 = alpha[1], alpha1 = alpha[2],
                     estimate = temp[1], variance = temp[2], 
                     CI_lb = temp[1] - qnorm(0.975)*sqrt(temp[2]),
                     CI_ub = temp[1] + qnorm(0.975)*sqrt(temp[2]))
  temp <- contrast(pos.0.alpha3, pos.0.alpha1, reg.estimates)
  row5 <- data.frame(estimand = "Indirect effect", alpha0 = alpha[1], alpha1 = alpha[3],
                     estimate = temp[1], variance = temp[2], 
                     CI_lb = temp[1] - qnorm(0.975)*sqrt(temp[2]),
                     CI_ub = temp[1] + qnorm(0.975)*sqrt(temp[2]))
  temp <- contrast(pos.0.alpha3, pos.0.alpha2, reg.estimates)
  row6 <- data.frame(estimand = "Indirect effect", alpha0 = alpha[2], alpha1 = alpha[3],
                     estimate = temp[1], variance = temp[2], 
                     CI_lb = temp[1] - qnorm(0.975)*sqrt(temp[2]),
                     CI_ub = temp[1] + qnorm(0.975)*sqrt(temp[2]))
  
  temp <- contrast(pos.1.alpha2, pos.0.alpha1, reg.estimates)
  row7 <- data.frame(estimand = "Total effect", alpha0 = alpha[1], alpha1 = alpha[2],
                     estimate = temp[1], variance = temp[2], 
                     CI_lb = temp[1] - qnorm(0.975)*sqrt(temp[2]),
                     CI_ub = temp[1] + qnorm(0.975)*sqrt(temp[2]))
  temp <- contrast(pos.1.alpha3, pos.0.alpha1, reg.estimates)
  row8 <- data.frame(estimand = "Total effect", alpha0 = alpha[1], alpha1 = alpha[3],
                     estimate = temp[1], variance = temp[2], 
                     CI_lb = temp[1] - qnorm(0.975)*sqrt(temp[2]),
                     CI_ub = temp[1] + qnorm(0.975)*sqrt(temp[2]))
  temp <- contrast(pos.1.alpha3, pos.0.alpha2, reg.estimates)
  row9 <- data.frame(estimand = "Total effect", alpha0 = alpha[2], alpha1 = alpha[3],
                     estimate = temp[1], variance = temp[2], 
                     CI_lb = temp[1] - qnorm(0.975)*sqrt(temp[2]),
                     CI_ub = temp[1] + qnorm(0.975)*sqrt(temp[2]))
  
  temp <- contrast(pos.alpha2, pos.alpha1, reg.estimates)
  row10 <- data.frame(estimand = "Overall effect", alpha0 = alpha[1], alpha1 = alpha[2],
                      estimate = temp[1], variance = temp[2], 
                      CI_lb = temp[1] - qnorm(0.975)*sqrt(temp[2]),
                      CI_ub = temp[1] + qnorm(0.975)*sqrt(temp[2]))
  temp <- contrast(pos.alpha3, pos.alpha1, reg.estimates)
  row11 <- data.frame(estimand = "Overall effect", alpha0 = alpha[1], alpha1 = alpha[3],
                      estimate = temp[1], variance = temp[2], 
                      CI_lb = temp[1] - qnorm(0.975)*sqrt(temp[2]),
                      CI_ub = temp[1] + qnorm(0.975)*sqrt(temp[2]))
  temp <- contrast(pos.alpha3, pos.alpha2, reg.estimates)
  row12 <- data.frame(estimand = "Overall effect", alpha0 = alpha[2], alpha1 = alpha[3],
                      estimate = temp[1], variance = temp[2], 
                      CI_lb = temp[1] - qnorm(0.975)*sqrt(temp[2]),
                      CI_ub = temp[1] + qnorm(0.975)*sqrt(temp[2]))
  
  table2 <- rbind(row1, row2, row3, row4, row5, row6, row7, row8, row9,
                  row10, row11, row12)
  
  return(list(table1, table2))
  
}


# Initial values for DR-BC estimator of the avg potential outcome under treatment
# coverage alpha ----------------------------------------

initial_values_drbc <- function(data, propensity_formula, outcome_formula, 
                                counterfactual, alpha){
  
  mod.ps <- glmer(propensity_formula,
                  data = data,
                  family = binomial(link = 'logit'))
  
  gamma.hat <- as.numeric(fixef(mod.ps))
  var.b <- as.numeric(VarCorr(mod.ps))
  design.matrix <- model.matrix(mod.ps)
  
  Z <- data$Z
  
  f <- sapply(data$id, FUN = int.GH, exposure.model = mod.ps, Z=Z)
  
  mod <- lm(outcome_formula,
            data = data)
  beta.hat <- as.numeric(coef(mod))
  
  Xm <- model.matrix(object = outcome_formula, data = data)
  
  Xm_df <- data.frame(Xm)
  Xm_df$id <- data$id
  Xm_df$component.id <-  data$component.id
  Xm_df <- left_join(Xm_df, counterfactual, by = "id")
  
  Xm1 <- Xm_df %>% mutate(prop_treated = counterfactual_prop_treated,
                          Z=1) %>%
    mutate(Z.prop_treated = Z*prop_treated) %>% 
    dplyr::select(-c(id, component.id, counterfactual_prop_treated, k_treated, degree)) %>% 
    as.matrix()
  colnames(Xm1) <- colnames(Xm)
  
  Xm0 <- Xm_df %>% mutate(prop_treated = counterfactual_prop_treated,
                          Z=0) %>%
    mutate(Z.prop_treated = Z*prop_treated) %>% 
    dplyr::select(-c(id, component.id, counterfactual_prop_treated, k_treated, degree)) %>% 
    as.matrix()
  colnames(Xm0) <- colnames(Xm)
  
  Xm_df <- Xm_df %>% mutate(weights1 = dbinom(x = k_treated, size = degree, 
                                              prob = alpha[1]),
                            weights2 = dbinom(x = k_treated, size = degree, 
                                              prob = alpha[2]),
                            weights3 = dbinom(x = k_treated, size = degree, 
                                              prob = alpha[3]))
  pred <- Xm %*% beta.hat
  pred.0 <- Xm0 %*% beta.hat
  pred.1 <- Xm1 %*% beta.hat
  Xm_df <- Xm_df %>% add_column(pred.0, pred.1)
  temp <- Xm_df %>% group_by(id) %>% dplyr::summarise(reg.0.alpha1 = sum(weights1*pred.0),
                                                      reg.1.alpha1 = sum(weights1*pred.1),
                                                      reg.alpha1 = alpha[1]*reg.1.alpha1 + (1-alpha[1])*reg.0.alpha1,
                                                      reg.0.alpha2 = sum(weights2*pred.0),
                                                      reg.1.alpha2 = sum(weights2*pred.1),
                                                      reg.alpha2 = alpha[2]*reg.1.alpha2 + (1-alpha[2])*reg.0.alpha2,
                                                      reg.0.alpha3 = sum(weights3*pred.0),
                                                      reg.1.alpha3 = sum(weights3*pred.1),
                                                      reg.alpha3 = alpha[3]*reg.1.alpha3 + (1-alpha[3])*reg.0.alpha3) 
  
  num1 <- (alpha[1]^data$k_treated)*((1-alpha[1])^(data$degree-data$k_treated))
  num2 <- (alpha[2]^data$k_treated)*((1-alpha[2])^(data$degree-data$k_treated))
  num3 <- (alpha[3]^data$k_treated)*((1-alpha[3])^(data$degree-data$k_treated))
  
  temp$bc.corrected.residual.0.alpha1 <- (f)^{-1}*(data$Y - pred)*as.numeric(data$Z==0)*num1
  temp$bc.corrected.residual.1.alpha1 <- (f)^{-1}*(data$Y - pred)*as.numeric(data$Z==1)*num1
  temp$bc.corrected.residual.alpha1 <- alpha[1]*temp$bc.corrected.residual.1.alpha1 + (1-alpha[1])*temp$bc.corrected.residual.0.alpha1
  temp$bc.corrected.residual.0.alpha2 <- (f)^{-1}*(data$Y - pred)*as.numeric(data$Z==0)*num2
  temp$bc.corrected.residual.1.alpha2 <- (f)^{-1}*(data$Y - pred)*as.numeric(data$Z==1)*num2
  temp$bc.corrected.residual.alpha2 <- alpha[2]*temp$bc.corrected.residual.1.alpha2 + (1-alpha[2])*temp$bc.corrected.residual.0.alpha2
  temp$bc.corrected.residual.0.alpha3 <- (f)^{-1}*(data$Y - pred)*as.numeric(data$Z==0)*num3
  temp$bc.corrected.residual.1.alpha3 <- (f)^{-1}*(data$Y - pred)*as.numeric(data$Z==1)*num3
  temp$bc.corrected.residual.alpha3 <- alpha[3]*temp$bc.corrected.residual.1.alpha3 + (1-alpha[3])*temp$bc.corrected.residual.0.alpha3
  
  
  temp <- temp %>% left_join(., data %>% dplyr::select(id, component.id), by = "id")
  
  component_summaries_drbc <- temp %>% group_by(component.id) %>%
    dplyr::summarise(drbc.0.alpha1.component = mean(reg.0.alpha1 + bc.corrected.residual.0.alpha1),
                     drbc.1.alpha1.component = mean(reg.1.alpha1+ bc.corrected.residual.1.alpha1),
                     drbc.alpha1.component = mean(reg.alpha1+ bc.corrected.residual.alpha1),
                     drbc.0.alpha2.component = mean(reg.0.alpha2 + bc.corrected.residual.0.alpha2),
                     drbc.1.alpha2.component = mean(reg.1.alpha2+ bc.corrected.residual.1.alpha2),
                     drbc.alpha2.component = mean(reg.alpha2+ bc.corrected.residual.alpha2),
                     drbc.0.alpha3.component = mean(reg.0.alpha3 + bc.corrected.residual.0.alpha3),
                     drbc.1.alpha3.component = mean(reg.1.alpha3 + bc.corrected.residual.1.alpha3),
                     drbc.alpha3.component = mean(reg.alpha3+ bc.corrected.residual.alpha3))
  
  drbc.0.alpha1 <- mean(component_summaries_drbc$drbc.0.alpha1.component)
  drbc.1.alpha1 <- mean(component_summaries_drbc$drbc.1.alpha1.component)
  drbc.alpha1 <- mean(component_summaries_drbc$drbc.alpha1.component)
  drbc.0.alpha2 <- mean(component_summaries_drbc$drbc.0.alpha2.component)
  drbc.1.alpha2 <- mean(component_summaries_drbc$drbc.1.alpha2.component)
  drbc.alpha2 <- mean(component_summaries_drbc$drbc.alpha2.component)
  drbc.0.alpha3 <- mean(component_summaries_drbc$drbc.0.alpha3.component)
  drbc.1.alpha3 <- mean(component_summaries_drbc$drbc.1.alpha3.component)
  drbc.alpha3 <- mean(component_summaries_drbc$drbc.alpha3.component)
  
  return(c(gamma.hat, var.b, beta.hat, drbc.0.alpha1, drbc.1.alpha1, 
         drbc.alpha1, drbc.0.alpha2, drbc.1.alpha2, drbc.alpha2,
         drbc.0.alpha3, drbc.1.alpha3, drbc.alpha3))
  
}

# Estimating function for the DRBC estimator ------------------------------

drbc_estfun <- function(data, models, alpha, counterfactual_prop_treated, N, m){ 
  
  id <- data$id
  Z <- data$Z
  Y <- data$Y
  X.N <- data$X.N
  Z.N <- data$Z.N
  k_treated <- data$k_treated
  degree <- data$degree
  Xe <- grab_design_matrix(data = data,
                           rhs_formula = grab_fixed_formula(models$e))
  
  Xm <- grab_design_matrix(data = data,
                           rhs_formula = grab_fixed_formula(models$m))
  
  # Positions for parameters
  gamma_pos <- 1:(ncol(Xe))
  varb_pos <- length(gamma_pos)+1
  beta_pos <- (varb_pos + 1):(varb_pos+ncol(Xm))
  
  # Create data sets for prediction 
  Xm_df <- data.frame(Xm)
  Xm_df$id <- data$id
  Xm_df <- left_join(Xm_df, counterfactual_prop_treated, by = "id")
  
  Xm1 <- Xm_df %>% mutate(prop_treated = counterfactual_prop_treated,
                          Z=1) %>%
    mutate(Z.prop_treated = Z*prop_treated) %>% 
    dplyr::select(-c(id, counterfactual_prop_treated, k_treated, degree)) %>% 
    as.matrix()
  colnames(Xm1) <- colnames(Xm)
  
  Xm0 <- Xm_df %>% mutate(prop_treated = counterfactual_prop_treated,
                          Z=0) %>%
    mutate(Z.prop_treated = Z*prop_treated) %>% 
    dplyr::select(-c(id, counterfactual_prop_treated, k_treated, degree)) %>% 
    as.matrix()
  colnames(Xm0) <- colnames(Xm)
  
  
  # Prepare data for IPW estimation
  
  g <- function(X.N, Z.N, b, coef){
    
    # Predict neighborhood level propensity scores
    p.n <- inv.logit(X.N %*% coef + b)
    
    return(prod(p.n^Z.N*(1-p.n)^(1-Z.N)))
  }
  
  int.GH <- function(X.N, Z.N, coef, var.b){
    w <- gauss.quad(10, kind="hermite")
    
    gx <- sapply(w$nodes*sqrt(2*var.b), FUN = g, X.N = X.N,
                 Z.N = Z.N, coef=coef)
    return(sum(w$weights*gx/sqrt(pi)))
  }
  
  
  # Grab estimating functions
  
  e_scores <- grab_psiFUN(models$e, data)
  m_scores <- grab_psiFUN(models$m, data)
  
  Xm_df <- Xm_df %>% mutate(weights1 = dbinom(x = k_treated, size = degree, 
                                              prob = alpha[1]),
                            weights2 = dbinom(x = k_treated, size = degree, 
                                              prob = alpha[2]),
                            weights3 = dbinom(x = k_treated, size = degree, 
                                              prob = alpha[3]))
  
  num1 <- (alpha[1]^k_treated)*((1-alpha[1])^(degree-k_treated))
  num2 <- (alpha[2]^k_treated)*((1-alpha[2])^(degree-k_treated))
  num3 <- (alpha[3]^k_treated)*((1-alpha[3])^(degree-k_treated))
  
  function(theta) {
    gamma <- theta[gamma_pos]
    var.b <- theta[varb_pos]
    e <- map2_dbl(X.N, Z.N, int.GH, coef = gamma, var.b = var.b)
    
    beta <- theta[beta_pos]
    pred.0 <- Xm0 %*% beta
    pred.1 <- Xm1 %*% beta
    Xm_df <- Xm_df %>% add_column(pred.0, pred.1)
    temp <- Xm_df %>% group_by(id) %>% dplyr::summarise(reg.0.alpha1 = sum(weights1*pred.0),
                                                        reg.1.alpha1 = sum(weights1*pred.1),
                                                        reg.alpha1 = alpha[1]*reg.1.alpha1 + (1-alpha[1])*reg.0.alpha1,
                                                        reg.0.alpha2 = sum(weights2*pred.0),
                                                        reg.1.alpha2 = sum(weights2*pred.1),
                                                        reg.alpha2 = alpha[2]*reg.1.alpha2 + (1-alpha[2])*reg.0.alpha2,
                                                        reg.0.alpha3 = sum(weights3*pred.0),
                                                        reg.1.alpha3 = sum(weights3*pred.1),
                                                        reg.alpha3 = alpha[3]*reg.1.alpha3 + (1-alpha[3])*reg.0.alpha3) 
    reg.0.alpha1 <- mean(temp$reg.0.alpha1)
    reg.1.alpha1 <- mean(temp$reg.1.alpha1)
    reg.0.alpha2 <- mean(temp$reg.0.alpha2)
    reg.1.alpha2 <- mean(temp$reg.1.alpha2)
    reg.0.alpha3 <- mean(temp$reg.0.alpha3)
    reg.1.alpha3 <- mean(temp$reg.1.alpha3)
    
    pred <- Xm %*% beta
    
    drbc.0.alpha1 <- reg.0.alpha1 + mean((e)^{-1}*(Y - pred)*as.numeric(Z==0)*num1)
    drbc.1.alpha1 <- reg.1.alpha1 + mean((e)^{-1}*(Y - pred)*as.numeric(Z==1)*num1)
    drbc.alpha1 <- alpha[1]*drbc.1.alpha1 + (1-alpha[1])*drbc.0.alpha1
    drbc.0.alpha2 <- reg.0.alpha2 + mean((e)^{-1}*(Y - pred)*as.numeric(Z==0)*num2)
    drbc.1.alpha2 <- reg.1.alpha2 + mean((e)^{-1}*(Y - pred)*as.numeric(Z==1)*num2)
    drbc.alpha2 <- alpha[2]*drbc.1.alpha2 + (1-alpha[2])*drbc.0.alpha2
    drbc.0.alpha3 <- reg.0.alpha3 + mean((e)^{-1}*(Y - pred)*as.numeric(Z==0)*num3)
    drbc.1.alpha3 <- reg.1.alpha3 + mean((e)^{-1}*(Y - pred)*as.numeric(Z==1)*num3)
    drbc.alpha3 <- alpha[3]*drbc.1.alpha3 + (1-alpha[3])*drbc.0.alpha3
    
    return(c(e_scores(theta[c(gamma_pos, varb_pos)]),
             m_scores(theta[beta_pos]),
             drbc.0.alpha1 - theta[length(theta)-8],
             drbc.1.alpha1 - theta[length(theta)-7],
             drbc.alpha1 - theta[length(theta)-6],
             drbc.0.alpha2 - theta[length(theta)-5],
             drbc.1.alpha2 - theta[length(theta)-4],
             drbc.alpha2 - theta[length(theta)-3],
             drbc.0.alpha3 - theta[length(theta)-2],
             drbc.1.alpha3 - theta[length(theta)-1],
             drbc.alpha3 - theta[length(theta)]))
  }
}

estimate_drbc <- function(data, propensity_formula, outcome_formula, init, alpha, 
                          counterfactual_prop_treated, m){
  e_model <- glmer(propensity_formula,
                   data = data,
                   family = binomial(link = 'logit'))
  m_model <- glm(outcome_formula,
                 data = data, family = "gaussian")
  models <- list(e = e_model, m = m_model)
  N <- nrow(data)
  m_estimate(estFUN = drbc_estfun, data = data, units="component.id",
             root_control = setup_root_control(start = init),
             outer_args = list(models = models,
                               alpha = alpha,
                               counterfactual_prop_treated = counterfactual_prop_treated,
                               m=m, 
                               N=N))
}

drbc_estimation <- function(data, propensity_formula, outcome_formula, alpha){
  
  #create a copy of the data set to alter it
  
  Z <- data$Z
  
  mod.ps <- glmer(propensity_formula,
                  data = data,
                  family = binomial(link = 'logit'))
  design.matrix <- model.matrix(mod.ps)
  
  data$X.N  <- sapply(data$id, FUN = x.alters, 
                      design.matrix=design.matrix)
  data$Z.N  <- sapply(data$id, FUN = z.alters, Z=Z)
  
  Xm <- model.matrix(object = outcome_formula, data = data)
  Xe <- model.matrix(object = propensity_formula, data = data)
  
  counterfactual <- data %>% uncount(degree+1) %>% 
    dplyr::select(-k_treated) %>%
    mutate(prop_treated = 0) %>% group_by(id) %>% 
    mutate_at(vars(prop_treated), ~. + 0:(first(degree)))  %>%
    rename(k_treated = prop_treated) %>% 
    mutate(prop_treated = ifelse(degree==0, 0, k_treated/degree)) %>%
    dplyr::select(id, prop_treated, k_treated, degree) %>% 
    rename(counterfactual_prop_treated = prop_treated)
  
  init_drbc <- initial_values_drbc(data=data, propensity_formula=propensity_formula,
                                   outcome_formula=outcome_formula, 
                                  counterfactual=counterfactual, alpha=alpha)
  
  # Call the M-estimation routine with the above initial values
  # This might take some time
  
  drbc.estimates <- estimate_drbc(data = data, 
                    propensity_formula = propensity_formula,
                    outcome_formula = outcome_formula,
                    alpha = alpha,
                    init = init_drbc,
                    counterfactual_prop_treated = counterfactual,
                    m=m)
  
  # Retrieve the quantities of interest and calculate causal contrasts
  
  nparm <- ncol(Xm)+ncol(Xe)
  pos.0.alpha1 <- nparm+1
  pos.1.alpha1 <- nparm+2
  pos.alpha1 <- nparm+3
  pos.0.alpha2 <- nparm+4
  pos.1.alpha2 <- nparm+5
  pos.alpha2 <- nparm+6
  pos.0.alpha3 <- nparm+7
  pos.1.alpha3 <- nparm+8
  pos.alpha3 <- nparm + 9
  
  # Table with estimates of average potential outcomes
  
  temp <- avg.po(pos = pos.0.alpha1, drbc.estimates)
  row1 <- data.frame(estimand = "Average potential outcome z=0", alpha = alpha[1],
                     estimate = temp[1], variance = temp[2], 
                     CI_lb = temp[1] - qnorm(0.975)*sqrt(temp[2]),
                     CI_ub = temp[1] + qnorm(0.975)*sqrt(temp[2]))
  temp <- avg.po(pos = pos.1.alpha1, drbc.estimates)
  row2 <- data.frame(estimand = "Average potential outcome z=1", alpha = alpha[1],
                     estimate = temp[1], variance = temp[2], 
                     CI_lb = temp[1] - qnorm(0.975)*sqrt(temp[2]),
                     CI_ub = temp[1] + qnorm(0.975)*sqrt(temp[2]))
  temp <- avg.po(pos = pos.alpha1, drbc.estimates)
  row3 <- data.frame(estimand = "Marginal average potential outcome", alpha = alpha[1],
                     estimate = temp[1], variance = temp[2], 
                     CI_lb = temp[1] - qnorm(0.975)*sqrt(temp[2]),
                     CI_ub = temp[1] + qnorm(0.975)*sqrt(temp[2]))
  
  temp <- avg.po(pos = pos.0.alpha2, drbc.estimates)
  row4 <- data.frame(estimand = "Average potential outcome z=0", alpha = alpha[2],
                     estimate = temp[1], variance = temp[2], 
                     CI_lb = temp[1] - qnorm(0.975)*sqrt(temp[2]),
                     CI_ub = temp[1] + qnorm(0.975)*sqrt(temp[2]))
  temp <- avg.po(pos = pos.1.alpha2, drbc.estimates)
  row5<- data.frame(estimand = "Average potential outcome z=1", alpha = alpha[2],
                    estimate = temp[1], variance = temp[2], 
                    CI_lb = temp[1] - qnorm(0.975)*sqrt(temp[2]),
                    CI_ub = temp[1] + qnorm(0.975)*sqrt(temp[2]))
  temp <- avg.po(pos = pos.alpha2, drbc.estimates)
  row6 <- data.frame(estimand = "Marginal average potential outcome", alpha = alpha[2],
                     estimate = temp[1], variance = temp[2], 
                     CI_lb = temp[1] - qnorm(0.975)*sqrt(temp[2]),
                     CI_ub = temp[1] + qnorm(0.975)*sqrt(temp[2]))
  
  temp <- avg.po(pos = pos.0.alpha3, drbc.estimates)
  row7 <- data.frame(estimand = "Average potential outcome z=0", alpha = alpha[3],
                     estimate = temp[1], variance = temp[2], 
                     CI_lb = temp[1] - qnorm(0.975)*sqrt(temp[2]),
                     CI_ub = temp[1] + qnorm(0.975)*sqrt(temp[2]))
  temp <- avg.po(pos = pos.1.alpha3, drbc.estimates)
  row8 <- data.frame(estimand = "Average potential outcome z=1", alpha = alpha[3],
                     estimate = temp[1], variance = temp[2], 
                     CI_lb = temp[1] - qnorm(0.975)*sqrt(temp[2]),
                     CI_ub = temp[1] + qnorm(0.975)*sqrt(temp[2]))
  temp <- avg.po(pos = pos.alpha3, drbc.estimates)
  row9 <- data.frame(estimand = "Marginal average potential outcome", alpha = alpha[3],
                     estimate = temp[1], variance = temp[2], 
                     CI_lb = temp[1] - qnorm(0.975)*sqrt(temp[2]),
                     CI_ub = temp[1] + qnorm(0.975)*sqrt(temp[2]))
  
  table1 <- rbind(row1, row2, row3, row4, row5, row6, row7, row8, row9)
  
  # Table with estimates of relevant contrasts
  # Direct effect
  
  temp <- contrast(pos.1.alpha1, pos.0.alpha1, drbc.estimates)
  row1 <- data.frame(estimand = "Direct effect", alpha0 = alpha[1], alpha1 = alpha[1],
                     estimate = temp[1], variance = temp[2], 
                     CI_lb = temp[1] - qnorm(0.975)*sqrt(temp[2]),
                     CI_ub = temp[1] + qnorm(0.975)*sqrt(temp[2]))
  temp <- contrast(pos.1.alpha2, pos.0.alpha2, drbc.estimates)
  row2 <- data.frame(estimand = "Direct effect", alpha0 = alpha[2], alpha1 = alpha[2],
                     estimate = temp[1], variance = temp[2], 
                     CI_lb = temp[1] - qnorm(0.975)*sqrt(temp[2]),
                     CI_ub = temp[1] + qnorm(0.975)*sqrt(temp[2]))
  temp <- contrast(pos.1.alpha3, pos.0.alpha3, drbc.estimates)
  row3 <- data.frame(estimand = "Direct effect", alpha0 = alpha[3], alpha1 = alpha[3],
                     estimate = temp[1], variance = temp[2], 
                     CI_lb = temp[1] - qnorm(0.975)*sqrt(temp[2]),
                     CI_ub = temp[1] + qnorm(0.975)*sqrt(temp[2]))
  
  temp <- contrast(pos.0.alpha2, pos.0.alpha1, drbc.estimates)
  row4 <- data.frame(estimand = "Indirect effect", alpha0 = alpha[1], alpha1 = alpha[2],
                     estimate = temp[1], variance = temp[2], 
                     CI_lb = temp[1] - qnorm(0.975)*sqrt(temp[2]),
                     CI_ub = temp[1] + qnorm(0.975)*sqrt(temp[2]))
  temp <- contrast(pos.0.alpha3, pos.0.alpha1, drbc.estimates)
  row5 <- data.frame(estimand = "Indirect effect", alpha0 = alpha[1], alpha1 = alpha[3],
                     estimate = temp[1], variance = temp[2], 
                     CI_lb = temp[1] - qnorm(0.975)*sqrt(temp[2]),
                     CI_ub = temp[1] + qnorm(0.975)*sqrt(temp[2]))
  temp <- contrast(pos.0.alpha3, pos.0.alpha2, drbc.estimates)
  row6 <- data.frame(estimand = "Indirect effect", alpha0 = alpha[2], alpha1 = alpha[3],
                     estimate = temp[1], variance = temp[2], 
                     CI_lb = temp[1] - qnorm(0.975)*sqrt(temp[2]),
                     CI_ub = temp[1] + qnorm(0.975)*sqrt(temp[2]))
  
  temp <- contrast(pos.1.alpha2, pos.0.alpha1, drbc.estimates)
  row7 <- data.frame(estimand = "Total effect", alpha0 = alpha[1], alpha1 = alpha[2],
                     estimate = temp[1], variance = temp[2], 
                     CI_lb = temp[1] - qnorm(0.975)*sqrt(temp[2]),
                     CI_ub = temp[1] + qnorm(0.975)*sqrt(temp[2]))
  temp <- contrast(pos.1.alpha3, pos.0.alpha1, drbc.estimates)
  row8 <- data.frame(estimand = "Total effect", alpha0 = alpha[1], alpha1 = alpha[3],
                     estimate = temp[1], variance = temp[2], 
                     CI_lb = temp[1] - qnorm(0.975)*sqrt(temp[2]),
                     CI_ub = temp[1] + qnorm(0.975)*sqrt(temp[2]))
  temp <- contrast(pos.1.alpha3, pos.0.alpha2, drbc.estimates)
  row9 <- data.frame(estimand = "Total effect", alpha0 = alpha[2], alpha1 = alpha[3],
                     estimate = temp[1], variance = temp[2], 
                     CI_lb = temp[1] - qnorm(0.975)*sqrt(temp[2]),
                     CI_ub = temp[1] + qnorm(0.975)*sqrt(temp[2]))
  
  temp <- contrast(pos.alpha2, pos.alpha1, drbc.estimates)
  row10 <- data.frame(estimand = "Overall effect", alpha0 = alpha[1], alpha1 = alpha[2],
                      estimate = temp[1], variance = temp[2], 
                      CI_lb = temp[1] - qnorm(0.975)*sqrt(temp[2]),
                      CI_ub = temp[1] + qnorm(0.975)*sqrt(temp[2]))
  temp <- contrast(pos.alpha3, pos.alpha1, drbc.estimates)
  row11 <- data.frame(estimand = "Overall effect", alpha0 = alpha[1], alpha1 = alpha[3],
                      estimate = temp[1], variance = temp[2], 
                      CI_lb = temp[1] - qnorm(0.975)*sqrt(temp[2]),
                      CI_ub = temp[1] + qnorm(0.975)*sqrt(temp[2]))
  temp <- contrast(pos.alpha3, pos.alpha2, drbc.estimates)
  row12 <- data.frame(estimand = "Overall effect", alpha0 = alpha[2], alpha1 = alpha[3],
                      estimate = temp[1], variance = temp[2], 
                      CI_lb = temp[1] - qnorm(0.975)*sqrt(temp[2]),
                      CI_ub = temp[1] + qnorm(0.975)*sqrt(temp[2]))
  
  table2 <- rbind(row1, row2, row3, row4, row5, row6, row7, row8, row9,
                  row10, row11, row12)
  
  return(list(table1, table2))
  
  
}

# Estimating functions for the DRWLS estimator ------------------------------

g.drwls <- function(b, i, dat, exposure.model, z, Z){
  
  gamma.hat <- as.numeric(fixef(exposure.model))
  
  var.b <- as.numeric(VarCorr(exposure.model))
  
  design.matrix <- model.matrix(exposure.model)
  
  # i-th row of the adjacency matrix
  x <- A[i,]
  X.N <- rbind(design.matrix[i,], design.matrix[x==1,])
  
  # Predict neighborhood level propensity scores
  temp <- inv.logit(as.numeric(X.N %*% as.matrix(gamma.hat)) + rep(b, nrow(X.N)))
  
  p <- temp[1]
  p.N <- temp[-1]
  
  # Neighborhood treatments
  Z.N <- dat$Z[x==1]
  
  return(p^z*(1-p)^(1-z)*prod(p.N^Z.N*(1-p.N)^(1-Z.N)))
}


int.GH.drwls <- function(i, data, exposure.model, z, Z){
  var.b <- as.numeric(VarCorr(exposure.model))
  
  w <- gauss.quad(10, kind="hermite")
  
  gx <- sapply(w$nodes*sqrt(2*var.b), FUN = g.drwls, i = i, z=z, Z=Z,
               exposure.model = exposure.model,
               dat=data)
  return(sum(w$weights*gx/sqrt(pi)))
}



drwls_estfun <- function(data, models, alpha, outcome_formula, 
                         counterfactual_prop_treated,
                         m, N){ 
  
  id <- data$id
  Z <- data$Z
  Y <- data$Y
  X.N <- data$X.N
  Z.N <- data$Z.N
  k_treated <- data$k_treated
  degree <- data$degree
  Xe <- grab_design_matrix(data = data,
                           rhs_formula = grab_fixed_formula(models$e))
  
  # Data for estimating function for the WLS estimators 
  Xm <- model.matrix(object = outcome_formula, data = data)
  Y <- model.response(model.frame(outcome_formula, data = data))
  
  
  # Positions for parameters 
  gamma_pos <- 1:(ncol(Xe))
  varb_pos <- max(gamma_pos)+1
  beta0_alpha1_pos <- (varb_pos + 1):(varb_pos+ncol(Xm)) 
  beta1_alpha1_pos <- (max(beta0_alpha1_pos) + 1):(max(beta0_alpha1_pos) + ncol(Xm))
  beta0_alpha2_pos <- (max(beta1_alpha1_pos) + 1):(max(beta1_alpha1_pos) + ncol(Xm))
  beta1_alpha2_pos <- (max(beta0_alpha2_pos) + 1):(max(beta0_alpha2_pos) + ncol(Xm))
  beta0_alpha3_pos <- (max(beta1_alpha2_pos) + 1):(max(beta1_alpha2_pos) + ncol(Xm))
  beta1_alpha3_pos <- (max(beta0_alpha3_pos) + 1):(max(beta0_alpha3_pos) + ncol(Xm))
  
  # Create data sets for prediction 
  
  Xm_df <- data.frame(Xm)
  Xm_df$id <- data$id
  Xm_df <- left_join(Xm_df, counterfactual_prop_treated, by = "id")
  
  Xm1 <- Xm_df %>% mutate(prop_treated = counterfactual_prop_treated) %>%
    dplyr::select(-c(id, counterfactual_prop_treated, k_treated, degree)) %>% 
    as.matrix()
  colnames(Xm1) <- colnames(Xm)
  
  Xm0 <- Xm_df %>% mutate(prop_treated = counterfactual_prop_treated) %>%
    dplyr::select(-c(id, counterfactual_prop_treated, k_treated, degree)) %>% 
    as.matrix()
  colnames(Xm0) <- colnames(Xm)
  
  
  # Prepare data for IPW estimation 
  num1 <- (alpha[1]^k_treated)*((1-alpha[1])^(degree-k_treated))
  num2 <- (alpha[2]^k_treated)*((1-alpha[2])^(degree-k_treated))
  num3 <- (alpha[3]^k_treated)*((1-alpha[3])^(degree-k_treated))
  
  g  <- function(X.N, Z.N, b, coef, z, Z){
    
    # Predict neighborhood level propensity scores
    temp <- inv.logit(X.N %*% coef + b)
    
    p <- temp[1]
    p.neigh <- temp[-1]
    Z.neigh <- Z.N[-1]
    
    return(p^z*(1-p)^(1-z)*prod(p.neigh^Z.neigh*(1-p.neigh)^(1-Z.neigh)))
  }
  
  int.GH <- function(X.N, Z.N, coef, var.b, z){
    w <- gauss.quad(10, kind="hermite")
    
    gx <- sapply(w$nodes*sqrt(2*var.b), FUN = g, X.N = X.N,
                 Z.N = Z.N, coef=coef, z=z)
    return(sum(w$weights*gx/sqrt(pi)))
  }
  
  
  # Grab estimating functions ----------------------------------------------
  
  e_scores <- grab_psiFUN(models$e, data)

  Xm_df <- Xm_df %>% mutate(weights1 = dbinom(x = k_treated, size = degree, 
                                             prob = alpha[1]),
                            weights2 = dbinom(x = k_treated, size = degree, 
                                              prob = alpha[2]),
                            weights3 = dbinom(x = k_treated, size = degree, 
                                              prob = alpha[3]))
  
  
  function(theta) {
    
    gamma <- theta[gamma_pos]
    var.b <- theta[varb_pos]
    
    # Compute propensity scores for each treatment arm 
    e0 <- map2_dbl(X.N, Z.N, int.GH, coef = gamma, var.b = var.b, z=0)
    e1 <- map2_dbl(X.N, Z.N, int.GH, coef = gamma, var.b = var.b, z=1)
    
    
    #Compute doubly robust estimator 
    beta0_alpha1 <- theta[beta0_alpha1_pos]
    beta1_alpha1 <- theta[beta1_alpha1_pos]
    beta0_alpha2 <- theta[beta0_alpha2_pos]
    beta1_alpha2 <- theta[beta1_alpha2_pos]
    beta0_alpha3 <- theta[beta0_alpha3_pos]
    beta1_alpha3 <- theta[beta1_alpha3_pos]
    pred.0.alpha1 <- Xm0 %*% beta0_alpha1
    pred.1.alpha1 <- Xm1 %*% beta1_alpha1
    pred.0.alpha2 <- Xm0 %*% beta0_alpha2
    pred.1.alpha2 <- Xm1 %*% beta1_alpha2
    pred.0.alpha3 <- Xm0 %*% beta0_alpha3 
    pred.1.alpha3 <- Xm1 %*% beta1_alpha3
    Xm_df <- Xm_df %>% add_column(pred.0.alpha1, pred.1.alpha1,
                                  pred.0.alpha2, pred.1.alpha2,
                                  pred.0.alpha3, pred.1.alpha3)
    temp <- Xm_df %>% group_by(id) %>% summarise(wls.0.alpha1.ind = sum(weights1*pred.0.alpha1),
                                                 wls.1.alpha1.ind = sum(weights1*pred.1.alpha1),
                                                 wls.alpha1.ind = alpha[1]*wls.1.alpha1.ind + (1-alpha[1])*wls.0.alpha1.ind,
                                                 wls.0.alpha2.ind = sum(weights2*pred.0.alpha2),
                                                 wls.1.alpha2.ind = sum(weights2*pred.1.alpha2),
                                                 wls.alpha2.ind = alpha[2]*wls.1.alpha2.ind + (1-alpha[2])*wls.0.alpha2.ind,
                                                 wls.0.alpha3.ind = sum(weights3*pred.0.alpha3),
                                                 wls.1.alpha3.ind = sum(weights3*pred.1.alpha3),
                                                 wls.alpha3.ind = alpha[3]*wls.1.alpha3.ind + (1-alpha[3])*wls.0.alpha3.ind) 
    drwls.0.alpha1 <- mean(temp$wls.0.alpha1.ind)
    drwls.1.alpha1 <- mean(temp$wls.1.alpha1.ind)
    drwls.alpha1 <- mean(temp$wls.alpha1.ind)
    drwls.0.alpha2 <- mean(temp$wls.0.alpha2.ind)
    drwls.1.alpha2 <- mean(temp$wls.1.alpha2.ind)
    drwls.alpha2 <- mean(temp$wls.alpha2.ind)
    drwls.0.alpha3 <- mean(temp$wls.0.alpha3.ind)
    drwls.1.alpha3 <- mean(temp$wls.1.alpha3.ind)
    drwls.alpha3 <- mean(temp$wls.alpha3.ind)
    
    
    # Variables for the EF of beta wls
    Lambda0.alpha1 <- diag(num1*(Z==0)/e0)
    Lambda1.alpha1 <- diag(num1*(Z==1)/e1)
    Lambda0.alpha2 <- diag(num2*(Z==0)/e0)
    Lambda1.alpha2 <- diag(num2*(Z==1)/e1)
    Lambda0.alpha3 <- diag(num3*(Z==0)/e0)
    Lambda1.alpha3 <- diag(num3*(Z==1)/e1)
    
    return(c(e_scores(theta[c(gamma_pos, varb_pos)]),
             t(Xm) %*% Lambda0.alpha1 %*% (Y-Xm %*% beta0_alpha1),
             t(Xm) %*% Lambda1.alpha1 %*% (Y-Xm %*% beta1_alpha1),
             t(Xm) %*% Lambda0.alpha2 %*% (Y-Xm %*% beta0_alpha2),
             t(Xm) %*% Lambda1.alpha2 %*% (Y-Xm %*% beta1_alpha2),
             t(Xm) %*% Lambda0.alpha3 %*% (Y-Xm %*% beta0_alpha3),
             t(Xm) %*% Lambda1.alpha3 %*% (Y-Xm %*% beta1_alpha3),
             drwls.0.alpha1 - theta[length(theta)-8],
             drwls.1.alpha1 - theta[length(theta)-7],
             drwls.alpha1 - theta[length(theta)-6],
             drwls.0.alpha2 - theta[length(theta)-5],
             drwls.1.alpha2 - theta[length(theta)-4],
             drwls.alpha2 - theta[length(theta)-3],
             drwls.0.alpha3 - theta[length(theta)-2],
             drwls.1.alpha3 - theta[length(theta)-1],
             drwls.alpha3 - theta[length(theta)]))
  }
}


estimate_drwls <- function(data, propensity_formula, outcome_formula, init, alpha, 
                           counterfactual_prop_treated, m){
  e_model <- glmer(propensity_formula,
                   data = data,
                   family = binomial(link = 'logit'))
  models <- list(e = e_model)
  N <- nrow(data)
  m_estimate(estFUN = drwls_estfun, data = data, units="component.id",
             root_control = setup_root_control(start = init),
             outer_args = list(models = models,
                               alpha = alpha,
                               outcome_formula = outcome_formula,
                               counterfactual_prop_treated = counterfactual_prop_treated,
                               m=m, N=N))
}


initial_values_drwls <- function(data, propensity_formula, outcome_formula, 
                                counterfactual, alpha){
  
  
  mod.ps <- glmer(propensity_formula,
                  data = data,
                  family = binomial(link = 'logit'))
  
  gamma.hat <- as.numeric(fixef(mod.ps))
  var.b <- as.numeric(VarCorr(mod.ps))
  design.matrix <- model.matrix(mod.ps)
  
  Z <- data$Z
  
  f <- sapply(data$id, FUN = int.GH, exposure.model = mod.ps, Z=Z)
  
  f.0 <- sapply(data$id, FUN = int.GH.drwls, data=data, exposure.model = mod.ps,
                            z=0, Z=Z)
  f.1 <- sapply(data$id, FUN = int.GH.drwls, data=data, exposure.model = mod.ps,
                            z=1, Z=Z)
  
  num1 <- (alpha[1]^data$k_treated)*((1-alpha[1])^(data$degree-data$k_treated))
  num2 <- (alpha[2]^data$k_treated)*((1-alpha[2])^(data$degree-data$k_treated))
  num3 <- (alpha[3]^data$k_treated)*((1-alpha[3])^(data$degree-data$k_treated))
  
  data$w0_1 <- (data$Z==0)*num1/f.0
  data$w1_1 <- (data$Z==1)*num1/f.1
  
  data$w0_2 <- (data$Z==0)*num2/f.0
  data$w1_2 <- (data$Z==1)*num2/f.1
  
  data$w0_3 <- (data$Z==0)*num3/f.0
  data$w1_3 <- (data$Z==1)*num3/f.1
  
  Xm <- model.matrix(object = outcome_formula, data = data)
  Xm_df <- data.frame(Xm)
  Xm_df$id <- data$id
  Xm_df <- left_join(Xm_df, counterfactual, by = "id")
  Xm_pred <- Xm_df %>% mutate(prop_treated = counterfactual_prop_treated) %>%
    dplyr::select(-c(id, counterfactual_prop_treated, k_treated, degree)) %>% 
    as.matrix()
  colnames(Xm_pred) <- colnames(Xm)
  Xm_df <- Xm_df %>% mutate(weights1 = dbinom(x = k_treated, size = degree, 
                                              prob = alpha[1]),
                            weights2 = dbinom(x = k_treated, size = degree, 
                                              prob = alpha[2]),
                            weights3 = dbinom(x = k_treated, size = degree, 
                                              prob = alpha[3]))
  reg.mod.0.1 <- lm(outcome_formula, data = data, weights = w0_1)
  reg.mod.1.1 <- lm(outcome_formula, data = data, weights = w1_1)
  reg.mod.0.2 <- lm(outcome_formula, data = data, weights = w0_2)
  reg.mod.1.2 <- lm(outcome_formula, data = data, weights = w1_2)
  reg.mod.0.3 <- lm(outcome_formula, data = data, weights = w0_3)
  reg.mod.1.3 <- lm(outcome_formula, data = data, weights = w1_3)
  beta0.1 <- coefficients(reg.mod.0.1)
  beta1.1 <- coefficients(reg.mod.1.1)
  beta0.2 <- coefficients(reg.mod.0.2)
  beta1.2 <- coefficients(reg.mod.1.2)
  beta0.3 <- coefficients(reg.mod.0.3)
  beta1.3 <- coefficients(reg.mod.1.3)
  
  #Compute doubly robust estimator -----------------------------------------
  pred.0.1 <- Xm_pred %*% beta0.1
  pred.1.1 <- Xm_pred %*% beta1.1
  pred.0.2 <- Xm_pred %*% beta0.2
  pred.1.2 <- Xm_pred %*% beta1.2
  pred.0.3 <- Xm_pred %*% beta0.3
  pred.1.3 <- Xm_pred %*% beta1.3
  Xm_df <- Xm_df %>% add_column(pred.0.1, pred.1.1, pred.0.2, pred.1.2, 
                                pred.0.3, pred.1.3)
  temp <- Xm_df %>% group_by(id) %>% summarise(wls.0.1.ind = sum(weights1*pred.0.1),
                                               wls.1.1.ind = sum(weights1*pred.1.1),
                                               wls.1.ind = alpha[1]*wls.1.1.ind + (1-alpha[1])*wls.0.1.ind,
                                               wls.0.2.ind = sum(weights2*pred.0.2),
                                               wls.1.2.ind = sum(weights2*pred.1.2),
                                               wls.2.ind = alpha[2]*wls.1.2.ind + (1-alpha[2])*wls.0.2.ind,
                                               wls.0.3.ind = sum(weights3*pred.0.3),
                                               wls.1.3.ind = sum(weights3*pred.1.3),
                                               wls.3.ind = alpha[3]*wls.1.3.ind + (1-alpha[3])*wls.0.3.ind) 
  
  temp <- temp %>% left_join(., data %>% dplyr::select(id, component.id), by = "id")
  
  component_summaries_drwls <- temp %>% group_by(component.id) %>%
    dplyr::summarise(drwls.0.alpha1.component = mean(wls.0.1.ind),
                     drwls.1.alpha1.component = mean(wls.1.1.ind),
                     drwls.alpha1.component = alpha[1]*drwls.1.alpha1.component + (1-alpha[1])*drwls.0.alpha1.component,
                     drwls.0.alpha2.component = mean(wls.0.2.ind),
                     drwls.1.alpha2.component = mean(wls.1.2.ind),
                     drwls.alpha2.component = alpha[2]*drwls.1.alpha2.component + (1-alpha[2])*drwls.0.alpha2.component,
                     drwls.0.alpha3.component = mean(wls.0.3.ind),
                     drwls.1.alpha3.component = mean(wls.1.3.ind),
                     drwls.alpha3.component = alpha[3]*drwls.1.alpha3.component + (1-alpha[3])*drwls.0.alpha3.component)
  
  drwls.0.alpha1 <- mean(component_summaries_drwls$drwls.0.alpha1.component)
  drwls.1.alpha1 <- mean(component_summaries_drwls$drwls.1.alpha1.component)
  drwls.alpha1 <- mean(component_summaries_drwls$drwls.alpha1.component)
  drwls.0.alpha2 <- mean(component_summaries_drwls$drwls.0.alpha2.component)
  drwls.1.alpha2 <- mean(component_summaries_drwls$drwls.1.alpha2.component)
  drwls.alpha2 <- mean(component_summaries_drwls$drwls.alpha2.component)
  drwls.0.alpha3 <- mean(component_summaries_drwls$drwls.0.alpha3.component)
  drwls.1.alpha3 <- mean(component_summaries_drwls$drwls.1.alpha3.component)
  drwls.alpha3 <- mean(component_summaries_drwls$drwls.alpha3.component)
  
  return(c(gamma.hat, var.b, 
           as.numeric(beta0.1), as.numeric(beta1.1),
           as.numeric(beta0.2), as.numeric(beta1.2),
           as.numeric(beta0.3), as.numeric(beta1.3),
           drwls.0.alpha1, drwls.1.alpha1, drwls.alpha1,
           drwls.0.alpha2, drwls.1.alpha2, drwls.alpha2,
           drwls.0.alpha3, drwls.1.alpha3, drwls.alpha3))
  
}

drwls_estimation <- function(data, propensity_formula, outcome_formula, alpha){
  
  #create a copy of the data set to alter it
  
  Z <- data$Z
  
  mod.ps <- glmer(propensity_formula,
                  data = data,
                  family = binomial(link = 'logit'))
  design.matrix <- model.matrix(mod.ps)
  
  data$X.N  <- sapply(data$id, FUN = x.alters, 
                      design.matrix=design.matrix)
  data$Z.N  <- sapply(data$id, FUN = z.alters, Z=Z)
  
  Xm <- model.matrix(object = outcome_formula, data = data)
  Xe <- model.matrix(object = propensity_formula, data = data)
  
  counterfactual <- data %>% uncount(degree+1) %>% 
    dplyr::select(-k_treated) %>%
    mutate(prop_treated = 0) %>% group_by(id) %>% 
    mutate_at(vars(prop_treated), ~. + 0:(first(degree)))  %>%
    rename(k_treated = prop_treated) %>% 
    mutate(prop_treated = ifelse(degree==0, 0, k_treated/degree)) %>%
    dplyr::select(id, prop_treated, k_treated, degree) %>% 
    rename(counterfactual_prop_treated = prop_treated)
  
  
  init_drwls <- initial_values_drwls(data=data, 
                                     propensity_formula = propensity_formula, 
                                     outcome_formula=outcome_formula, 
                                     counterfactual=counterfactual, 
                                     alpha = alpha)
  
  
  
  # Call the M-estimation routine with the above initial values
  # This might take some time
  
  drwls.estimates <- estimate_drwls(data = data, 
                                  propensity_formula = propensity_formula,
                                  outcome_formula = outcome_formula,
                                  alpha = alpha,
                                  init = init_drwls,
                                  counterfactual_prop_treated = counterfactual,
                                  m=m)
  
  # Retrieve the quantities of interest and calculate causal contrasts
  
  nparm <- ncol(Xe) + 6*ncol(Xm) 
  pos.0.alpha1 <- nparm+1
  pos.1.alpha1 <- nparm+2
  pos.alpha1 <- nparm+3
  pos.0.alpha2 <- nparm+4
  pos.1.alpha2 <- nparm+5
  pos.alpha2 <- nparm+6
  pos.0.alpha3 <- nparm+7
  pos.1.alpha3 <- nparm+8
  pos.alpha3 <- nparm + 9
  
  # Table with estimates of average potential outcomes
  
  temp <- avg.po(pos = pos.0.alpha1, drwls.estimates)
  row1 <- data.frame(estimand = "Average potential outcome z=0", alpha = alpha[1],
                     estimate = temp[1], variance = temp[2], 
                     CI_lb = temp[1] - qnorm(0.975)*sqrt(temp[2]),
                     CI_ub = temp[1] + qnorm(0.975)*sqrt(temp[2]))
  temp <- avg.po(pos = pos.1.alpha1, drwls.estimates)
  row2 <- data.frame(estimand = "Average potential outcome z=1", alpha = alpha[1],
                     estimate = temp[1], variance = temp[2], 
                     CI_lb = temp[1] - qnorm(0.975)*sqrt(temp[2]),
                     CI_ub = temp[1] + qnorm(0.975)*sqrt(temp[2]))
  temp <- avg.po(pos = pos.alpha1, drwls.estimates)
  row3 <- data.frame(estimand = "Marginal average potential outcome", alpha = alpha[1],
                     estimate = temp[1], variance = temp[2], 
                     CI_lb = temp[1] - qnorm(0.975)*sqrt(temp[2]),
                     CI_ub = temp[1] + qnorm(0.975)*sqrt(temp[2]))
  
  temp <- avg.po(pos = pos.0.alpha2, drwls.estimates)
  row4 <- data.frame(estimand = "Average potential outcome z=0", alpha = alpha[2],
                     estimate = temp[1], variance = temp[2], 
                     CI_lb = temp[1] - qnorm(0.975)*sqrt(temp[2]),
                     CI_ub = temp[1] + qnorm(0.975)*sqrt(temp[2]))
  temp <- avg.po(pos = pos.1.alpha2, drwls.estimates)
  row5<- data.frame(estimand = "Average potential outcome z=1", alpha = alpha[2],
                    estimate = temp[1], variance = temp[2], 
                    CI_lb = temp[1] - qnorm(0.975)*sqrt(temp[2]),
                    CI_ub = temp[1] + qnorm(0.975)*sqrt(temp[2]))
  temp <- avg.po(pos = pos.alpha2, drwls.estimates)
  row6 <- data.frame(estimand = "Marginal average potential outcome", alpha = alpha[2],
                     estimate = temp[1], variance = temp[2], 
                     CI_lb = temp[1] - qnorm(0.975)*sqrt(temp[2]),
                     CI_ub = temp[1] + qnorm(0.975)*sqrt(temp[2]))
  
  temp <- avg.po(pos = pos.0.alpha3, drwls.estimates)
  row7 <- data.frame(estimand = "Average potential outcome z=0", alpha = alpha[3],
                     estimate = temp[1], variance = temp[2], 
                     CI_lb = temp[1] - qnorm(0.975)*sqrt(temp[2]),
                     CI_ub = temp[1] + qnorm(0.975)*sqrt(temp[2]))
  temp <- avg.po(pos = pos.1.alpha3, drwls.estimates)
  row8 <- data.frame(estimand = "Average potential outcome z=1", alpha = alpha[3],
                     estimate = temp[1], variance = temp[2], 
                     CI_lb = temp[1] - qnorm(0.975)*sqrt(temp[2]),
                     CI_ub = temp[1] + qnorm(0.975)*sqrt(temp[2]))
  temp <- avg.po(pos = pos.alpha3, drwls.estimates)
  row9 <- data.frame(estimand = "Marginal average potential outcome", alpha = alpha[3],
                     estimate = temp[1], variance = temp[2], 
                     CI_lb = temp[1] - qnorm(0.975)*sqrt(temp[2]),
                     CI_ub = temp[1] + qnorm(0.975)*sqrt(temp[2]))
  
  table1 <- rbind(row1, row2, row3, row4, row5, row6, row7, row8, row9)
  
  # Table with estimates of relevant contrasts
  # Direct effect
  
  temp <- contrast(pos.1.alpha1, pos.0.alpha1, drwls.estimates)
  row1 <- data.frame(estimand = "Direct effect", alpha0 = alpha[1], alpha1 = alpha[1],
                     estimate = temp[1], variance = temp[2], 
                     CI_lb = temp[1] - qnorm(0.975)*sqrt(temp[2]),
                     CI_ub = temp[1] + qnorm(0.975)*sqrt(temp[2]))
  temp <- contrast(pos.1.alpha2, pos.0.alpha2, drwls.estimates)
  row2 <- data.frame(estimand = "Direct effect", alpha0 = alpha[2], alpha1 = alpha[2],
                     estimate = temp[1], variance = temp[2], 
                     CI_lb = temp[1] - qnorm(0.975)*sqrt(temp[2]),
                     CI_ub = temp[1] + qnorm(0.975)*sqrt(temp[2]))
  temp <- contrast(pos.1.alpha3, pos.0.alpha3, drwls.estimates)
  row3 <- data.frame(estimand = "Direct effect", alpha0 = alpha[3], alpha1 = alpha[3],
                     estimate = temp[1], variance = temp[2], 
                     CI_lb = temp[1] - qnorm(0.975)*sqrt(temp[2]),
                     CI_ub = temp[1] + qnorm(0.975)*sqrt(temp[2]))
  
  temp <- contrast(pos.0.alpha2, pos.0.alpha1, drwls.estimates)
  row4 <- data.frame(estimand = "Indirect effect", alpha0 = alpha[1], alpha1 = alpha[2],
                     estimate = temp[1], variance = temp[2], 
                     CI_lb = temp[1] - qnorm(0.975)*sqrt(temp[2]),
                     CI_ub = temp[1] + qnorm(0.975)*sqrt(temp[2]))
  temp <- contrast(pos.0.alpha3, pos.0.alpha1, drwls.estimates)
  row5 <- data.frame(estimand = "Indirect effect", alpha0 = alpha[1], alpha1 = alpha[3],
                     estimate = temp[1], variance = temp[2], 
                     CI_lb = temp[1] - qnorm(0.975)*sqrt(temp[2]),
                     CI_ub = temp[1] + qnorm(0.975)*sqrt(temp[2]))
  temp <- contrast(pos.0.alpha3, pos.0.alpha2, drwls.estimates)
  row6 <- data.frame(estimand = "Indirect effect", alpha0 = alpha[2], alpha1 = alpha[3],
                     estimate = temp[1], variance = temp[2], 
                     CI_lb = temp[1] - qnorm(0.975)*sqrt(temp[2]),
                     CI_ub = temp[1] + qnorm(0.975)*sqrt(temp[2]))
  
  temp <- contrast(pos.1.alpha2, pos.0.alpha1, drwls.estimates)
  row7 <- data.frame(estimand = "Total effect", alpha0 = alpha[1], alpha1 = alpha[2],
                     estimate = temp[1], variance = temp[2], 
                     CI_lb = temp[1] - qnorm(0.975)*sqrt(temp[2]),
                     CI_ub = temp[1] + qnorm(0.975)*sqrt(temp[2]))
  temp <- contrast(pos.1.alpha3, pos.0.alpha1, drwls.estimates)
  row8 <- data.frame(estimand = "Total effect", alpha0 = alpha[1], alpha1 = alpha[3],
                     estimate = temp[1], variance = temp[2], 
                     CI_lb = temp[1] - qnorm(0.975)*sqrt(temp[2]),
                     CI_ub = temp[1] + qnorm(0.975)*sqrt(temp[2]))
  temp <- contrast(pos.1.alpha3, pos.0.alpha2, drwls.estimates)
  row9 <- data.frame(estimand = "Total effect", alpha0 = alpha[2], alpha1 = alpha[3],
                     estimate = temp[1], variance = temp[2], 
                     CI_lb = temp[1] - qnorm(0.975)*sqrt(temp[2]),
                     CI_ub = temp[1] + qnorm(0.975)*sqrt(temp[2]))
  
  temp <- contrast(pos.alpha2, pos.alpha1, drwls.estimates)
  row10 <- data.frame(estimand = "Overall effect", alpha0 = alpha[1], alpha1 = alpha[2],
                      estimate = temp[1], variance = temp[2], 
                      CI_lb = temp[1] - qnorm(0.975)*sqrt(temp[2]),
                      CI_ub = temp[1] + qnorm(0.975)*sqrt(temp[2]))
  temp <- contrast(pos.alpha3, pos.alpha1, drwls.estimates)
  row11 <- data.frame(estimand = "Overall effect", alpha0 = alpha[1], alpha1 = alpha[3],
                      estimate = temp[1], variance = temp[2], 
                      CI_lb = temp[1] - qnorm(0.975)*sqrt(temp[2]),
                      CI_ub = temp[1] + qnorm(0.975)*sqrt(temp[2]))
  temp <- contrast(pos.alpha3, pos.alpha2, drwls.estimates)
  row12 <- data.frame(estimand = "Overall effect", alpha0 = alpha[2], alpha1 = alpha[3],
                      estimate = temp[1], variance = temp[2], 
                      CI_lb = temp[1] - qnorm(0.975)*sqrt(temp[2]),
                      CI_ub = temp[1] + qnorm(0.975)*sqrt(temp[2]))
  
  table2 <- rbind(row1, row2, row3, row4, row5, row6, row7, row8, row9,
                  row10, row11, row12)
  
  return(list(table1, table2))
  
}


