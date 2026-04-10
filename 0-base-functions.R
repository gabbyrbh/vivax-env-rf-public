###################################
# Base functions for analysis
###################################

# [TEMPLATE FOR NEW FUNCTIONS]
# Documentation: 
# Usage: 
# Description: 
# Args/Options: 
# 
# Returns: 
# Output: ...
# [TEMPLATE FOR NEW FUNCTIONS]

###############################################
# Make asset-based PCA wealth index
###############################################
# Usage: 
# Description: Takes in a list of assets and:
#   1) replaces blank characters with NA & creates a factor level for missing data
#   2) drops rows (observations/subjects) with no asset data
#   3) drops assets with great missingness or almost no variance 
#   4) computes a PCA and takes loadings from 1st component as wealth index
#   5) computes quartiles of wealth index (HHwealth_quart)
#
# Args/Options: 
#   dfull: dataframe with all assets to be included
#   varlist: character vector of variable names from dfull to indicate which assets to include
#   reorder: T/F to indicate if wealth quartiles should be re-ordered 
# 
# Returns: A data.frame with dataid, union, and HHwealth_quart
# Output: none


assetPCA<-function(dfull, varlist, reorder=F){
  
  varlist<-c("dataid", "union", varlist)
  
  #Subset to only needed variables for subgroup analysis
  ret <- dfull %>%
    subset(select = c(varlist)) 
  
  
  #Select assets
  ret<-as.data.frame(ret) 
  id<-subset(ret, select=c("dataid", "union")) #drop subjectid
  ret<-ret[,which(!(colnames(ret) %in% c("dataid", "union")))]
  
  #Replace character blank with NA
  for(i in 1:ncol(ret)){
    ret[,i]<-ifelse(ret[,i]=="",NA,ret[,i])
  } 
  
  #drop rows with no asset data
  id<-id[rowSums(is.na(ret[,3:ncol(ret)])) != ncol(ret)-2,]  
  ret<-ret[rowSums(is.na(ret[,3:ncol(ret)])) != ncol(ret)-2,]  
  
  
  #Drop assets with great missingness
  for(i in 1:ncol(ret)){
    cat(colnames(ret)[i],"\n")
    print(table(is.na(ret[,i])))
    print(class((ret[,i])))
  }
  
  #create level for for missing factor levels
  table(is.na(ret))
  for(i in 1:ncol(ret)){
    ret[,i]<-as.character(ret[,i])
    ret[is.na(ret[,i]),i]<-"miss"
    ret[,i]<-as.factor(ret[,i])
    
  }
  table(is.na(ret))
  
  #Convert factors into indicators
  ret<-droplevels(ret)
  ret<-design_matrix(ret)
  
  
  #Remove columns with almost no variance
  if(length(nearZeroVar(ret))>0){
    ret<-ret[,-nearZeroVar(ret)]
  }
  
  ## Convert the data into matrix ##
  ret<-as.matrix(ret)
  
  ##Computing the principal component using eigenvalue decomposition ##
  princ.return <- princomp(ret) 
  
  ## To get the first principal component in a variable ##
  load <- loadings(princ.return)[,1]   
  
  pr.cp <- ret %*% load  ## Matrix multiplication of the input data with the loading for the 1st PC gives us the 1st PC in matrix form. 
  
  HHwealth <- as.numeric(pr.cp) ## Gives us the 1st PC in numeric form in pr.
  
  #Create 4-level household wealth index
  quartiles<-quantile(HHwealth, probs=seq(0, 1, 0.25))
  print(quartiles)
  ret<-as.data.frame(ret)
  ret$HHwealth_quart<-rep(1, nrow(ret))
  ret$HHwealth_quart[HHwealth>=quartiles[2]]<-2
  ret$HHwealth_quart[HHwealth>=quartiles[3]]<-3
  ret$HHwealth_quart[HHwealth>=quartiles[4]]<-4
  table(ret$HHwealth_quart)
  ret$HHwealth_quart<-factor(ret$HHwealth_quart)
  
  if(reorder==T){
    levels(ret$HHwealth_quart)<-c("Wealth Q4","Wealth Q3","Wealth Q2","Wealth Q1")
    ret$HHwealth_quart<-factor(ret$HHwealth_quart, levels=c("Wealth Q1", "Wealth Q2","Wealth Q3","Wealth Q4"))
  }else{
    levels(ret$HHwealth_quart)<-c("Wealth Q1", "Wealth Q2","Wealth Q3","Wealth Q4")
  }
  
  #Table assets by pca quartile to identify wealth/poverty levels
  d<-data.frame(id, ret)
  wealth.tab <- d %>% subset(., select=-c(dataid)) %>%
    group_by(HHwealth_quart) %>%
    summarise_all(list(mean = mean)) %>% as.data.frame()
  print(wealth.tab)
  
  #Save just the wealth data
  pca.wealth<-d %>% subset(select=c(dataid, union, HHwealth_quart))
  
  # pca.wealth$dataid<-as.character(pca.wealth$dataid)
  
  d <-dfull %>% subset(., select=c("dataid","union"))
  # d$dataid<-as.numeric(d$dataid)
  d<-left_join(d, pca.wealth, by=c("dataid","union"))
  return(d)
}


###############################################
# sparsity check
###############################################
# Check whether the number of cases per variable is sufficient to perform analysis
# y:                    outcome variable name, as numeric or factor
# a:                    risk factor variable name, as numeric or factor
# w:                    covariate variable name(s), as a data frame
# y_per_variable:       Default 10. Number of cases per variable required to perform analysis 

# returns T if number of cases per variable is sufficient to perform analysis 
# returns F if not 

check_sparsity = function(y, a, w, y_per_variable = 10){
  
  assert_that(length(y) == length(a))
  assert_that(length(y) == nrow(w))
  assert_that(length(a) == nrow(w))
  
  df = data.frame(y = y, a = a)
  
  df = bind_cols(df, w)
  
  # subset to complete cases
  df_complete = df[complete.cases(df),]
  
  # subset to categorical factor variables 
  df_numeric = df_complete %>% dplyr::select_if(is.factor)
  
  # count number of observations within each stratum
  dfN = df_numeric %>% group_by_all() %>% count()
  
  # check whether cases per variable meets required number
  check = ifelse(min(dfN$n) >= y_per_variable, T, F) 
  
  return(check)
}

###################################
# check correlation function 
###################################
# Documentation: check_cor
# Usage:         check_cor(y, wlist, data, threshold)
# Description:   check whether covariates are correlated with the outcome
#                using a spearman's correlation coefficient. 
# Args/Options:  
# y:              outcome variable name, as character
# wlist:          covariate variable name(s), as character string
# data:           data frame including outcome and covariates
# threshold:      required value of spearman correlation coefficient
#                 in order for the covariate to be considered associated
#                 with the outcome 
# Returns:        list of covariates associated with the outcome, as character string
# Output:         none

check_cor = function(y, wlist, data, threshold){
  
  cor_values = map_dbl(wlist, function(x) cor(data[[y]], data[[x]], method = "spearman"))
  cor_df = data.frame(
    w = wlist,
    cor = cor_values
  ) %>%
    filter(cor>threshold) %>%
    dplyr::select(w) %>% pull() %>% as.character()
  
  return(cor_df)
  
}


###############################################
# check for residual spatial autocorrelation using Moran's coefficients
###############################################
# Usage:        check_autocorr(gam_fit$gam, lldata = d[,c("qgpslong", "qgpslat")], nbc = 8)
# Description:  checks for residual spatial autocorrelation using Moran's coefficients
# Args/Options: 
#   fit:    The gam model fit, from which residuals will be extracted
#   lldata: A 2-column data.frame of spatial coordinates
#   nbc:    Default 10, number of bins passed to nbclass of the pgirmess::correlog function

# Returns: A list with two elements 1) a logical value if there is residual spatial autocorrelation and 
#           2) the correlogram
# Output: A statement whether residual spatial autocorrelation was detected and whether non-spatial model 
#           is saved or proceeding with spatial model.
check_autocorr <- function(fit, lldata, nbc=10){
  # Compute correlogram of the residuals
  lldata = lldata %>% 
    mutate(res = residuals(fit)) %>% 
    filter(!is.na(res)) %>%
    select(-res)
  
  cor_r <- correlog(coords=lldata,
                    z=residuals(fit)[!is.na(residuals(fit))],
                    method="Moran", nbclass=nbc) #na.action = na.omit
  
  correlograms <- as.data.frame(cor_r)

  # define boolean for whether any spatial autocorrelation is present 
  spatial_autocorr <- any(correlograms$p.value <= 0.05)
  
  if(spatial_autocorr) print("Residual spatial autocorrelation is present. Proceeding with spatial model.")
  if(!spatial_autocorr) print("No residual spatial autocorrelation. Saving non-spatial model.")
  
  return(list(spatial_autocorr = spatial_autocorr,
              correlograms = correlograms))
}


###############################################
# get_covariate_list function
###############################################
# Documentation:    get_covariate_list
# Usage:            get_covariate_list(outcome_cat = "sth", risk_factor = "vpd")
# Description:      Will determine the list of covariates and confounders to include in adjusted analyses 
#                     for the supplied risk factor.  These were defined in the analysis plan.
# Args/Options: 
#   outcome:      The pathogen outcome being evaluated 
#   risk_factor:  The meteorological or environmental risk factor being evaluated
# Returns: A list of adjustment covariates and confounders.  These will be pre-screened with a Likelihood ratio test
# Output: A statement whether residual spatial autocorrelation was detected and whether non-spatial model 
#           is saved or proceeding with spatial model.

get_covariate_list = function(outcome_cat = c("giardia", "pathogens", "sth", "diarrhea"), risk_factor) {

  cov_list = c("sex", "aged_C", "HHwealth_quart") #"any_holiday", "season" 
  
  ## Add in tr covariate  
  #  Diarrhea & giardia datasets include samples collected at baseline before interventions were delivered.  So all baseline samples will be coded as Control (as their household had not received an intervention at the time of sample collection)
  #  sth dataset only includes children that have received interventions already, so it's the assigned intervention arm
  #  pathogen dataset has bias in season and control arm data collection, so we run into sparsity issues when adjusting.  Will make binary WASH and nutrition variables to replace tr.
  if (outcome_cat %in% c("giardia", "diarrhea")) {
    cov_list = c(cov_list, "tr_received")
  } else if (outcome_cat == "sth") {
    cov_list = c(cov_list, "tr")
  } else if (outcome_cat == "pathogens") {
    cov_list = c(cov_list, "wsh", "nutrition")
  }
  
  if (outcome_cat %in% c("pathogens")) { 
    cov_list = c(cov_list, "abx_7d")
  } else if (outcome_cat == "sth") {
    cov_list = c(cov_list, "dw_2wks")
  }
  
  if (risk_factor == "vpd_C") {
    cov_list = c(cov_list, "monthly_precipitation_C")
  } else if (risk_factor == "flow_accum_median") {
    cov_list = c(cov_list, "ncow_c_C", "ngoat_c_C", "nchicken_c_C", "monthly_precipitation_C")
  } else if (risk_factor %in% c("distance_from_any_surface_water_C", "distance_from_seasonal_surface_water_C", "distance_from_ephemeral_surface_water_C")) {
    cov_list = c(cov_list, "flow_accum_median", "ncow_c_C", "ngoat_c_C", "nchicken_c_C") # "flow_accumulation_C"
  } else if (grepl("prop_detected", risk_factor) | grepl("max_months_water", risk_factor)) {
    cov_list = c(cov_list, "monthly_precipitation_C")
  } else if (risk_factor == "land_use_type") {
    cov_list = c(cov_list, "momeduy_C", "dadeduy_C")
  } else if (risk_factor == "evi_C") {
    cov_list = c(cov_list, "momeduy_C", "dadeduy_C", "ncow_c_C", "ngoat_c_C", "nchicken_c_C", "vpd_C", "flow_accum_median") # "flow_accumulation_C"
  } else if (risk_factor == "pop_density_C") {
    cov_list = c(cov_list, "momeduy_C", "dadeduy_C", "Nhh_median")
  }
  
  return(cov_list)
}


###############################################
# check_collinearity function
###############################################
# Documentation:    check_collinearity
# Usage:            check_collinearity(risk_factor = "avgtemp_7_days", covariate = "aged", df = d_diarrhea)
# Description:    Depending on whether the risk_factor and covariate are continuous or categorical, will
#                 evaluate collinearity with a Pearson correlation, ANOVA, or chi-squared test. 
# Args/Options: 
#   risk_factor:  The meterological or environmental risk factor being evaluated
#   covariate:    The covariate being screened for collinearity
#   df:           The data.frame in which collinearity will be evaluated
# Returns: A list of correlation result and a plot of the association between risk_factor and covariate.
# Output: 

check_collinearity <- function(risk_factor, covariate, df) {

  if (is.numeric(df[,risk_factor]) & is.numeric(df[,covariate])) {
    res = cor.test(df[,risk_factor], df[,covariate], method = "pearson", alternative = "two.sided")
    plot = ggplot(df, aes_string(x = covariate, y = risk_factor)) + 
      geom_point(shape = 1, alpha = 0.5)
    out = list(res = res, 
               plot = plot) 
  } else if (is.numeric(df[,risk_factor]) & !is.numeric(df[,covariate])) {
    fm <- as.formula(paste0(risk_factor, "~", covariate))
    res = aov(fm, df)
    plot = ggplot(df, aes_string(x = covariate, y = risk_factor)) + 
      geom_boxplot(outlier.shape = NA) + 
      geom_jitter(shape = 1, width = 0.25, alpha = 0.5)
    out = list(res = summary(res), 
               plot = plot)
  } else if (is.numeric(df[,covariate]) & !is.numeric(df[,risk_factor])) {
    fm <- as.formula(paste0(covariate, "~", risk_factor))
    res = aov(fm, df)    
    plot = ggplot(df, aes_string(x = risk_factor, y = covariate)) + 
      geom_boxplot(outlier.shape = NA) + 
      geom_jitter(shape = 1, width = 0.25, alpha = 0.5)
    out = list(res = summary(res), 
               plot = plot)
  } else if (!is.numeric(df[,risk_factor]) & !is.numeric(df[,covariate])) {
    res = chisq.test(df[,risk_factor], df[,covariate])
    plot = NA
    out = list(res = res, 
               plot = plot)
  }
  return(out)
}


###############################################
# washb_prescreen_mod function
###############################################
# Documentation:    washb_prescreen_mod
# Usage:            washb_prescreen_mod()
# Description:    This function modifies the wasb_prescreen function to return the list of covariates in
#                 order of increasing p-value from the likelihood ratio test. This is useful when sparsity 
#                 is a problem, so that covariates can be systematically dropped based on those with
#                 the least strong association with the outcome.  
# Args/Options: 
#   Y:      Outcome variable (continuous, such as LAZ, or binary, such as diarrhea)
#   Ws:     data frame that includes candidate adjustment covariates to screen
#   family: GLM model family (gaussian, binomial, poisson, or negative binomial). Use "neg.binom" 
#             for Negative binomial.
#   pval:   The p-value threshold: any variables with a p-value from the lielihood ratio test 
#             below this threshold will be returned. Defaults to 0.2
#   print:  Logical for whether to print function output, defaults to TRUE.
# Returns: A list which includes covariates that passed the LR test ranked in order of increasing p-value.
# Output: 
#   1) Result of likelihood ratio test for pre-screening covariates
#   2) List of covariates passing LR screening, in ranked order of increasing p-value

washb_prescreen_mod = function (Y, Ws, family = c("gaussian", "binomial", "neg.binom", "poisson"), pval = 0.2, print = TRUE) {

  require(lmtest)
  if (family[[1]] == "neg.binom") {
    require(MASS)
  }
  if (pval > 1 | pval < 0) {
    stop("P-value threshold not set between 0 and 1.")
  }
  Ws <- as.data.frame(Ws)
  dat <- data.frame(Ws, Y)
  dat <- dat[complete.cases(dat), ]
  nW <- ncol(Ws)
  LRp <- matrix(rep(NA, nW), nrow = nW, ncol = 1)
  rownames(LRp) <- names(Ws)
  colnames(LRp) <- "P-value"
  if (family[[1]] != "neg.binom") {
    for (i in 1:nW) {
      # i = 3
      dat$W <- dat[, i]
      if (class(dat$W) == "factor" & dim(table(dat$W)) ==
          1) {
        fit1 <- fit0 <- glm(Y ~ 1, data = dat, family = family)
      }
      else {
        fit1 <- glm(Y ~ W, data = dat, family = family)
        fit0 <- glm(Y ~ 1, data = dat, family = family)
      }
      LRp[i] <- lrtest(fit1, fit0)[2, 5]
    }
  }
  else {
    if (!requireNamespace("MASS", quietly = TRUE)) {
      stop("Pkg needed for this function to work. Please install it.",
           call. = FALSE)
    }
    else {
      for (i in 1:nW) {
        dat$W <- dat[, i]
        if (class(dat$W) == "factor" & dim(table(dat$W)) ==
            1) {
          fit1 <- fit0 <- glm(Y ~ 1, data = dat, family = family)
        }
        else {
          fit1 <- glm.nb(Y ~ W, data = dat, family = family)
          fit0 <- glm.nb(Y ~ 1, data = dat, family = family)
        }
        LRp[i] <- lrtest(fit1, fit0)[2, 5]
      }
    }
  }
  p20 <- ifelse(LRp < pval, 1, 0)
  if (print == TRUE) {
    cat("\nLikelihood Ratio Test P-values:\n")
    print(round(LRp, 5))
    if (sum(p20) > 0) {
      LRps <- matrix(LRp[p20 == 1, ], ncol = 1)
      rownames(LRps) <- names(Ws)[p20 == 1]
      colnames(LRps) <- "P-value"
      cat(paste("\n\nCovariates selected (P<", pval, "):\n",
                sep = ""))
      print(LRps)
    }
    else {
      cat(paste("\nNo covariates were associated with the outcome at P<",
                pval))
    }
  }
  out = data.frame(LRp) %>%
    rownames_to_column("covariate") %>%
    filter(P.value < pval) %>%
    arrange(P.value)
  # return(names(Ws)[p20 == 1])
  return(out$covariate)
}

###############################################
# fit_gam function
###############################################
# Documentation:    fit_gam
# Usage:            fit_gam(df = d_sth, 
#                           y = as.character(sth_tbl_row[["outcome"]]), 
#                           a = as.character(sth_tbl_row[["risk_factor"]]), 
#                           w = c("sex", "aged"), 
#                           family = "binomial")
# Description:      Function that 1) checks there are sufficient observations for analysis, 
#                                 2) screens covariates for association with the outcome, 
#                                 3) subsets the data to include complete observations,
#                                 4) fits a GAM ignoring spatial variation,
#                                 5) verifies that the model converged and throws an error if not,
#                                 6) checks for residual spatial autocorrelation,
#                                 7) fits a spatial GAM if residual spatial autocorrelation is present,
#                                 8) verifies the spatial GAM converged and throws and error if not.
#
# Args/Options: 
#   df:               Data.frame containing all variables needed for analysis
#   y:                Character argument that specifies the outcome variable
#   a:                Character argument that specifies the exposure variable
#   w:                Default NULL. A list or single character vector that specifies the names of covariates or other 
#                       variables to adjust for in the gam function. This list will be screened with a likelihood ratio 
#                       test for association with the outcome before being included in the model.
#   y_per_variable:   Defaut 10. Number of cases per variable required to perform analysis, passed to check_sparsity function.
#   random_intercept: Character argument that specifies the variable to be used as a random intercept. Default: "dataid" 
#                       (for all main cohort analyses - diarrhea, sth, giardia) and "clusterid" for EE cohort diaarrhea and pathogen analysis.
#   family:           The family object specifying the distribution and link to use in fitting the gam.  Acceptable responses  
#                       are currently "binomial" and "gaussian", although in practice any allowed by family.mgcv could be possible.
# 
# Returns: A list with the following elements:
#   gam_fit:                          gam model fit object 
#   spatial_gam_fit:                  spatial gam model fit if spatial autocorrelation was detected, else NA
#   spatial_autocorr:                 True/False if spatial autocorrelation was detected in the gam_fit
#   spatial_autocorr_coeff:           the correlogram of the Moran's coefficients from the spatial autocorrelation check on the gam_fit model
#   residual_spatial_autocorr:		    True/False if spatial autocorrelation was detected still in the spatial_gam_fit 
#   residual_spatial_autocorr_coeff:  the correlogram of the Moran's coefficients from the spatial autocorrelation check on the spatial_gam_fit model
#   y:                                the outcome variable specified
#   a:                                the exposure variable specified
#   w_included:                       the list of covariates included in the model (i.e. those that passed screening with the likelihood ratio test)
#   formula:                          the model formula used in the gam.
#   family:                           the family specified by the user indicating the distribution used to fit the gam. 
#   model_input_data:                 the data.frame used for model fitting (with observations dropped that had missing variables)
#   drop_outcome:				              the number of observations dropped due to missing outcome
#   drop_predictors:				          the number of observations dropped due to missing predictors
#   message:                          One of 5 messages which indicate the outcome of the attempted gam model fit. Possible messages include:
#                                       -"Adjusted model not run because no covariates were associated with the outcome at P< 0.1" [no model was fit]
#                                       -"Adjusted/interaction model could not be run due to data sparsity" [no model was fit]
#                                       -"Non-spatial GAM failed to converge" [no model was fit]
#                                       -"Non-spatial gam fit successfully, no spatial autocorrelation detected so no spatial gam fit" [a non-spatial model was fit]
#                                       -"Spatial gam fit successfully" [both a non-spatial and a spatial model were fit]
#
# Output: 
#   1) Result of likelihood ratio test for pre-screening covariates 
#   2) Error message if sparsity check fails
#   3) Error message if gam fails to converge 
#   4) Message reporting if any residual spatial autocorrelation remains
#   5) If residual spatial autocorrelation remains, output of gam model fit including spatial autocorrelation
#   6) Error message if spatial gam fails to converge 
#   7) Message reporting if spatial GAM has larger AIC than non-spatial GAM
#
fit_gam = function(df, y, a, w = NULL, interaction_term = NULL, y_per_variable = 10, random_intercept = c("dataid", "clusterid"), family = c("binomial", "gaussian", "poisson")) {

  #---------------------------------------------
  # define new variable for message and set to null (this will trigger function to exit at certain spots when no longer null)
  #---------------------------------------------# 
  message = NULL
  #---------------------------------------------
  # check that the family is one of the allowed responses (no typos, no Poisson, etc)
  #---------------------------------------------
  random_intercept <- match.arg(random_intercept)
  family <- match.arg(family)
  y = sym(y)
  #---------------------------------------------
  # drop rows with missing outcome
  #---------------------------------------------
  d = df %>% filter(!is.na({{y}}))
  drop_outcome = dim(df)[1] - dim(d)[1]
  message(paste0("Observations lost due to missing outcome: ", drop_outcome))
  
  #---------------------------------------------
  # covariate screening with likelihood ratio test, keeping p<=0.1
  #---------------------------------------------
  if (!is.null(w)) {
    keep_w = washb_prescreen_mod(d %>% pull(y), Ws = d %>% select(all_of(w)), family = family, pval = 0.1, print = T)
    if (length(keep_w) == 0) message = "Adjusted model not run because no covariates were associated with the outcome at P< 0.1"
  } else {
    keep_w = NULL
  }
  
  #---------------------------------------------
  # check there are sufficient observations for analysis
  #---------------------------------------------
  sparsity_check = try(validate_that(check_sparsity(y = d %>% pull(y),
                                                  a = d[, a],
                                                  w = data.frame(d[,all_of(keep_w)]),
                                                  y_per_variable = y_per_variable),
                                   msg = "Data sparsity check failed. Choose new covariate set."))
  
  if (sparsity_check != TRUE & "season" %in% keep_w) {
    keep_w = keep_w[!keep_w %in% "season"] 
    sparsity_check = try(validate_that(check_sparsity(y = d %>% pull(y),
                                                    a = d[, a],
                                                    w = data.frame(d[,all_of(keep_w)]),
                                                    y_per_variable = y_per_variable),
                                     msg = "Data sparsity check failed. Choose new covariate set."))
  }
  while (sparsity_check != TRUE & length(keep_w) > 1) {
    keep_w = head(keep_w, -1)
    sparsity_check = try(validate_that(check_sparsity(y = d %>% pull(y),
                                                    a = d[, a],
                                                    w = data.frame(d[,all_of(keep_w)]),
                                                    y_per_variable = y_per_variable),
                                     msg = "Data sparsity check failed. Choose new covariate set."))
  }
  if (sparsity_check != TRUE & length(keep_w) == 1) {
    message = "Adjusted/interaction model could not be run due to data sparsity"
  }
  
  if (!is.null(message)) {
    return(list(gam_fit = NA,
                spatial_gam_fit = NA,
                spatial_autocorr = NA,
                spatial_autocorr_coeff = NA,
                residual_spatial_autocorr = NA,
                residual_spatial_autocorr_coeff = NA,
                y = y, 
                a = a, 
                w_included = keep_w,
                formula = NA,
                family = family,
                model_input_data = d, 
                drop_outcome = drop_outcome, 
                drop_predictors = NA,
                message = message))
  } 
  #---------------------------------------------
  # build formula
  #---------------------------------------------
  if (is.numeric(d[,a])) {
    formula_start = paste(y, paste0("s(", a, ", bs=\"cs\")"), sep = " ~ ")
  } else {
    formula_start = paste(y, a, sep = " ~ ")
  }
  keep_w_terms = NULL
  if (!is.null(keep_w)) {
    ## function to write s() for numeric variables
    keep_w_terms = map_chr(keep_w, function(x) {
      ifelse(is.numeric(d[,x]), paste0("s(", x, ", bs=\"cs\")"), x)
    })
  } 
  if (!is.null(keep_w_terms)) {
    formula = as.formula(paste(formula_start,
                               paste(keep_w_terms, collapse = " + "),
                               # "s(clusterid, bs = \"re\", by = dummy)",
                               sep = " + "))
    
    formula_spatial =  as.formula(paste(formula_start,
                                        paste(keep_w_terms, collapse = " + "),    
                                        # thin plate smoother for lat & long
                                        "s(qgpslong, qgpslat, bs = \"tp\")",
                                        # "s(clusterid, bs = \"re\", by = dummy)",
                                        sep = " + "))
  } else {
    formula = as.formula(formula_start)
    
    formula_spatial = as.formula(paste(formula_start,
                                       # thin plate smoother for lat & long
                                       "s(qgpslong, qgpslat, bs = \"tp\")", 
                                       # "s(clusterid, bs = \"re\", by = dummy)",
                                       sep = " + "))
  }
  
  
  #---------------------------------------------
  # drop rows with missing risk factor or covariates
  #---------------------------------------------
  d = d %>% 
    select(all_of(y), all_of(a), all_of(keep_w), qgpslat, qgpslong, all_of(random_intercept)) %>%
    na.omit() 
  
  drop_predictors = dim(df)[1] - drop_outcome - dim(d)[1]
  message(paste0("Observations lost due to missing predictor or covariates: ", drop_predictors))
  
  random = as.formula(paste0("~(1|", random_intercept, ")"))
  
  
  #---------------------------------------------
  # fit GAM ignoring spatial variation
  #---------------------------------------------
  gam_fit = try(gamm4(formula, 
                      family = family,
                      data = d,
                      REML = T,
                      na.action = na.exclude, 
                      random = random
  ))

  #---------------------------------------------
  # confirm the model converged
  #---------------------------------------------
  print(paste0("Non-spatial gam fit: ", y, " ~ ", a))
  if(length(gam_fit) > 1 & !class(gam_fit) %in% c("try-error", "try-warning", "try-message")){
    print(summary(gam_fit$gam))
  } else {
    print(gam_fit)
    message = "Non-spatial GAM failed to converge"
    return(list(gam_fit = NA,
                spatial_gam_fit = NA,
                spatial_autocorr = NA,
                spatial_autocorr_coeff = NA,
                residual_spatial_autocorr = NA,
                residual_spatial_autocorr_coeff = NA,
                y = y, 
                a = a, 
                w_included = keep_w,
                formula = formula,
                family = family,
                model_input_data = d, 
                drop_outcome = drop_outcome, 
                drop_predictors = drop_predictors,
                message = message))
  }
  #---------------------------------------------
  # check for residual spatial autocorrelation
  #---------------------------------------------
  check_spatial_autocorr <- check_autocorr(gam_fit$gam, lldata = d[,c("qgpslong", "qgpslat")], nbc = 8)
  
  #---------------------------------------------
  # if residual spatial autocorrelation is present,
  # proceed with spatial GAM
  #---------------------------------------------
  if(!check_spatial_autocorr$spatial_autocorr){
    return(list(gam_fit = gam_fit,
                spatial_gam_fit = NA,
                spatial_autocorr = check_spatial_autocorr$spatial_autocorr,
                spatial_autocorr_coeff = check_spatial_autocorr$correlograms,
                residual_spatial_autocorr = NA,
                residual_spatial_autocorr_coeff = NA,
                y = y, 
                a = a, 
                w_included = keep_w,
                formula = formula,
                family = family,
                model_input_data = d, 
                drop_outcome = drop_outcome, 
                drop_predictors = drop_predictors,
                message = "Non-spatial gam fit successfully, no spatial autocorrelation detected so no spatial gam fit"))
  }else{
    
    #---------------------------------------------
    # If g is a factor then s(g,bs="re") produces a
    # random coefficient for each level of g, with
    # the random coefficients all modelled as i.i.d. normal.
    #---------------------------------------------
    
    
    spatial_gam_fit = try(gamm4(formula_spatial, 
                                family = family,
                                data = d,
                                REML = T,
                                na.action = na.exclude, 
                                random = random
    ))
    
    #---------------------------------------------
    # confirm the model converged
    #---------------------------------------------
    print("Spatial gam fit")
    summary(spatial_gam_fit)
    if(length(spatial_gam_fit) > 1 & !class(spatial_gam_fit) %in% c("try-error", "try-warning", "try-message")){
      print(summary(spatial_gam_fit$gam))
    } else {
      print(spatial_gam_fit)
      stop("Spatial GAM failed to converge")
    }
    
    # compare AIC for model with and without spatial smoother
    if (AIC(gam_fit$mer) < AIC(spatial_gam_fit$mer)) {
      print("Even though residual spatial autocorrelation was present, spatial GAM has larger AIC than non-spatial GAM. Check model.")
    }
    
    #---------------------------------------------
    # check for residual spatial autocorrelation 
    #   **after** adjusting for lat and long
    #---------------------------------------------
    check_residual_spatial_autocorr <- check_autocorr(spatial_gam_fit$gam, lldata = d[,c("qgpslong", "qgpslat")], nbc = 8)
    
    
    return(list(gam_fit = gam_fit,
                spatial_gam_fit = spatial_gam_fit,
                spatial_autocorr = check_spatial_autocorr$spatial_autocorr,
                spatial_autocorr_coeff = check_spatial_autocorr$correlograms,
                residual_spatial_autocorr = check_residual_spatial_autocorr$spatial_autocorr,
                residual_spatial_autocorr_coeff = check_residual_spatial_autocorr$correlograms,
                y = y, 
                a = a, 
                w_included = c(keep_w, "qgpslong", "qgpslat"),
                formula = formula_spatial,
                family = family,
                model_input_data = d, 
                drop_outcome = drop_outcome, 
                drop_predictors = drop_predictors,
                message = "Spatial gam fit successfully"))
    # end of spatial autocorrelation if statement
  }
  
}

###############################################
# fit_gam function for models with interaction terms
# Should eventually be merged with the above, but I had insufficient time
#  to test this and ensure it didn't break any of the above code so it's still separate.
#   Main differenes include that 
#   1) the interaction term is added to the sparsity check
#   2) an interaction term is added to the formula, 
#   3) the function output includes a list item for the interaction term.  
###############################################
fit_gam_interaction = function(df, y, a, w = NULL, interaction_term = NULL, y_per_variable = 10, random_intercept = c("dataid", "clusterid"), family = c("binomial", "gaussian", "poisson")) {

  #---------------------------------------------
  # define new variable for message and set to null (this will trigger function to exit at certain spots when no longer null)
  #---------------------------------------------# 
  message = NULL
  
  #---------------------------------------------
  # check that the family is one of the allowed responses (no typos, no Poisson, etc)
  #---------------------------------------------
  random_intercept <- match.arg(random_intercept)
  family <- match.arg(family)
  y = sym(y)
  #---------------------------------------------
  # drop rows with missing outcome
  #---------------------------------------------
  d = df %>% filter(!is.na({{y}}))
  drop_outcome = dim(df)[1] - dim(d)[1]
  message(paste0("Observations lost due to missing outcome: ", drop_outcome))
  
  #---------------------------------------------
  # covariate screening with likelihood ratio test, keeping p<=0.1
  #---------------------------------------------
  if (!is.null(w)) {
    keep_w = washb_prescreen_mod(d %>% pull(y), Ws = d %>% select(all_of(w)), family = family, pval = 0.1, print = T)
    if (length(keep_w) == 0) message = "Adjusted model not run because no covariates were associated with the outcome at P< 0.1"
  } else {
    keep_w = NULL
  }
  
  
  #---------------------------------------------
  # check there are sufficient observations for analysis
  #---------------------------------------------
  if (is.null(interaction_term)) {
    check_w = keep_w
  } else {
    check_w = c(keep_w, interaction_term)
  }
  
  sparsity_check = try(validate_that(check_sparsity(y = d %>% pull(y),
                                                  a = d[, a],
                                                  w = data.frame(d[,all_of(check_w)]),
                                                  y_per_variable = y_per_variable),
                                   msg = "Data sparsity check failed. Choose new covariate set."))
  
  if (sparsity_check != TRUE & "season" %in% check_w) {
    check_w = check_w[!check_w %in% "season"] 
    sparsity_check = try(validate_that(check_sparsity(y = d %>% pull(y),
                                                    a = d[, a],
                                                    w = data.frame(d[,all_of(check_w)]),
                                                    y_per_variable = y_per_variable),
                                     msg = "Data sparsity check failed. Choose new covariate set."))
  }
  while (sparsity_check != TRUE & length(check_w) > 1) {
    check_w = head(check_w, -1)
    sparsity_check = try(validate_that(check_sparsity(y = d %>% pull(y),
                                                    a = d[, a],
                                                    w = data.frame(d[,all_of(check_w)]),
                                                    y_per_variable = y_per_variable),
                                     msg = "Data sparsity check failed. Choose new covariate set."))
  }
  if (sparsity_check != TRUE & length(check_w) == 1) {
    message = "Adjusted/interaction model could not be run due to data sparsity"
    message("Adjusted/interaction model could not be run due to data sparsity")
  }
  
  if (!is.null(message)) {
    return(list(gam_fit = NA,
                spatial_gam_fit = NA,
                spatial_autocorr = NA,
                spatial_autocorr_coeff = NA,
                residual_spatial_autocorr = NA,
                residual_spatial_autocorr_coeff = NA,
                y = y, 
                a = a, 
                w_included = keep_w,
                interaction_term = interaction_term,
                formula = NA,
                family = family,
                model_input_data = d, 
                drop_outcome = drop_outcome, 
                drop_predictors = NA,
                message = message))
  } 
  
  #---------------------------------------------
  # build formula
  #---------------------------------------------
  ## Formula start with outcome and risk factor
  if (is.numeric(d[,a])) {
    # formula_start = paste(y, paste0("s(",a,", bs=\"cs\")"), sep = " ~ ")
    formula_start = paste(y, "", sep = " ~ ") #, paste0("s(",a,")")
  } else {
    # formula_start = paste(y, a, sep = " ~ ")
    formula_start = paste(y, "", sep = " ~ ") #a,
  }
  
  ## Format adjustment covariates
  keep_w_terms = NULL
  if (!is.null(keep_w)) {
    ## function to write s() for numeric variables
    keep_w_terms = map_chr(keep_w, function(x) {
      ifelse(is.numeric(d[,x]), paste0("s(", x, ", bs=\"cs\")"), x)
    })
  } 
  
  ## Format interaction, if applicable
  interaction_formula = NULL
  if (!is.null(interaction_term)) {
    if (is.numeric(d[,a]) & is.factor(d[,interaction_term])) {
      interaction_formula = paste0("s(", a, ", by = ", interaction_term, ", bs=\"cs\", m = 1)", " + ", interaction_term, " + s(", a,", bs=\"cs\")")
    } else if (is.numeric(d[,a]) & is.numeric(d[,interaction_term])) {
      interaction_formula = paste0("t2(", a, ",", interaction_term, ", bs=\"cs\")")
    } else if (is.factor(d[,a]) & is.numeric(d[,interaction_term])) {
      print("There should be no factor risk_factors with continuous interaction_terms")
    } else if (is.factor(d[,a]) & is.factor(d[,interaction_term])) {
      interaction_formula = paste0(a, "*", interaction_term)
    } 
  }
  ## Make complete formula
  if (!is.null(keep_w_terms) & is.null(interaction_formula)) {
    print("here1")
    formula = as.formula(paste(formula_start,
                               paste0(ifelse(is.numeric(d[,a]), paste0("s(",a,", bs=\"cs\") +"), paste0(a, " + ")), ifelse(length(keep_w_terms) > 1, paste(keep_w_terms, collapse = " + "), keep_w_terms)),
                               sep = ""))
    
    formula_spatial =  as.formula(paste(paste0(formula_start,
                                               paste0(ifelse(is.numeric(d[,a]), paste0("s(",a,", bs=\"cs\") +"), paste0(a, " + ")), ifelse(length(keep_w_terms) > 1, paste(keep_w_terms, collapse = " + "), keep_w_terms))),
                                        # thin plate smoother for lat & long
                                        "s(qgpslong, qgpslat, bs = \"tp\")",
                                        sep = " + "))
  } else if (!is.null(keep_w_terms) & !is.null(interaction_formula)) {
    print("here2")
    formula = as.formula(paste(paste0(formula_start, interaction_formula), 
                               paste(ifelse(length(keep_w_terms) > 1, paste(keep_w_terms, collapse = " + "), keep_w_terms)),
                               sep = " + "))
    
    formula_spatial =  as.formula(paste(paste0(formula_start, interaction_formula), 
                                        paste(ifelse(length(keep_w_terms) > 1, paste(keep_w_terms, collapse = " + "), keep_w_terms)),
                                        # thin plate smoother for lat & long
                                        "s(qgpslong, qgpslat, bs = \"tp\")",
                                        sep = " + "))
    
  } else if (is.null(interaction_formula)) {
    print("here3")
    formula = as.formula(paste(formula_start,
                               paste0(ifelse(is.numeric(d[,a]), paste0("s(",a,", bs=\"cs\")"), a))))
    
    formula_spatial = as.formula(paste(formula_start,
                                       paste0(ifelse(is.numeric(d[,a]), paste0("s(",a,", bs=\"cs\") +"), paste0(a, " + ")),
                                              # thin plate smoother for lat & long
                                              "s(qgpslong, qgpslat, bs = \"tp\")"), 
                                       # "s(clusterid, bs = \"re\", by = dummy)",
                                       sep = ""))
  } else {
    print("here4")
    formula = as.formula(paste(formula_start, 
                               interaction_formula,
                               sep = ""))
    
    formula_spatial = as.formula(paste(paste0(formula_start,
                                              interaction_formula), 
                                       # thin plate smoother for lat & long
                                       "s(qgpslong, qgpslat, bs = \"tp\")", 
                                       # "s(clusterid, bs = \"re\", by = dummy)",
                                       sep = " + "))
  }
  
  #---------------------------------------------
  # drop rows with missing risk factor or covariates
  #---------------------------------------------
  d = d %>% 
    select(all_of(y), all_of(a), all_of(check_w), qgpslat, qgpslong, all_of(random_intercept)) %>%
    na.omit() 
  
  drop_predictors = dim(df)[1] - drop_outcome - dim(d)[1]
  message(paste0("Observations lost due to missing predictor or covariates: ", drop_predictors))
  
  random = as.formula(paste0("~(1|", random_intercept, ")"))
  #---------------------------------------------
  # fit GAM ignoring spatial variation
  #---------------------------------------------
  print("gamfit")
  gam_fit = try(gamm4(formula, 
                      family = family,
                      data = d,
                      REML = T,
                      na.action = na.exclude, 
                      random = random
  ))
  
  
  #---------------------------------------------
  # confirm the model converged
  #---------------------------------------------
  print(paste0("Non-spatial gam fit: ", formula[3]))
  if(length(gam_fit) > 1 & !class(gam_fit) %in% c("try-error", "try-warning", "try-message")){
    print(summary(gam_fit$gam))
  } else {
    print("here2")
    print(gam_fit)
    #stop("Non-spatial GAM failed to converge")
    message = "Non-spatial GAM failed to converge"
    return(list(gam_fit = NA,
                spatial_gam_fit = NA,
                spatial_autocorr = NA,
                spatial_autocorr_coeff = NA,
                residual_spatial_autocorr = NA,
                residual_spatial_autocorr_coeff = NA,
                y = y, 
                a = a, 
                w_included = keep_w,
                interaction_term = interaction_term,
                formula = formula,
                family = family,
                model_input_data = d, 
                drop_outcome = drop_outcome, 
                drop_predictors = drop_predictors,
                message = message))
  }
  #---------------------------------------------
  # check for residual spatial autocorrelation
  #---------------------------------------------
  # print(gam_fit$gam)
  check_spatial_autocorr <- check_autocorr(gam_fit$gam, lldata = d[,c("qgpslong", "qgpslat")], nbc = 4)
  
  #---------------------------------------------
  # if residual spatial autocorrelation is present,
  # proceed with spatial GAM
  #---------------------------------------------
  if(!check_spatial_autocorr$spatial_autocorr){
    return(list(gam_fit = gam_fit,
                spatial_gam_fit = NA,
                spatial_autocorr = check_spatial_autocorr$spatial_autocorr,
                spatial_autocorr_coeff = check_spatial_autocorr$correlograms,
                residual_spatial_autocorr = NA,
                residual_spatial_autocorr_coeff = NA,
                y = y, 
                a = a, 
                w_included = keep_w,
                interaction_term = interaction_term,
                formula = formula,
                family = family,
                model_input_data = d, 
                drop_outcome = drop_outcome, 
                drop_predictors = drop_predictors,
                message = "Non-spatial gam fit successfully, no spatial autocorrelation detected so no spatial gam fit"))
  }else{
    
    #---------------------------------------------
    # If g is a factor then s(g,bs="re") produces a
    # random coefficient for each level of g, with
    # the random coefficients all modelled as i.i.d. normal.
    #---------------------------------------------
    
    
    spatial_gam_fit = try(gamm4(formula_spatial, 
                                family = family,
                                data = d,
                                REML = T,
                                na.action = na.exclude, 
                                random = random
    ))
    
   
    #---------------------------------------------
    # confirm the model converged
    #---------------------------------------------
    print("Spatial gam fit")
    # summary(spatial_gam_fit$gam)
    if(length(spatial_gam_fit) > 1 & !class(spatial_gam_fit) %in% c("try-error", "try-warning", "try-message")){
      print(summary(spatial_gam_fit$gam))
    } else {
      print(spatial_gam_fit)
      stop("Spatial GAM failed to converge")
    }
    
    # compare AIC for model with and without spatial smoother
    
    if (AIC(gam_fit$mer) < AIC(spatial_gam_fit$mer)) {
      print("Even though residual spatial autocorrelation was present, spatial GAM has larger AIC than non-spatial GAM. Check model.")
    }
    
    #---------------------------------------------
    # check for residual spatial autocorrelation 
    #   **after** adjusting for lat and long
    #---------------------------------------------
    check_residual_spatial_autocorr <- check_autocorr(spatial_gam_fit$gam, lldata = d[,c("qgpslong", "qgpslat")], nbc = 4)
    
    return(list(gam_fit = gam_fit,
                spatial_gam_fit = spatial_gam_fit,
                spatial_autocorr = check_spatial_autocorr$spatial_autocorr,
                spatial_autocorr_coeff = check_spatial_autocorr$correlograms,                
                residual_spatial_autocorr = check_residual_spatial_autocorr$spatial_autocorr,
                residual_spatial_autocorr_coeff = check_residual_spatial_autocorr$correlograms,
                y = y, 
                a = a, 
                w_included = c(keep_w, "qgpslong", "qgpslat"),
                interaction_term = interaction_term,
                formula = formula_spatial,
                family = family,
                model_input_data = d, 
                drop_outcome = drop_outcome, 
                drop_predictors = drop_predictors,
                message = "Spatial gam fit successfully"))
    # end of spatial autocorrelation if statement
  }
  
}

###############################################
# check_gam_fit
###############################################

# Documentation:    check_gam_fit
# Usage:            check_gam_fit(outcome, risk_factor, output_type, analysis, plot, results_dir)
# Description:      check model fit for a given outcome and risk factor pair, 
#                   provides output as either a .pdf file or inline plots/summary 
# Args/Options: 
#    outcome          character argument specifying the outcome
#    risk_factor      character argument specifying the risk factor
#    interaction_term character argument specifying the interaction term
#    output_type      character argument specifying the output type, accepted values are "pdf" (Default) or "in-line"
#    analysis         character argument specifying the outcome category (sth, diarrhea, giardia, pathogens), 
#                       adjusted/unadjusted, and model type (binary, continuous, count) with each element separated by a  
#                       hyphen (e.g. "sth-adjusted-binary"). NOTE: the diarrhea dataset can have -0 or -1 appended 
#                       to the end to signify which dataset, control (0) or intervention (1), and effect modification 
#                       analyses can also have -EM appended to the end.
#    plot             character argument specifying whether to plot only the risk factor (Default) or 
#                       risk factor and all adjustment covariates.
#    results_dir      character argument specifying the output directory, only used if outcome_type == "pdf"
# Returns:          The gam model fit summary, plot, and gam.check results, either as a pdf file or in-line
# Prints:           The gam model fit summary, plot, and gam.check results; nothing if output_type == "pdf"

check_gam_fit <- function(outcome, 
                          risk_factor, 
                          interaction_term = NULL,
                          output_type = c("pdf", "in-line"), 
                          analysis = NULL,
                          plot = c("only risk factor", "all"),
                          results_dir = results_path) {
  output_dir = paste0(results_dir, "gam_check/", analysis)
  if (output_type == "pdf") {
    render(input = gam_check_template_rmd,
           params = list(outcome = outcome, risk_factor = risk_factor, interaction_term = interaction_term, analysis = analysis, plot = plot),
           output_file = ifelse(is.null(interaction_term), paste0("gam_check_", outcome, "_", risk_factor, ".pdf"), paste0("gam_check_", outcome, "_", risk_factor, "_by_", interaction_term, ".pdf")),
           output_dir = output_dir)
  } else {
    if (is.null(interaction_term)) {
      x = readRDS(paste0(here::here(), 
                         "/results/gam_outputs/", analysis, "/gam_",
                         outcome, 
                         "_",
                         risk_factor,
                         ".RDS"))
    } else {
      x = readRDS(paste0(here::here(), 
                         "/results/gam_outputs/", analysis, "/gam_",
                         outcome, 
                         "_",
                         risk_factor,
                         "_by_", 
                         interaction_term,
                         ".RDS"))
      
      ################################################
      #### JADE, for testing, you can use this call
      ################################################
      # x = box_read(file_id = 952170420135) ## This is for an interaction model that was fit with formula: diar7d ~ t2(ppt_week_sum_1weeklag_C, absmintemp_30_days_C) + s(qgpslong, qgpslat, bs = "tp") 
    }
    
    message(paste0("Observations lost due to missing outcome: ", x$drop_outcome))
    message(paste0("Observations lost due to missing predictor or covariates: ", x$drop_predictors))
    if (is.na(x$spatial_gam_fit)) {
      fit.gam = x$gam_fit$gam
      fit.mer = x$gam_fit$mer
      message("Spatial GAM not fit, plotting non-spatial GAM results")
    } else {
      fit.gam = x$spatial_gam_fit$gam
      fit.mer = x$spatial_gam_fit$mer
    }
    print(summary(fit.mer))
    print(summary(fit.gam))
    simulationOut = simulateResiduals(fittedModel = fit.mer)
    plot(simulationOut)
    
    res2 = recalculateResiduals(simulationOut, group = x$model_input_data$dataid)
    plot(res2)

    
    if(plot == "only risk factor") {
      plot(x = fit.gam, trans = plogis, shift = coef(fit.gam)[1],
           select = 1,
           seWithMean = TRUE, residuals = T, pch = 1, cex = 1)#,
      # ylab = paste("Prevalence of", outcome))
    } else {
      plot(x = fit.gam, trans = plogis, shift = coef(fit.gam)[1],
           seWithMean = TRUE, residuals = TRUE, pch = 1, cex = 1)#,
      # ylab = paste("Prevalence of", outcome))
    }
    if (x$spatial_autocorr) { # Only run this chunk if spatial gam was used
      vis.gam(x = fit.gam,
              view = c("qgpslong", "qgpslat"),
              plot.type = "persp",
              se = 2)
      vis.gam(x = fit.gam,
              view = c("qgpslong", "qgpslat"),
              plot.type = "contour")
    }
 
  }
}



###############################################
# predict_gam
###############################################

# Documentation:     predict_gam
# Usage:             predict_gam(fit, data, riskfactorName)
# Description:       get marginal predicted probabilities from 
#                    logistic model for the observed range 
#                    of a given risk factor 
# Args/Options: 
#    fit             glm or gam model fit object
#    data            dataset used to fit model
#    riskfactorName  name of risk factor of interest, as character string
# Returns:           data frame with a column for the risk factor values
#                    over which the prediction was done and the predicted
#                    probabilities from the input model
# Prints:            nothing
predict_gam = function(fit, data, riskfactorName){
  
  
  wrapr::let(
    alias=list(riskfactor = riskfactorName),
    expr={
      
      # create vector of observed unique values
      riskfactor_vals = seq(min(data$riskfactor, na.rm = T),
                            max(data$riskfactor, na.rm = T),
                            0.01)
      
      pred_wrapper = function(x){
        mean(predict(fit, 
                     newdata = data %>% mutate(riskfactor = x), 
                     type= "response"))
      }
      
      preds = map_dbl(riskfactor_vals, 
                      pred_wrapper)
      
      df = data.frame(riskfactor = riskfactor_vals,
                      pred = preds)
      
    })
  
  return(df)
}




#----------------------------------
# adapted from: https://github.com/ben-arnold/mbita-schisto

# simultaneous CIs for GAMs
# estimated by resampling the 
# Bayesian posterior estimates of
# the variance-covariance matrix
# assuming that it is multivariate normal
# the function below also estimates 
# the unconditional variance-covariance
# matrix, Vb=vcov(x,unconditional=TRUE), 
# which allows for undertainty in the actual
# estimated mean as well 
# (Marra & Wood 2012 Scandinavian Journal of Statistics, 
#  Vol. 39: 53-74, 2012, doi: 10.1111/j.1467-9469.2011.00760.x )
# simultaneous CIs provide much better coverage than pointwise CIs
# see: http://www.fromthebottomoftheheap.net/2016/12/15/simultaneous-interval-revisited/
#
# @param       m : GAM model fit object from mgcv gam()
# @param newdata : data.frame on which to make predictions from the model fit.
#                  Must include all variables used to fit the model m
# @param  nreps  : number of replications to sample from the Bayesian posterior
#
# @param exclude_terms : model terms to exclude from prediction. Term names should be given 
#                        as they appear in the model summary (for example, "s(x0,x1)").
# @returns : gamCI returns a data.frame identical to newdata with 6 new variables:
#            NOTE: ALL PREDICTIONS ARE ON THE SCALE OF LINEAR PREDICTIONS FROM THE MODEL 
#                 (i.e., log-odds for logit model)
#            fit    : marginal spline prediction for each observation
#            se_fit : approximate standard error of the spline prediction
#            uprP   : upper pointwise 95% confidence interval
#            lwrP   : lower pointwise 95% confidence interval
#            uprS   : upper simultaneous 95% confidence interval
#            lwrS   : lower simultaneous 95% confidence interval
#----------------------------------
gamCI <- function(m,newdata,exclude_terms=NULL, nreps=10000) {
  require(mgcv)
  require(dplyr)
  Vb <- vcov(m,unconditional = TRUE)
  if(is.null(exclude_terms)) pred <- predict(m, newdata, se.fit = TRUE)
  if(!is.null(exclude_terms)) pred <- predict(m, newdata, exclude = exclude_terms, se.fit = TRUE)
  
  fit <- pred$fit
  se.fit <- pred$se.fit
  BUdiff <- MASS::mvrnorm(n=nreps, mu = rep(0, nrow(Vb)), Sigma = Vb)
  Cg <- predict(m, newdata, type = "lpmatrix")
  simDev <- Cg %*% t(BUdiff)
  absDev <- abs(sweep(simDev, 1, se.fit, FUN = "/"))
  masd <- apply(absDev, 2L, max)
  crit <- quantile(masd, prob = 0.95, type = 8)
  pred <- data.frame(newdata,fit=pred$fit,se_fit=pred$se.fit)
  pred <- mutate(pred,
                 uprP = fit + (2 * se.fit),
                 lwrP = fit - (2 * se.fit),
                 uprS = fit + (crit * se.fit),
                 lwrS = fit - (crit * se.fit)
  )
  return(pred)
}




###############################################
# fit_mgcv_gam_predict_contrasts function
# function to fit mgcv::gam; we will do prediciton using cluster-level bootstrap resampling

###############################################
# Documentation:    fit_mgcv_gam_predict_contrasts
# Usage:            fit_mgcv_gam_predict_contrasts(formula = as.formula(),
#                                                  family = "binomial",
#                                                  a = "pop_density_C",
#                                                  data = d_diarr %>%
#                                                           filter(season == "dry"),
#                                                  season = "dry",
#                                                  percentiles_df = risk_factor_mean_values)

# Description:      Function that uses mgcv::gam function to predict counterfactual outcomes for all children 
#                   using the 20th and 80th percentile risk factor (a) values and then calculates the risk 
#                   difference and risk ratio for having an 80th percentile value compared to a 20th percentile.
#                   The function executes these steps:
#                       1) generates a new model fit using mgcv::gam using the given formula and data, 
#                       2) predicts outcome for each observation using the 20th percentile value of the exposure, 
#                       3) predicts outcome for each observation using the 80th percentile value of the exposure,
#                       4) calculates the risk difference (RD) and risk ratio (RR) for 80th vs 20th percentile exposures
#                         
# Args/Options: 
#   formula:        Model formula (for consistency between this analysis and the main analysis, we want to provide the formula  
#                   from the original gamm4 model that was fit previously)
#   a:              Character argument that specifies the exposure variable (risk factor)
#   family:         The family object specifying the distribution and link to use in fitting the gam. Acceptable responses  
#                   are currently "binomial" and "gaussian", although in practice any allowed by family.mgcv could be possible.
#   data:           Data.frame containing all variables needed for analysis (for consistency between this analysis and the 
#                   main analysis, we want to provide the data.frame used as model input for the original gamm4 model, which 
#                   can be obtained from the model fit object output $model_input_data element. This data.frame has already
#                   screened covariates with likelihood ratio test, checked for sparsity, etc.)
#   season:         Character argument specifying which season to do the analysis on. Acceptable responses are "all", "rainy", 
#                   or "dry." Defaults to "all."
#   predict_low: Data.frame containing the mean, 20th and 80th percentile values for all risk factors in the dry and rainy season.
# 
# Returns: A list with the following elements:
#   RD:      risk difference
#   RR:      risk ratio
#   OR:      odds ratio



fit_mgcv_gam_predict_contrasts <- function(formula, a, family, data, predict_low, predict_high) { #season = c("all", "rainy", "dry"),
  new_fit = gam(formula,
                family = family,
                data = data,
                method = "REML",
                select = TRUE,
                na.action=na.exclude)
  
  # pct20 = quantile(data[,a], probs = 0.2, na.rm= T)
  # pct80 = quantile(data[,a], probs = 0.8, na.rm= T)
  
  predict.low = predict(new_fit, 
                        newdata = data %>% mutate(!!sym(a) := predict_low),
                        se.fit = TRUE)
  predict.high = predict(new_fit, 
                         newdata = data %>% mutate(!!sym(a) := predict_high),
                         se.fit = TRUE)
  
  odds.low = exp(predict.low$fit)
  odds.high = exp(predict.high$fit)
  prev.low = odds.low/(1 + odds.low)
  prev.high = odds.high/(1 + odds.high)
  
  # Calculate Risk Difference, Risk Ratio, and Odds ratio for binary outcomes
  RD = mean(prev.high) - mean(prev.low) 
  RR = mean(prev.high)/mean(prev.low)
  OR = mean(odds.high)/mean(odds.low)
  # prev.df = data.frame(prev.low = prev.low,
  #                      prev.high = prev.high,
  #                      odds.low = odds.low,
  #                      odds.high = odds.high)
  out = list(RD = RD, 
             RR = RR,
             OR = OR)#,
  # predicted.df = prev.df)
  ##############################################################################
  ## Will need to add code later to handle non-binary outcomes
  ##############################################################################
  
  return(out)
}



##############################################
##############################################
# Documentation: plot_gam_int_hist
# Usage: plot_gam_int_hist(file_path, file_name, model_subdir, outcome_subdir, 
# risk_factor1, risk_factor1_label, 
# risk_factor2, risk_factor2_label, outcome, outcome_label,
# scale, x_lower, x_upper)
# Description: Make a plot with the smooth gam fits for interactions for 
# continuous x categorical risk factor levels and include a histogram
# for the continuous risk factor 

# Args/Options:
# file_path:          file path for interaction output, as a string
# file_name:          file name with interaction output, as a string
# model_subdir:       model subdirectory for saved figures, as a string 
# risk_factor1:       variable name for continuous risk factor, as a string 
# risk_factor1_label: label for continuous risk factor, as a string 
# risk_factor2:       variable name for categorical risk factor, as a string 
# risk_factor2_label: label for categorical risk factor, as a string 
# outcome:            variable name for outcome, as a string 
# outcome_label:      label for outcome, as a string 
# scale:              number to multiply outcome by (e.g., 100 for prevalence)
# x_lower:
# x_upper:

# Returns: saves png with figure in appropriate directory 
# Output: none

plot_gam_int_hist <- function(file_path, file_name, model_subdir,
                              risk_factor1, risk_factor1_label, 
                              risk_factor2, risk_factor2_label, outcome, outcome_label,
                              scale, x_lower, x_upper) {
  #------------------------------------------------------
  # load in output
  #------------------------------------------------------
  x = readRDS(paste0(file_path, file_name))
  data = x$model_input_data
  
  #------------------------------------------------------
  # check whether model was fit 
  # if not, exit function
  #------------------------------------------------------
  if(!(all(is.na(x$gam_fit)) & all(is.na(x$spatial_gam_fit)))){

  
  #------------------------------------------------------
  # Outcomes with no spatial autocorrelation will use a separate gam fit object
  #------------------------------------------------------
  if (!all(is.na(x$spatial_gam_fit))) {
    gamfit = x$spatial_gam_fit$gam
  } else {
    gamfit = x$gam_fit$gam
  }
  
  
  #------------------------------------------------------
  # Create data frame with all unique observed values of risk factors
  #------------------------------------------------------
  ## Add means values back so centered risk factors to put back on their original scale
  rf_mean = box_search(paste0("\"","risk_factor_mean_values.RDS", "\"")) %>%
    box_read() %>%
    filter(variable == gsub("_C", "", risk_factor1)) %>%
    pull(mean_value)
  
  # number of levels of second risk factor 
  rf2_levels <- levels(data %>% pull(!!sym(risk_factor2)))
  
  newd_list <- list()
  fitp_ci_list <-list()
  
  # for each level of risk factor 2, create dataset filtered to those
  # values, then fit CIs
  for(i in 1:length(rf2_levels)){
    # filter to risk factor 2 level i
    newd_list[[i]] <- data %>% filter(!sym(risk_factor2) == rf2_levels[i])
    # Obtain simultaneous CIs
    fitp_ci_list[[i]] <- gamCI(m=gamfit,newdata=newd_list[[i]],exclude_terms = "s(qgpslong,qgpslat)", nreps=10000)
    
    #------------------------------------------------------
    # Obtain model family
    #------------------------------------------------------
    family <- gamfit$family[[1]]
    
    #------------------------------------------------------
    # Calculate y-value transformations and y-axis labels based on model family
    #------------------------------------------------------
    if(family == "binomial") {
      fitp_ci_list[[i]] <- fitp_ci_list[[i]] %>% 
        mutate(fit = exp(fit)/(1+exp(fit)) * scale,
               lwrS = exp(lwrS)/(1+exp(lwrS)) * scale,
               uprS = exp(uprS)/(1+exp(uprS)) * scale)
      ylab_name <- "Prevalence (%)"
    }
    
    if(family == "poisson") {
      fitp_ci_list[[i]] <- fitp_ci_list[[i]] %>% 
        mutate(fit = exp(fit) * scale,
               lwrS = exp(lwrS) * scale,
               uprS = exp(uprS) * scale)
      ylab_name <- "Number"
    }  
    
    # ### NOT WORKING next two chunks
    # fitp_ci_list[[i]] <- fitp_ci_list[[i]] %>% 
    #   mutate(!!sym(risk_factor1) := !!sym(risk_factor1) + rf_mean) %>%
    #   filter(!!sym(risk_factor1) > x_lower & !!sym(risk_factor1) < x_upper)
    # 
    # data <- data %>% 
    #   mutate(!!sym(risk_factor1) := !!sym(risk_factor1) + rf_mean) %>%
    #   filter(!!sym(risk_factor1) > x_lower & !!sym(risk_factor1) < x_upper)
    
    # xlim <- if (!is.null(x_lower) & !is.null(x_upper)) {
    #   xlim(x_lower, x_upper)
    # }  else if (!is.null(x_lower) & is.null(x_upper)) {
    #   xlim(x_lower, max(fitp_ci[[risk_factor]]))
    # } else if (is.null(x_lower) & !is.null(x_upper)) {
    #   xlim(min(fitp_ci[[risk_factor]]), x_upper)
    # } else if (is.null(x_lower) & is.null(x_upper)) {
    #   xlim(min(fitp_ci[[risk_factor]]), max(fitp_ci[[risk_factor]]))
    # }
    
  }
  
  fitp_ci <- bind_rows(fitp_ci_list)
  
  #------------------------------------------------------
  # Make plots
  #------------------------------------------------------
  plot_smooth <- ggplot(fitp_ci, aes(x = !!sym(risk_factor1), y = fit, 
                                     group = !!sym(risk_factor2))) + 
    geom_ribbon(aes(ymin = lwrS, ymax= uprS, fill = !!sym(risk_factor2)), alpha =0.25) +
    geom_line(aes(col = !!sym(risk_factor2)))  +
    scale_color_discrete(name=risk_factor2_label) + 
    scale_fill_discrete(name=risk_factor2_label) + 
    ylab(ylab_name) +
    # xlim +
    theme_bw() +
    theme(
      axis.title.x = element_blank(),
      legend.position = "bottom"
    )
  
  plot_histogram <- ggplot(data, aes(x = !!sym(risk_factor1))) +
    geom_histogram(bins = 100) +
    xlab(risk_factor1_label) + 
    ylab("Count") + 
    # xlim +
    theme_bw() +
    theme(
      panel.grid.major = element_blank(), 
      panel.grid.minor = element_blank()
    )
  
  combined_plot <- grid.arrange(plot_smooth, plot_histogram, heights = c(4,2),
                                top = textGrob(paste0(outcome_label,"\n",
                                                      "(", risk_factor1_label, " x\n",
                                                      risk_factor2_label,")"), gp=gpar(fontsize=14)))
  
  
  #------------------------------------------------------
  # Save plots
  #------------------------------------------------------
  figure_directory_box <- paste0(local_root_path, local_box_path,
         "figures-continuous-risk-factors-intrxn/")
  
  save_name <- paste0(figure_directory_box, model_subdir, 
                      outcome,"_",risk_factor1, "_", risk_factor2, ".png")
  ggsave(save_name, combined_plot, bg = "white",
         width = 14, height = 11.5, units = "cm")

  
  } else{
    # if model wasn't fit
    print(paste0("No model output for ", outcome, 
                 " - ", risk_factor1, " - ", risk_factor2))
  }
}



crosspred_inc <- function (basis, model = NULL, coef = NULL, vcov = NULL, model.link = NULL, 
                           at = NULL, from = NULL, to = NULL, by = NULL, lag, bylag = 1, 
                           ci.level = 0.95, cumul = FALSE) 
{
  type <- if (any(class(basis) %in% "crossbasis")) 
    "cb"
  else if (any(class(basis) %in% "onebasis")) 
    "one"
  else "gam"
  errormes <- "arguments 'basis' and 'model' not consistent. See help(crosspred)"
  if (type == "gam") {
    if (!is.character(basis) || length(basis) > 1L) 
      stop(errormes)
    if (is.null(model) || !any(class(model) %in% "gam")) 
      stop(errormes)
    name <- basis
    sterms <- sapply(model$smooth, function(x) x$term[1])
    if (name %in% sterms) 
      basis <- model$smooth[[which(sterms == name)[1]]]
    else stop(errormes)
    if (length(which(sterms == name)) > 1) 
      warning(paste(name, "included in multiple smoothers, only the first one taken"))
    if (!"cb.smooth" %in% class(basis) && basis$dim > 1L) 
      stop("predictions not provided for multi-dimensional smoothers other than 'cb'")
  }
  else name <- deparse(substitute(basis))
  origlag <- switch(type, cb = attr(basis, "lag"), one = c(0, 
                                                           0), gam = if (is.null(basis$lag)) c(0, 0) else basis$lag)
  lag <- if (missing(lag)) 
    origlag
  else dlnm:::mklag(lag)
  if (!all(lag == origlag) && cumul) 
    stop("cumulative prediction not allowed for lag sub-period")
  lagfun <- switch(type, cb = attr(basis, "arglag")$fun, one = NULL, 
                   gam = if (basis$dim == 1L) NULL else basis$margin[[2]]$fun)
  if (bylag != 1L && !is.null(lagfun) && lagfun == "integer") 
    stop("prediction for non-integer lags not allowed for type 'integer'")
  if (is.null(model) && (is.null(coef) || is.null(vcov))) 
    stop("At least 'model' or 'coef'-'vcov' must be provided")
  if (!is.numeric(ci.level) || ci.level >= 1 || ci.level <= 
      0) 
    stop("'ci.level' must be numeric and between 0 and 1")
  cond <- if (type == "gam") 
    with(basis, first.para:last.para)
  else if (ncol(basis) == 1L) 
    name
  else if (type == "one") 
    paste(name, "[[:print:]]*b[0-9]{1,2}", sep = "")
  else paste(name, "[[:print:]]*v[0-9]{1,2}\\.l[0-9]{1,2}", 
             sep = "")
  if (!is.null(model)) {
    model.class <- class(model)
    coef <- dlnm:::getcoef(model, model.class)
    vcov <- dlnm:::getvcov(model, model.class)
    indcoef <- if (type == "gam") 
      cond
    else grep(cond, names(coef))
    indvcov <- if (type == "gam") 
      cond
    else grep(cond, rownames(vcov))
    coef <- coef[indcoef]
    vcov <- vcov[indvcov, indvcov, drop = FALSE]
    model.link <- dlnm:::getlink(model, model.class, model.link)
  }
  else model.class <- NA
  npar <- if (type == "gam") 
    length(indcoef)
  else ncol(basis)
  if (length(coef) != npar || length(coef) != dim(vcov)[1] || 
      any(is.na(coef)) || any(is.na(vcov))) 
    stop("coef/vcov not consistent with basis matrix. See help(crosspred)")
  range <- if (type == "gam") 
    range(model$model[[basis$term[1]]])
  else attr(basis, "range")
  at <- dlnm:::mkat(at, from, to, by, range, lag, bylag)
  predvar <- if (is.matrix(at)) 
    rownames(at)
  else at
  predlag <- dlnm:::seqlag(lag, bylag)
  Xpred <- dlnm:::mkXpred(type, basis, at, predvar, predlag)
  matfit <- matrix(Xpred %*% coef, length(predvar), length(predlag))
  matse <- matrix(sqrt(pmax(0, rowSums((Xpred %*% vcov) * Xpred))), 
                  length(predvar), length(predlag))
  rownames(matfit) <- rownames(matse) <- predvar
  colnames(matfit) <- colnames(matse) <- outer("lag", predlag, 
                                               paste, sep = "")
  predlag <- dlnm:::seqlag(lag)
  Xpred <- dlnm:::mkXpred(type, basis, at, predvar, predlag)
  Xpredall <- 0
  if (cumul) {
    cumfit <- cumse <- matrix(0, length(predvar), length(predlag))
  }
  for (i in seq(length(predlag))) {
    ind <- seq(length(predvar)) + length(predvar) * (i - 
                                                       1)
    Xpredall <- Xpredall + Xpred[ind, , drop = FALSE]
    if (cumul) {
      cumfit[, i] <- Xpredall %*% coef
      cumse[, i] <- sqrt(pmax(0, rowSums((Xpredall %*% 
                                            vcov) * Xpredall)))
    }
  }
  allfit <- as.vector(Xpredall %*% coef)
  allse <- sqrt(pmax(0, rowSums((Xpredall %*% vcov) * Xpredall)))
  names(allfit) <- names(allse) <- predvar
  if (cumul) {
    rownames(cumfit) <- rownames(cumse) <- predvar
    colnames(cumfit) <- colnames(cumse) <- outer("lag", dlnm:::seqlag(lag), 
                                                 paste, sep = "")
  }
  list <- list(predvar = predvar)
  list <- c(list, list(lag = lag, bylag = bylag, coefficients = coef, 
                       vcov = vcov, matfit = matfit, matse = matse, allfit = allfit, 
                       allse = allse))
  if (cumul) 
    list <- c(list, list(cumfit = cumfit, cumse = cumse))
  z <- qnorm(1 - (1 - ci.level)/2)
  list$matlow <- matfit - z * matse
  list$mathigh <- matfit + z * matse
  list$alllow <- allfit - z * allse
  names(list$alllow) <- names(allfit)
  list$allhigh <- allfit + z * allse
  names(list$allhigh) <- names(allfit)
  if (cumul) {
    list$cumlow <- cumfit - z * cumse
    list$cumhigh <- cumfit + z * cumse
  }
  list$ci.level <- ci.level
  list$model.class <- model.class
  list$model.link <- model.link
  class(list) <- "crosspred"
  return(list)
}

fit_dlnm <- function(data, predictor, predictor_type, var_value, lag_value, covariate, predictor_name, color_code, plot_tags, outcome = "malaria_wk_inc", offset_var = "population",
                     sens_analysis) {
  
  # Compute knots for var and lag
  predictor_varknots <- equalknots(data[[predictor]], fun = "ns", df = 2)
  covariate_varknots <- equalknots(data[[covariate]], fun = "ns", df = 2)
  lagknots <- equalknots(x = c(2, 16), nk = 1)
  
  # Create crossbasis for predictor
  cb_predictor <- crossbasis(data[[predictor]], lag = c(2, 16), bylag = 1,
                             argvar = list(fun = "ns", knots = predictor_varknots),
                             arglag = list(knots = lagknots))
  
  # Create crossbasis for covariate
  cb_covariate <- crossbasis(data[[covariate]], lag = c(2, 16), bylag = 1,
                             argvar = list(fun = "ns", knots = covariate_varknots),
                             arglag = list(knots = lagknots))
  
  if (predictor == "temp_wk_max" & sens_analysis == 1) {
    
    cb_covariate_extra <- crossbasis(data[["temp_wk_min"]], lag = c(2, 16), bylag = 1,
                                     argvar = list(fun = "ns", 
                                                   knots = equalknots(data[["temp_wk_min"]], 
                                                                      fun = "ns", df = 2)),
                                     arglag = list(knots = lagknots))
    
    formula <- as.formula(paste(outcome, "~ cb_predictor + cb_covariate + cb_covariate_extra + year + offset(log(", offset_var, "))"))
    
  } else if (predictor == "temp_wk_max" & sens_analysis == 2) {
    cb_covariate <- crossbasis(data[[covariate]], lag = c(6, 20), bylag = 1,
                               argvar = list(fun = "ns", knots = covariate_varknots),
                               arglag = list(knots = equalknots(x = c(6, 20), nk = 2)))
    cb_covariate_extra <- crossbasis(data[["temp_wk_min"]], lag = c(2, 16), bylag = 1,
                                     argvar = list(fun = "ns", 
                                                   knots = equalknots(data[["temp_wk_min"]], 
                                                                      fun = "ns", df = 2)),
                                     arglag = list(knots = lagknots))
    
    formula <- as.formula(paste(outcome, "~ cb_predictor + cb_covariate + year + offset(log(", offset_var, "))"))
    
  } else {
    formula <- as.formula(paste(outcome, "~ cb_predictor + cb_covariate + year + offset(log(", offset_var, "))"))
    
  }
  
  model <- gam(data = data, formula = formula, family = "poisson")
  
  if (predictor_type == "temp") {
    if (predictor == "temp_wk_max") {
      pred <- crosspred(cb_predictor, model, 
                        cen = min(data[[predictor]], na.rm = TRUE), by = 0.1, cumul = T)
    } else {
      pred <- crosspred(cb_predictor, model, 
                        cen = min(data[[predictor]], na.rm = TRUE), by = 0.1, cumul = T)
    }
  } else {
    pred <- crosspred(cb_predictor, model, cen = 0, by = 50, cumul = T)
  }
  
  plot_data <- data.frame(pred["predvar"], 
                          pred[["matRRfit"]]) %>% 
    pivot_longer(cols = starts_with("lag"),
                 names_to = "lag",
                 names_prefix = "lag",
                 values_to = "fit")
  
  lower_ci_data <- data.frame(pred["predvar"], 
                              pred["matRRlow"]) %>% 
    pivot_longer(cols = starts_with("matRRlow."),
                 names_to = "lag",
                 names_prefix = "matRRlow.lag",
                 values_to = "lower_ci")
  plot_data <- left_join(plot_data, lower_ci_data, by = c("predvar", "lag"))
  
  upper_ci_data <- data.frame(pred["predvar"], 
                              pred["matRRhigh"]) %>% 
    pivot_longer(cols = starts_with("matRRhigh."),
                 names_to = "lag",
                 names_prefix = "matRRhigh.lag",
                 values_to = "upper_ci")
  if (predictor == "temp_wk_max") {
    plot_data <- left_join(plot_data, upper_ci_data, by = c("predvar", "lag")) %>% 
      mutate(lag = as.numeric(lag)) %>% 
      mutate(predvar = round(predvar, 1))
  } else {
    plot_data <- left_join(plot_data, upper_ci_data, by = c("predvar", "lag")) %>% 
      mutate(lag = as.numeric(lag)) %>% 
      mutate(predvar = round(predvar, 1)) %>% 
      filter(predvar > quantile(data[[predictor]], p = 0.05) &
               predvar < quantile(data[[predictor]], p = 0.95)) # truncate limits to data range!
  }
  
  plot_list <- list()
  result_data <- data.frame()
  this_plot_data <- plot_data %>% filter(predvar == var_value)
  
  if (predictor_name == "Total Precipitation (mm)") {
    point_title = paste0("Association at Median ", predictor_name)
  } else {
    point_title = paste0("Association at Median ", predictor_name)
  }
  
  plot <- ggplot(this_plot_data, aes(x = lag, y = fit)) +
    geom_errorbar(aes(ymin = lower_ci, ymax = upper_ci), 
                  width = 0.2, color = color_code) +
    geom_point(color = color_code) +
    geom_hline(yintercept = 1, linetype = "longdash", color = "black") +
    labs(title = point_title,
         tag = plot_tags[1],
         x = "Lag (weeks)", 
         y = "Incidence Ratio") +
    scale_y_continuous(trans = "log", 
                       breaks = scales::breaks_extended()) +
    theme_minimal() +
    theme(legend.position = "none",
          plot.title = element_text(size = 10, hjust = 0.5),
          plot.margin = margin(10, 10, 10, 10),
          plot.tag = element_text(size = 12, face = "bold", hjust = 0.1, vjust = 0.1),
          axis.line = element_line(color = "black", linewidth = 0.2),
          axis.title = element_text(size = 10),  # Smaller axis title size
          axis.text = element_text(size = 8))
  result_data <- rbind(result_data, this_plot_data)
  plot_list[[1]] <- plot
  
  this_plot_data <- plot_data %>% filter(lag == lag_value)
  lag_plot <- ggplot(this_plot_data, mapping = aes(x = predvar, y = fit)) +
    geom_ribbon(mapping = aes(ymin = lower_ci, ymax = upper_ci), 
                alpha = 0.2, color = NA, fill = color_code) +
    geom_line(color = color_code) +
    geom_hline(yintercept = 1, linetype = "longdash", color = "black") +
    labs(title = paste0("Association at ", lag_value, "-week lag"),
         tag = plot_tags[2],
         x = predictor_name, 
         y = "Incidence Ratio") +
    scale_y_continuous(trans = "log", 
                       breaks = scales::breaks_extended()) +
    theme_minimal() +
    theme(legend.position = "none",
          plot.title = element_text(size = 10, hjust = 0.5),
          plot.margin = margin(10, 10, 10, 10),
          plot.tag = element_text(size = 12, face = "bold", hjust = 0.1, vjust = 0.1),
          axis.line = element_line(color = "black", linewidth = 0.2),
          axis.title = element_text(size = 10),  # Smaller axis title size
          axis.text = element_text(size = 8))
  result_data <- rbind(result_data, this_plot_data)
  plot_list[[2]] <- lag_plot
  
  # create contour plot
  contour_plot <- ggplot(plot_data, aes(x = predvar, y = lag, z = fit)) +
    geom_raster(aes(fill = fit)) +
    geom_contour(color = "white", alpha = 0.5) +
    scale_fill_gradient2(
      low = "#4DAF4A",        
      mid = "white",       
      high = "#C51B7D",        
      midpoint = 1,        
      breaks = labeling::extended(m = 15, dmin = 0, dmax = max(plot_data$fit, na.rm = T)),  
      labels = scales::number_format(accuracy = 0.1)  # Format the labels
    ) +
    labs(title = "Predictor-Lag Surface",
         tag = plot_tags[3],
         x = predictor_name,
         y = "Lag (weeks)",
         fill = "Incidence\nRatio"
    ) +
    theme_minimal() +
    theme(plot.title = element_text(size = 10, hjust = 0.5),
          plot.margin = margin(10, 10, 10, 10),
          plot.tag = element_text(size = 12, face = "bold", hjust = 0.2, vjust = 0.1),
          axis.title = element_text(size = 10),  # Smaller axis title size
          axis.text = element_text(size = 8),
          legend.text = element_text(size = 10),
          legend.title = element_text(size = 10),
          legend.key.size = unit(1.0, "cm"),
          legend.key.width = unit(0.3, "cm"),
          panel.grid = element_blank())
  plot_list[[3]] <- contour_plot
  
  result_data <- result_data %>% mutate(Measure = predictor_name,
                                        Covariate = covariate,
                                        `Predictor Value` = predvar,
                                        `Lag` = lag,
                                        `Incidence Ratio (95% CI)` = paste0(sprintf("%.02f", fit), " (", sprintf("%.02f", lower_ci), " - ", sprintf("%.02f", upper_ci), ")")) %>% select(Measure, Covariate, Lag, `Predictor Value`, `Incidence Ratio (95% CI)`)
  
  plot_list[[4]] <- result_data
  print(summary(model))
  return(plot_list)
}
