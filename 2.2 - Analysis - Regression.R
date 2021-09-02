library(tidyverse)
library(gtsummary)

# LOADING DATASET

db <- 
  readRDS("brain.rds")

covs <- paste(c('age',
                'pmi',
                'area_cw',
                'apoe_allele_grp1',
                'cognitive_status_num',
                'sex'), collapse = "+")

form <- as.formula(paste("dens_avg ~ factor(braak06)", covs, sep = "+"))
form_lin <- as.formula(paste("dens_avg ~ braak06", covs, sep = "+"))
form_log <- as.formula(paste("log(dens_avg) ~ factor(braak06)", covs, sep = "+"))
form_log_lin <- as.formula(paste("log(dens_avg) ~ braak06", covs, sep = "+"))

label <- list(age ~ 'Age (years)',
              pmi ~ 'PMI (minutes)',
              area_cw ~ 'Area (10^6)',
              #apoe_allele_grp1 ~ 'APOE ALLELE',
              apoe_allele_grp2 ~ 'Apoe Allele',
              cognitive_status_num ~ 'Cognitive Status',
              sex ~ 'Sex',
              braak06 ~ "Braak6 Score")

##########################
# LINEAR MODEL - OUTCOME #
##########################

# UNIVARIATE ANALYSIS #

fit_lm_univ <- lm(dens_avg ~ factor(braak06), data = db); summary(fit_lm_univ) # univariate - categorical
fit_lm_univ2 <- lm(dens_avg ~ braak06, data = db); summary(fit_lm_univ2) # univariate - linear

#anova(fit_lm_univ2, fit_lm_univ) # fail to reject simpler
#hnp::hnp(fit_lm_univ2)

tab_lm_univ <- 
  tbl_regression(fit_lm_univ2) %>% 
  bold_labels() %>% 
  add_global_p()

# MULTIVARIATE ANALYSIS #

fit_lm_mult <- lm(form, data = db); summary(fit_lm_mult) # categorical effect
fit_lm_mult2 <- lm(form_lin, data = db); summary(fit_lm_mult2) # linear effect
#anova(fit_lm_mult2, fit_lm_mult) # reject linear effect

fit_lm_mult3 <- update(fit_lm_mult, ~.-apoe_allele_grp1 + apoe_allele_grp2); summary(fit_lm_mult3) # linear effect
#anova(fit_lm_mult3, fit_lm_mult) # fail to reject group2 allele

#hnp::hnp(fit_lm_mult3)

tab_lm_mult <- 
  tbl_regression(fit_lm_mult3) %>% 
  bold_labels() %>% 
  add_global_p()

##############################
# LINEAR MODEL - LOG OUTCOME #
##############################

# UNIVARIATE ANALYSIS

fit_lm_log_univ <- lm(log(dens_avg) ~ factor(braak06), data = db); summary(fit_lm_log_univ) # univariate - categorical
fit_lm_log_univ2 <- lm(log(dens_avg) ~ braak06, data = db); summary(fit_lm_log_univ2) # univariate - linear

#anova(fit_lm_log_univ2, fit_lm_log_univ) # fail to reject simpler
#hnp::hnp(fit_lm_log_univ2)

tab_lm_log_univ <- 
  tbl_regression(fit_lm_log_univ2) %>% 
  bold_labels()

# MULTIVARIATE ANALYSIS

fit_lm_log_mult <- lm(form_log, data = db); summary(fit_lm_log_mult) # categorical effect
fit_lm_log_mult2 <- lm(form_log_lin, data = db); summary(fit_lm_log_mult2) # linear effect
#anova(fit_lm_log_mult2, fit_lm_log_mult) # fail to reject linear effect

fit_lm_log_mult3 <- update(fit_lm_log_mult2, ~.-apoe_allele_grp1 + apoe_allele_grp2); summary(fit_lm_log_mult3) # linear effect
#anova(fit_lm_log_mult3, fit_lm_log_mult2) # fail to reject group2 allele

#hnp::hnp(fit_lm_log_mult3)

tab_lm_log_mult <- 
  tbl_regression(fit_lm_log_mult3) %>% 
  bold_labels() %>% 
  add_global_p()

#################################### 
# GENERALIZED LINEAR MODEL - GAMMA #
####################################

# UNIVARIATE ANALYSIS

fit_glm_univ <- glm(dens_avg ~ factor(braak06), 
                    data = db, 
                    family = Gamma(link = log)); summary(fit_glm_univ)
fit_glm_univ2 <- glm(dens_avg ~ braak06, 
                     data = db, 
                     family = Gamma(link = log)); summary(fit_glm_univ2)
#anova(fit_glm_univ2, fit_glm_univ, test = "LRT") # reject linear effect

#hnp::hnp(fit_glm_univ2)

tab_glm_univ <- 
  tbl_regression(fit_glm_univ2,
                  label = list(braak06 ~ 'Braak6 Score'),
                 exponentiate = TRUE) %>% 
  bold_labels()

# MULTIVARIATE ANALYSIS

fit_glm_mult <- glm(form, 
                    data = db, 
                    family = Gamma(link = log)); summary(fit_glm_mult) # categorical effect
fit_glm_mult2 <- glm(form_lin, 
                     data = db, 
                     family = Gamma(link = log)); summary(fit_glm_mult2) # linear effect
#anova(fit_glm_mult2, fit_glm_mult, test = "LRT") # fail to reject linear effect

fit_glm_mult3 <- update(fit_glm_mult2, ~.-apoe_allele_grp1 + apoe_allele_grp2); summary(fit_glm_mult3) # linear effect
#anova(fit_glm_mult3, fit_glm_mult2, test = "LRT") # fail to reject group2 allele

#hnp::hnp(fit_glm_mult3)

tab_glm_mult <- 
  tbl_regression(fit_glm_mult3,
                 label = label,
                 exponentiate = TRUE) %>% 
  bold_labels() %>% 
  add_global_p()

##############################################
# GENERALIZED LINEAR MODEL - GAMMA - INVERSE #
##############################################

# UNIVARIATE ANALYSIS

fit_glm_univ_INV <- glm(dens_avg ~ factor(braak06), 
                    data = db, 
                    family = Gamma); summary(fit_glm_univ_INV)
fit_glm_univ2_INV <- glm(dens_avg ~ braak06, 
                     data = db, 
                     family = Gamma); summary(fit_glm_univ2_INV)
#anova(fit_glm_univ2_INV, fit_glm_univ_INV, test = "LRT") # reject linear effect

#hnp::hnp(fit_glm_univ_INV)

tab_glm_univ_INV <- 
  tbl_regression(fit_glm_univ_INV) %>% 
  bold_labels() %>% 
  add_global_p()

# MULTIVARIATE ANALYSIS

fit_glm_mult_INV <- glm(form, 
                    data = db, 
                    family = Gamma); summary(fit_glm_mult_INV) # categorical effect
fit_glm_mult2_INV <- glm(form_lin, 
                     data = db, 
                     family = Gamma); summary(fit_glm_mult2_INV) # linear effect
#anova(fit_glm_mult2_INV, fit_glm_mult_INV, test = "LRT") # fail to reject linear effect

fit_glm_mult3_INV <- update(fit_glm_mult2_INV, ~.-apoe_allele_grp1 + apoe_allele_grp2); summary(fit_glm_mult3_INV) # linear effect
#anova(fit_glm_mult3_INV, fit_glm_mult2_INV, test = "LRT") # fail to reject group2 allele

#hnp::hnp(fit_glm_mult3_INV)

tab_glm_mult_INV <- 
  tbl_regression(fit_glm_mult3_INV,
                 label = label) %>% 
  bold_labels() %>% 
  add_global_p()


