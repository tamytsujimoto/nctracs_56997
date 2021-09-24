library(tidyverse)
library(gtsummary)
library(multcomp)
library(gt)

# LOADING DATASET

db <- 
  readRDS("brain.rds") %>% 
  mutate(braak06_grp = ifelse(braak06 > 3, 1, 0),
         braak06_fac = factor(braak06))

covs <- paste(c('age',
                'pmi',
                #'area_cw',
                'apoe_allele_grp2',
                'cognitive_status_num',
                'sex'), collapse = "+")

form1 <- paste("braak06_fac", covs, sep = "+")
form2 <- paste("braak06 + I(braak06-3):braak06_grp", covs, sep = "+")
form3 <- paste("braak06", covs, sep = "+")

##############################
# LINEAR MODEL - LOG OUTCOME #
##############################

# UNIVARIATE #
fit_univ_lm1 <- lm(log(dens_avg) ~ braak06_fac, data = db); 
fit_univ_lm2 <- lm(log(dens_avg) ~ braak06 + I(braak06-3):braak06_grp, data = db); 
fit_univ_lm3 <- lm(log(dens_avg) ~ braak06, data = db); 

# anova(fit_univ_lm2, fit_univ_lm1) # Reduced model
# anova(fit_univ_lm3, fit_univ_lm2) # Reduced model
# 
# hnp::hnp(fit_univ_lm1)
# hnp::hnp(fit_univ_lm2)
# hnp::hnp(fit_univ_lm3)

# MULTIVARIATE #
fit_mult_lm1 <- lm(as.formula(paste0("log(dens_avg) ~", form1)), data = db); summary(fit_mult_lm1)
fit_mult_lm2 <- lm(as.formula(paste0("log(dens_avg) ~", form2)), data = db); summary(fit_mult_lm2)
fit_mult_lm3 <- lm(as.formula(paste0("log(dens_avg) ~", form3)), data = db); summary(fit_mult_lm3)

# anova(fit_mult_lm2, fit_mult_lm1) # Reduced model
# anova(fit_mult_lm3, fit_mult_lm2) # Reduced model
# 
# hnp::hnp(fit_mult_lm1)
# hnp::hnp(fit_mult_lm2)
# hnp::hnp(fit_mult_lm3)

#################################### 
# GENERALIZED LINEAR MODEL - GAMMA #
####################################

# UNIVARIATE #
fit_univ_glm1 <- glm(dens_avg ~ factor(braak06), data = db, family = Gamma); summary(fit_univ_glm1)
fit_univ_glm2 <- glm(dens_avg ~ braak06 + I(braak06-3):braak06_grp, data = db, family = Gamma); summary(fit_univ_glm2)
fit_univ_glm3 <- glm(dens_avg ~ braak06, data = db, family = Gamma); summary(fit_univ_glm3)

# anova(fit_univ_glm2, fit_univ_glm1, test = "LRT") # Reduced model
# anova(fit_univ_glm3, fit_univ_glm2, test = "LRT") # Reduced model
# 
# hnp::hnp(fit_univ_glm1)
# hnp::hnp(fit_univ_glm2)
# hnp::hnp(fit_univ_glm3)

# MULTIVARIATE #
fit_mult_glm1 <- glm(as.formula(paste0("dens_avg ~", form1)), data = db, family = Gamma); summary(fit_mult_glm1)
fit_mult_glm2 <- glm(as.formula(paste0("dens_avg ~", form2)), data = db, family = Gamma); summary(fit_mult_glm2)
fit_mult_glm3 <- glm(as.formula(paste0("dens_avg ~", form3)), data = db, family = Gamma); summary(fit_mult_glm3)

# anova(fit_mult_glm2, fit_mult_glm1, test = "LRT") # Reduced model
# anova(fit_mult_glm3, fit_mult_glm2, test = "LRT") # Reduced model
# 
# hnp::hnp(fit_mult_glm1)
# hnp::hnp(fit_mult_glm2)
# hnp::hnp(fit_mult_glm3)

########################################## 
# GENERALIZED LINEAR MODEL - GAMMA - LOG #
##########################################

# UNIVARIATE #
fit_univ_glm_log1 <- glm(dens_avg ~ factor(braak06), data = db, family = Gamma(link=log)); summary(fit_univ_glm_log1)
fit_univ_glm_log2 <- glm(dens_avg ~ braak06 + I(braak06-3):braak06_grp, data = db, family = Gamma(link=log)); summary(fit_univ_glm_log2)
fit_univ_glm_log3 <- glm(dens_avg ~ braak06, data = db, family = Gamma(link=log)); summary(fit_univ_glm_log3)

# anova(fit_univ_glm_log2, fit_univ_glm_log1, test = "LRT") # Reduced model
# anova(fit_univ_glm_log3, fit_univ_glm_log2, test = "LRT") # Reduced model
# 
# hnp::hnp(fit_univ_glm_log1)
# hnp::hnp(fit_univ_glm_log2)
# hnp::hnp(fit_univ_glm_log3)

# MULTIVARIATE #
fit_mult_glm_log1 <- glm(as.formula(paste0("dens_avg ~", form1)), data = db, family = Gamma(link = log)); summary(fit_mult_glm_log1)
fit_mult_glm_log2 <- glm(as.formula(paste0("dens_avg ~", form2)), data = db, family = Gamma(link = log)); summary(fit_mult_glm_log2)
fit_mult_glm_log3 <- glm(as.formula(paste0("dens_avg ~", form3)), data = db, family = Gamma(link = log)); summary(fit_mult_glm_log3)

# anova(fit_mult_glm_log2, fit_mult_glm_log1, test = "LRT") # Reduced model
# anova(fit_mult_glm_log3, fit_mult_glm_log2, test = "LRT") # Reduced model
# 
# hnp::hnp(fit_mult_glm_log1)
# hnp::hnp(fit_mult_glm_log2)
# hnp::hnp(fit_mult_glm_log3)

##########
# TABLES #
##########

# AIC Models #

tab_aic <-
  data.frame(Model = c("Categorical", "Segmented Linear", "Linear"),
           AIC_lm = c(AIC(fit_mult_lm1),
                      AIC(fit_mult_lm2),
                      AIC(fit_mult_lm3)),
           AIC_glm = c(AIC(fit_mult_glm1),
                       AIC(fit_mult_glm2),
                       AIC(fit_mult_glm3)),
           AIC_glm_log = c(AIC(fit_mult_glm_log1),
                           AIC(fit_mult_glm_log2),
                           AIC(fit_mult_glm_log3))) %>% 
  gt() %>% 
  fmt_number(
    columns = vars(AIC_lm, AIC_glm, AIC_glm_log),
    decimals = 2) %>% 
  cols_label(Model = md("**Braak06 Effect**"),
             AIC_lm = md("**LM - log(outcome)**"),
             AIC_glm = md("**Gamma - inverse link**"),
             AIC_glm_log = md("**Gamma - log link**")) %>% 
  tab_spanner(label = md("**AIC**"),
              columns = 2:4) 

# Regression Results #

label <- list(age ~ 'Age (years)',
              pmi ~ 'PMI (minutes)',
              #area_cw ~ 'Area (10^6)',
              #apoe_allele_grp1 ~ 'APOE ALLELE',
              apoe_allele_grp2 ~ 'Apoe Allele',
              cognitive_status_num ~ 'Cognitive Status',
              sex ~ 'Sex',
              braak06 ~ "Braak6 Stage")

tab_lm <-
  tbl_regression(fit_mult_lm3,
                 label = label) %>% 
  bold_labels() %>% 
  add_global_p(keep = TRUE)

tab_glm <-
  tbl_regression(fit_mult_glm2,
                 label = label) %>% 
  bold_labels() %>% 
  add_global_p(keep = TRUE)

tab_glm_log <-
  tbl_regression(fit_mult_glm_log1,
                 exponentiate = TRUE,
                 label = list(age ~ 'Age (years)',
                              pmi ~ 'PMI (minutes)',
                              apoe_allele_grp2 ~ 'Apoe Allele',
                              cognitive_status_num ~ 'Cognitive Status',
                              sex ~ 'Sex',
                              braak06_fac ~ "Braak6 Stage")) %>% 
  bold_labels() %>% 
  add_global_p(keep = TRUE)

# AD-HOC COMPARISON #

comp <- glht(fit_mult_glm_log1, linfct = mcp(braak06_fac = "Tukey"))
comp_sum <- summary(comp, test = adjusted(type = "bonferroni"))

tab_comp <-
  data.frame(pval = comp_sum$test$pvalues) %>% 
  rownames_to_column(var = "Comparison") %>% 
  gt() %>% 
  fmt_number(
    columns = vars(pval),
    decimals = 2) %>% 
  cols_label(pval = md("**Adjusted p-value**"),
             Comparison = md("**Comparison**")) 
