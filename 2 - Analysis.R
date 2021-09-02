library(tidyverse)
library(gtsummary)

# LOADING DATASET

db <- 
  readRDS("brain.rds")

##########################
# TABLE 1 - DEMOGRAPHICS #
##########################

tab1_dem <-
  db %>%
  select(-c(patient_id, 
            other_notes, 
            outlier_criteria, 
            ca_density, 
            e4_allele_dose, 
            cognitive_status, 
            apoe_allele_grp1,
            apoe_allele_grp2)) %>%
  tbl_summary(label = list(area_cw ~ "Area (10^6)",
                           ca_number ~ "CA Number",
                           #ca_density ~ "DENSITY",
                           cw_adj ~ "CW Adjusted Density (mm2)",
                           sm_adj ~ "SM Adjusted Density (mm2)",
                           vb_adj ~ "VB Adjusted Density (mm2)",
                           dens_avg ~ "All Density Avg (mm2)",
                           apoe_allele ~ "Apoe Allele",
                           #e4_allele_dose ~ "E4 Allele Dose",
                           age ~ "Age (years)",
                           sex ~ "Sex",
                           ethnicity ~ "Ethnicity",
                           #cognitive_status ~ "COGNITIVE STATUS",
                           cognitive_status_num ~ "Cognitive Status",
                           braak06 ~ "Braak6 Score",
                           pmi ~ "PMI (minutes)",
                           ad_npath_score ~ "AD Neuropath score",
                           amy_angi ~ "Amyloid Angiopathy",
                           athero ~ "Atherosclerosis",
                           arteriosc ~ "Arteriosclerosis",
                           hipp_scle ~ "Hippocampal Sclerosis (0 = not present, 1 = present)",
                           mck_dlb ~ "McKeith DLB",
                           infarction ~ "Infarction ( 0 = not  present, 1 = present)",
                           hemorr ~ "Hemorrhage ( 0 = not present, 1 = present )",
                           artag ~ "ARTAG (0 = not present , 1 = Mild ARTAG )",
                           late ~ "LATE: Limbic-PAR TDP-43  (0 =  not present, 1 = present )",
                           a_score ~ "A Score",
                           a_thal_score ~ "A Score (Thal)",
                           b_score ~ "B Score (Braak)",
                           c_score ~ "C Score (CERAD)")
                           ) %>% 
  bold_labels()

tab_braak06 <-
  db %>% 
  select(braak06,
         dens_avg) %>%
  tbl_summary(by = braak06,
              label = list(dens_avg ~ "All Density Avg (mm2)",
                           braak06 ~ "BRAAK 06")) %>% 
  bold_labels() %>% 
  add_p()

##########
# GRAPHS #
##########

# HISTOGRAM #

p_hist <- 
  db %>% 
  ggplot(aes(x = dens_avg)) +
  geom_histogram(color="black", fill = "gray75") +
  theme_bw() +
  scale_x_continuous(name ="All Density Avg (mm2)") + 
  scale_y_continuous(name ="Count") 

p_hist_log <- 
  db %>% 
  ggplot(aes(x = log(dens_avg))) +
  geom_histogram(color="black", fill = "gray75") +
  theme_bw() +
  scale_x_continuous(name ="Log(All Density Avg (mm2))") + 
  scale_y_continuous(name ="Count") 

ggpubr::ggarrange(p_hist,p_hist_log, ncol=2, nrow=1)
ggsave('Output/hist.png', height=6, width=10)

# SCATTER PLOT #

p_scat <- 
  db %>% 
  ggplot(aes(x = braak06, y = dens_avg)) +
  geom_point() +
  geom_smooth() +
  theme_bw() +
  scale_x_continuous(name ="Braak6 Score") + 
  scale_y_continuous(name ="All Density Avg (mm2)") 

p_scat_log <- 
  db %>% 
  ggplot(aes(x = braak06, y = log(dens_avg))) +
  geom_point() +
  geom_smooth() +
  theme_bw() +
  scale_x_continuous(name ="Braak6 Score") + 
  scale_y_continuous(name ="Log(All Density Avg (mm2))") 

ggpubr::ggarrange(p_scat,p_scat_log, ncol=2, nrow=1)
ggsave('Output/scatter.png', height=6, width=10)

# BOX PLOT #

p_box <- 
  db %>% 
  filter(!is.na(braak06)) %>% 
  ggplot(aes(x = factor(braak06), y = dens_avg)) +
  geom_boxplot() +
  theme_bw() +
  scale_x_discrete(name ="Braak6 Score") + 
  scale_y_continuous(name ="All Density Avg (mm2)") 

p_box_log <- 
  db %>% 
  filter(!is.na(braak06)) %>% 
  ggplot(aes(x = factor(braak06), y = log(dens_avg))) +
  geom_boxplot() +
  geom_smooth() +
  theme_bw() +
  scale_x_discrete(name ="Braak6 Score") + 
  scale_y_continuous(name ="Log(All Density Avg (mm2))") 

ggpubr::ggarrange(p_box,p_box_log, ncol=2, nrow=1)
ggsave('Output/box.png', height=6, width=10)





  











