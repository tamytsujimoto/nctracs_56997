require(tidyverse)
library(gtsummary)

# Loaging datasets

db = 
  read.csv("../../1 - From PI/APOE Dataset Formula Scrubbed 7-10-21.csv",
           na.strings = c("--", "n/a")) 

db2 <-
  db %>% 
  mutate(apoe_allele_grp1 = ifelse(apoe_allele %in% c("E23", "E33"), "E23/E33", apoe_allele),
         apoe_allele_grp2 = ifelse(apoe_allele %in% c("E23", "E33"), "E23/E33", "E34/E44"),
         area_cw = area_cw/1000000,
         cognitive_status_num = factor(cognitive_status_num, 
                                       labels = c("Cognitively Normal", "Mild Impairment", "Dementia")),
         across(ad_npath_score:arteriosc, ~factor(.x, labels = c("None", "Mild", "Moderate", "Severe")))) %>% 
  filter(!is.na(braak06))

saveRDS(db2, file = "brain.rds")






