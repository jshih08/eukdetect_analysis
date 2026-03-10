library(purrr)
library(tidyverse)

source("eukdetect_analysis_functions.R")

cohorts <- c("predict1","predict2","predict3_uk",
             "predict3_us1","predict3us22a","hadza")


# formats and prepares for model
# walk(cohorts, ~
#   format_eukdetect_metaphlan(
#     str_c("data/", .x, "_eukdetect_normalized_rpksms.txt"),
#     str_c("data/", .x, "_metaphlan4_merged_abundances.txt")
#   ) %>%
#     saveRDS(str_c("data/cleaned_", .x, "_dat_clr.rds"))
# )

# # runs blast model and saves model outputs
# walk(cohorts, ~ {

#   input <- str_c("data/cleaned_", .x ,"_dat_clr.rds")
#   out <- str_c("data/", .x, "_blast_interactions.csv")

#   dat <- readRDS(input) %>%
#     mutate(Blastocystis = rowSums(across(starts_with("Euk_Blast"))))

#   write_csv(blast_model_run(dat), out)
  
# })

# runs dient model and saves model outputs
# walk(cohorts, ~ {

#   input <- str_c("data/cleaned_", .x ,"_dat_clr.rds")
#   out <- str_c("data/", .x, "_dient_interactions.csv")

#   dat <- readRDS(input) %>%
#     mutate(Dientamoeba = rowSums(across(starts_with("Euk_Dient"))))

#   write_csv(dient_model_run(dat), out)
  
# })

walk(cohorts, ~ {

  input <- str_c("data/cleaned_", .x ,"_dat_clr.rds")
  out <- str_c("data/", .x, "_itxn_interactions.csv")

  dat <- readRDS(input) %>%
    mutate(Blastocystis = rowSums(across(starts_with("Euk_Blast")))) %>%
    mutate(Dientamoeba = rowSums(across(starts_with("Euk_Dient"))))

  write_csv(itxn_model_run(dat), out)
  
})



