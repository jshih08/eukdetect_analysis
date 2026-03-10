library(purrr)
library(tidyverse)

#function to format the eukdetect and metaphlan4 files for analysis
format_eukdetect_metaphlan <- function(eukdetect_file, metaphlan_file){
  euk_out <- read_delim(eukdetect_file)

    eukdetect_protist_presence <-
    euk_out %>%
      filter(Rank == "species") %>% #RPKSM values are only valid at species level
      filter(str_detect(Name, "Blastocystis") |  str_detect(Name, "Dient")) %>%
      select(-Lineage, -Rank, -TaxID) %>%
      pivot_longer(-Name, values_to = "rkpsm", names_to = "sample") %>%
      mutate(Name = str_c("Euk_", Name)) %>% 
      filter(str_detect(sample, "RPKSM")) %>%
      pivot_wider(values_from = "rkpsm", names_from = Name) %>% 
      mutate(sample = str_remove(sample, "_.*"))

  metaphlan_out <- read_delim(metaphlan_file, skip = 1) %>%  # Only keeping species level abundances
      filter(str_detect(clade_name, "s__")) %>%
      filter(!str_detect(clade_name, "t__"))

  #check that metaphaln abunadnces are sum to 100
  metaphlan_out %>% 
    summarise(across(where(is.numeric), sum, na.rm = TRUE))

  cleaned_dat <- 
    metaphlan_out %>%
    mutate(clade_name = str_replace(clade_name, ".*\\|s__", "")) %>%
    column_to_rownames("clade_name") %>% 
    t() %>%
    microbiome::transform("clr") %>% # CLR transform relative abundances
    as_tibble(rownames = "sample") %>%
    mutate(sample = str_replace(sample, "d_metagenome", "")) %>%
    left_join(eukdetect_protist_presence) %>% 
    pivot_longer(-c(sample, starts_with("Euk_")) , names_to = "taxa", values_to = "abundance")

  return(cleaned_dat)
}

blast_model_run <- function(cleaned_dat){
  blast_results_by_taxon <-  cleaned_dat %>%
    group_by(taxa) %>%
    nest() %>%
    mutate(
      mod = map(data, ~aov(abundance ~ Blastocystis, data = .x)),
      anova_summary = map(mod, ~summary(.x)[[1]]),
      
      # df = map_dbl(anova_summary, ~.x["Blastocystis", "Df"]),
      # sum_sq = map_dbl(anova_summary, ~.x["Blastocystis", "Sum Sq"]),
      # mean_sq = map_dbl(anova_summary, ~.x["Blastocystis", "Mean Sq"]),
      f_statistic = map_dbl(anova_summary, ~.x["Blastocystis", "F value"]),
      p.value = map_dbl(anova_summary, ~.x["Blastocystis", "Pr(>F)"]),
      
      coefficient = map_dbl(mod, ~coef(.x)["Blastocystis"]),
      conf_int = map(mod, ~confint(.x, "Blastocystis", level = 0.95)),
      cilower = map_dbl(conf_int, ~.x[1]),
      ciupper = map_dbl(conf_int, ~.x[2])
    ) %>%
    ungroup() %>%
    mutate(
      p.adj.hochberg = p.adjust(p.value, method = "hochberg")
    ) %>%
    select(taxa, coefficient, p.adj.hochberg, cilower, ciupper,f_statistic)
  return(blast_results_by_taxon)
}

dient_model_run <- function(cleaned_dat){
  dient_results_by_taxon <-  cleaned_dat %>%
    group_by(taxa) %>%
    nest() %>%
    mutate(
      mod = map(data, ~aov(abundance ~ Dientamoeba, data = .x)),
      anova_summary = map(mod, ~summary(.x)[[1]]),
      
      # df = map_dbl(anova_summary, ~.x["Dientamoeba", "Df"]),
      # sum_sq = map_dbl(anova_summary, ~.x["Dientamoeba", "Sum Sq"]),
      # mean_sq = map_dbl(anova_summary, ~.x["Dientamoeba", "Mean Sq"]),
      f_statistic = map_dbl(anova_summary, ~.x["Dientamoeba", "F value"]),
      p.value = map_dbl(anova_summary, ~.x["Dientamoeba", "Pr(>F)"]),
      
      coefficient = map_dbl(mod, ~coef(.x)["Dientamoeba"]),
      conf_int = map(mod, ~confint(.x, "Dientamoeba", level = 0.95)),
      cilower = map_dbl(conf_int, ~.x[1]),
      ciupper = map_dbl(conf_int, ~.x[2])
    ) %>%
    ungroup() %>%
    mutate(
      p.adj.hochberg = p.adjust(p.value, method = "hochberg")
    ) %>%
    select(taxa, coefficient, p.adj.hochberg, cilower, ciupper,f_statistic)
  return(dient_results_by_taxon)
}



itxn_model_run <- function(cleaned_dat){
  inter_results_by_taxon <- cleaned_dat %>%
    group_by(taxa) %>%
    nest() %>%
    mutate(
      mod = map(data, ~aov(abundance ~ Blastocystis + Dientamoeba + Blastocystis * Dientamoeba, data = .x)),
      anova_summary = map(mod, ~summary(.x)[[1]]),
      # df = map_dbl(anova_summary, ~.x[1, "Df"]),
      # sum_sq = map_dbl(anova_summary, ~.x[1, "Sum Sq"]),
      # mean_sq = map_dbl(anova_summary, ~.x[1, "Mean Sq"]),
      # f_statistic = map_dbl(anova_summary, ~.x[1, "F value"]),
      p.value_blast = map_dbl(anova_summary, ~.x[1, "Pr(>F)"]),
      p.value_dient = map_dbl(anova_summary, ~.x[2, "Pr(>F)"]),
      p.value_itxn = map_dbl(anova_summary, ~.x[3, "Pr(>F)"]),

      coefficient_blast = map_dbl(mod, ~coef(.x)["Blastocystis"]),
      coefficient_dient = map_dbl(mod, ~coef(.x)["Dientamoeba"]),
      coefficient_itxn = map_dbl(mod, ~coef(.x)["Blastocystis:Dientamoeba"]),

      conf_int_blast = map(mod, ~confint(.x, "Blastocystis", level = 0.95)),
      cilower_blast = map_dbl(conf_int_blast, ~.x[1]),
      ciupper_blast = map_dbl(conf_int_blast, ~.x[2]),

      conf_int_dient = map(mod, ~confint(.x, "Dientamoeba", level = 0.95)),
      cilower_dient = map_dbl(conf_int_dient, ~.x[1]),
      ciupper_dient = map_dbl(conf_int_dient, ~.x[2]),

      conf_int_itxn = map(mod, ~confint(.x, "Blastocystis:Dientamoeba", level = 0.95)),
      cilower_itxn = map_dbl(conf_int_itxn, ~.x[1]),
      ciupper_itxn = map_dbl(conf_int_itxn, ~.x[2]),

    ) %>%
    ungroup() %>% 
    mutate(
      p.adj.hochberg_blast = p.adjust(p.value_blast, method = "hochberg"),
      p.adj.hochberg_dient = p.adjust(p.value_dient, method = "hochberg"),
      p.adj.hochberg_itxn = p.adjust(p.value_itxn, method = "hochberg")
    ) %>%
    select(taxa, starts_with("coefficient"), starts_with("p.adj.hochberg"), starts_with("cilower"), starts_with("ciupper"))

  return(inter_results_by_taxon)
}