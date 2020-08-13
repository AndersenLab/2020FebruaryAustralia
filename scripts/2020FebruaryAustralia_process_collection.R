#devtools::install_github("AndersenLab/easyfulcrum")
library(easyfulcrum)
library(tidyverse)

# Set working directory
setwd(glue::glue("{dirname(rstudioapi::getActiveDocumentContext()$path)}/.."))

# process Fulcrum data and check
raw_data <- readFulcrum('data/raw/fulcrum')

proc_data <- procFulcrum(raw_data)

checkParameters(proc_data)

checkProc(proc_data)

join_data <- joinFulcrum(proc_data)

checkJoin(join_data)

anno_data <- annotateFulcrum(join_data)

# load genotyping data. 
geno_data <- readGenotypes(gsKey = "1_6u4sk_Zj-Hm5d_058Lg8WYWLe7BZHGTWxXcH6EsDUI")

# test procGenotypes function
proc_geno_data <- procGenotypes(geno_data = geno_data, fulc_data = anno_data)

checkGenotypes(proc_geno_data)

# test joinGenoFulc function to join genotype data to fulcrum data
joingeno_data <- joinGenoFulc(geno = proc_geno_data, fulc = anno_data)

# test the procPhotos function, output is final dataframe
final_data1 <- procPhotos(dir = "data/raw/fulcrum/photos", data = joingeno_data, max_dim = 500, overwrite = T)

# make initial species sheet
species_sheet1 <- makeSpSheet(joingeno_data)

# correct email addresses to Names
species_sheet2 <- species_sheet1 %>%
  dplyr::mutate(sampled_by = ifelse(sampled_by == "dec@u.northwestern.edu", "D. Cook", "wtf?"),
                isolated_by = ifelse(isolated_by == "claire.buchanan@northwestern.edu", "C. Buchanan")) %>%
  dplyr::select(-flag_sampled_by_is_email_address, -flag_isolated_by_is_email_address, -flag_species_not_in_target_species)

# save species sheet
current_date <- stringr::str_extract(Sys.time(), pattern = regex("^.{10}")) %>%
  stringr::str_replace_all(., pattern = "-", replacement = "")

# export species sheet
rio::export(species_sheet2, glue::glue("reports/{current_date}_2020JanuaryHawaii_species_sheet.csv"))

