#devtools::install_github("AndersenLab/easyfulcrum")
library(easyfulcrum)
library(tidyverse)
library(rebus)

# Set working directory
setwd(glue::glue("{dirname(rstudioapi::getActiveDocumentContext()$path)}/.."))

# process Fulcrum data and check
raw_data <- readFulcrum('data/raw/fulcrum')

proc_data <- procFulcrum(raw_data)

# There are rows with high temps that are converted from F to C
par_issues <- checkParameters(proc_data, return = TRUE)

# # Run fixParameters
par_issues_to_fix <- par_issues$ambient_temperature_run %>%
  dplyr::filter(flag_ambient_temperature_run == TRUE) %>%
  dplyr::pull(fulcrum_id)

proc_data2 <- fixParameters(proc_data, ambient_temperature_run_ids = par_issues_to_fix)
checkParameters(proc_data2)

# need to run fixParameters function but no advice on how to do that in documentation.
# Need to add to easyfulcrum in printed message (include return = TRUE, if necessary) and also in checkParameters documentation
issues_fulcrum <- checkProc(proc_data2, return = TRUE)

# Running joinFulcrum function to assess how to deal with flags generated in procFulcrum
join_data <- joinFulcrum(proc_data2)

########################################################
### correcting duplicated c labels in join_data. 
########################################################
join_data2 <- join_data %>%
  dplyr::mutate(remove_dup_c_labels = ifelse((flag_duplicated_c_label_field_sampling == TRUE & 
                                                   flag_missing_isolation_record == TRUE), 1, 0)) %>%
  dplyr::filter(remove_dup_c_labels == 0) %>%
  dplyr::select(-remove_dup_c_labels) %>%
  dplyr::mutate(flag_duplicated_c_label_field_sampling = FALSE)
# All duplicated C-labels have one instance with a missing isolation and another without a missing isolation. In each case
# the C-label with the missing isolation record is removed.

########################################################
#### correcting falg_duplicated_isolation_for_c_label
########################################################
join_data3 <- join_data2 %>%
  dplyr::group_by(c_label) %>%
  dplyr::mutate(best_iso_record = min(isolation_datetime_UTC)) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(remove_dup_isolations = ifelse((flag_duplicated_isolation_for_c_label == TRUE & as.character(isolation_datetime_UTC) != best_iso_record), 1, 0)) %>%
  dplyr::select(-best_iso_record, -remove_dup_isolations)
# All the duplications are instances where a C-label was isolated from twice.
# All instances were for C-label with  'Tracks only' for the first and  second isolation.
# For each instance the first isolation record will be retained and the later isolation removed.

################################################################
#### correcting flag_missing_s_label_isolation_s_labeled_plates
################################################################
# correcting flag_missing_s_label_isolation_s_labeled_plates.
# Two C-labels have isolation records but not s-labels attached to them. It's not
# clear how to resolve this issue yet. Possibly there are S-labels in the genotyping sheet that
# are not present in Fulcrum and those could help.

checkJoin(join_data3)

################################################################
#### correcting flag_missing_isolation_records
################################################################
# These records might just be lost C-labels
# leaving them for now, but they may have to be filtered out

# Annotate collection with location and trail data
anno_data <- annotateFulcrum(join_data3)

################################################################
### Joining Matt Rockmann genotype dataset
################################################################
rock_dat <- readr::read_csv("data/raw/genotypes/Rockman_OzFunnelStrains.csv") %>%
  dplyr::mutate(s_label = ifelse(s_label == "5537.2b", "5437.2b", s_label)) %>% # fix mislabelled S-label. CHECK WITH ROCKMAN THAT HE IS NOT USING THIS LABEL.
  dplyr::mutate(s_label_prefix = "S-",
                s_label = paste0(s_label_prefix, s_label)) %>% # add appropriate S-label prefix
  dplyr::select(-s_label_prefix) %>%
  dplyr::mutate(c_label = stringr::str_replace_all(s_label, pattern = "S-", replacement = "C-"),
                c_label = stringr::str_replace_all(c_label, pattern = "\\..*", replacement = "")) %>% # add C-label
  dplyr::mutate(s_label = ifelse(c_label == "C-5250", stringr::str_replace(s_label, pattern = "S-", replacement = ""), s_label)) # fix special case EDITED FULCRUM SO RELOAD THESE THIS LINE ISN"T NEEDED

###############################################################
### DO NOT REMOVE Fulcrum Nemtatode Isolation edits 20200814
###############################################################
# The S-labels S-labels associated with C-5250 were edited  from 5250.1 - 5250.6 to S-5250.1 - S-5250.6 by TAC 20200814
# S-5421.4a and S-5421.4b S-labels were added to Fulcrum manually by TAC 20200814
# The S-label S-5437.2a and S-5437.2b were added to Fulcrum manually by TAC 20200814. Also, S-5537.2b was renamed S-5437.2b above. CHECK WITH M. ROCKMAN
# The S-labels S-5456.1a and S-5456.1b were added to Fulcrum manually by TAC 20200814.
# The S-label 5460.2 was edited to S-5460.2 in Fulcrum by TAC 20200814.
# The S-label S-5508.1a was added to Fulcrum by TAC 20200814.


rock_dat_in_fulcrum <- rock_dat %>%
  dplyr::filter(s_label %in% anno_data$s_label)

rock_dat_not_in_fulcrum <- rock_dat %>%
  dplyr::filter(!(s_label %in% anno_data$s_label))



anno_data2 <- anno_data %>%
  dplyr::left_join(rock_dat, )




# load genotyping data. 
geno_data <- readGenotypes(gsKey = "1XqR-HaQOEcCc6f-SDw-SX3_8V9aNqBpW61HdgcikEws")

# use procGenotypes function to include flags in genotyping data
proc_geno_data <- procGenotypes(geno_data = geno_data, fulc_data = anno_data)

# use checkGenotypes fucntion to return records with certain flags
issues_geno <- checkGenotypes(proc_geno_data, return = TRUE)

######################################################
### Correcting flag_s_label_not_in_fulcrum
######################################################


######################################################
### Correcting flag_unusual_target_species_name
######################################################

# Use joinGenoFulc function to join genotype data to fulcrum data
joingeno_data <- joinGenoFulc(geno = proc_geno_data, fulc = anno_data)

# test the procPhotos function, output is final dataframe
final_data1 <- procPhotos(dir = "data/raw/fulcrum/photos", data = joingeno_data, max_dim = 500, overwrite = T)

# make initial species sheet
species_sheet1 <- makeSpSheet(joingeno_data)

# correct email addresses to Names
species_sheet2 <- species_sheet1 %>%
  dplyr::mutate(sampled_by = ifelse(sampled_by == "erik.andersen@northwestern.edu", "E. Andersen", "wtf?"),
                isolated_by = ifelse(isolated_by == "erik.andersen@northwestern.edu", "E. Andersen", "wtf?")) %>%
  dplyr::select(-flag_sampled_by_is_email_address, -flag_isolated_by_is_email_address, -flag_species_not_in_target_species)

# save species sheet
current_date <- stringr::str_extract(Sys.time(), pattern = regex("^.{10}")) %>%
  stringr::str_replace_all(., pattern = "-", replacement = "")

# export species sheet
rio::export(species_sheet2, glue::glue("reports/{current_date}_2020FebruaryUCLA_species_sheet.csv"))





1XqR-HaQOEcCc6f-SDw-SX3_8V9aNqBpW61HdgcikEws