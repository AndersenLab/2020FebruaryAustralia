library(easyfulcrum)
library(tidyverse)
library(rebus)

# Set working directory
setwd(glue::glue("{dirname(rstudioapi::getActiveDocumentContext()$path)}/.."))

#=====================================================#
# PART 1: process Fulcrum data
#=====================================================#
raw_data <- readFulcrum('data/raw/fulcrum')

proc_data <- procFulcrum(raw_data)

# There are rows with high temps that are converted from F to C
par_issues <- checkParameters(proc_data, return = TRUE)

# Run fixParameters to set run to NA
par_issues_to_fix <- par_issues$ambient_temperature_run %>%
  dplyr::filter(flag_ambient_temperature_run == TRUE) %>%
  dplyr::pull(fulcrum_id)

proc_data2 <- fixParameters(proc_data, ambient_temperature_run_ids = par_issues_to_fix)

# Recheck proc_data2
checkParameters(proc_data2)

# Check rest of data
issues_fulcrum <- checkProc(proc_data2, return = TRUE)

# Running joinFulcrum function to assess how to deal with flags generated in procFulcrum
join_data <- joinFulcrum(proc_data2)

#=====================================================#
#  correcting duplicated c labels in join_data. 
#=====================================================#
join_data2 <- join_data %>%
  dplyr::mutate(remove_dup_c_labels = ifelse((flag_duplicated_c_label_field_sampling == TRUE & 
                                                   flag_missing_isolation_record == TRUE), 1, 0)) %>%
  dplyr::filter(remove_dup_c_labels == 0) %>%
  dplyr::select(-remove_dup_c_labels) %>%
  dplyr::mutate(flag_duplicated_c_label_field_sampling = FALSE)
# All duplicated C-labels have one instance with a missing isolation and another without a missing isolation. In each case
# the C-label with the missing isolation record is removed.

#=====================================================#
#  correcting flag_duplicated_isolation_for_c_label
#=====================================================#
join_data3 <- join_data2 %>%
  dplyr::group_by(c_label) %>%
  dplyr::mutate(best_iso_record = min(isolation_datetime_UTC)) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(remove_dup_isolations = ifelse((flag_duplicated_isolation_for_c_label == TRUE & as.character(isolation_datetime_UTC) != best_iso_record), 1, 0)) %>%
  dplyr::select(-best_iso_record, -remove_dup_isolations)
# All the duplications are instances where a C-label was isolated from twice.
# All instances were for C-label with  'Tracks only' for the first and  second isolation.
# For each instance the first isolation record will be retained and the later isolation removed.

#=================================================================#
#  correcting flag_missing_s_label_isolation_s_labeled_plates
#=================================================================#
# correcting flag_missing_s_label_isolation_s_labeled_plates.
# Two C-labels have isolation records but not s-labels attached to them. It's not
# clear how to resolve this issue yet. Possibly there are S-labels in the genotyping sheet that
# are not present in Fulcrum and those could help.

# Check the joined data
checkJoin(join_data3)

#=====================================================#
#  correcting flag_missing_isolation_records
#=====================================================#
# These records might just be lost C-labels
# leaving them for now, but they may have to be filtered out if the isolation records cannot be found

# Annotate collection with location and trail data
anno_data <- annotateFulcrum(join_data3)

#=====================================================#
# PART 2: load and check genotyping data
#=====================================================#
geno_data <- readGenotypes(gsKey = "1XqR-HaQOEcCc6f-SDw-SX3_8V9aNqBpW61HdgcikEws")

# Add Rockman genotying data
rock_dat <- readr::read_csv("data/processed/genotypes/20200828_Rockman_OzFunnel_processed.csv")

rock_collection_info <- rock_dat %>%
  dplyr::select(strain_name, strain_lost:c_label)

rock_genotypes <- rock_dat %>%
  dplyr::select(project_id:pair_name)

# Check for rockman data not in fulcrum
rock_dat_not_in_fulcrum <- rock_dat %>%
  dplyr::filter(!(s_label %in% anno_data$s_label))

# append Rockman s_labels and genotypes to our sheet
geno_data2 <- geno_data %>%
  dplyr::bind_rows(rock_genotypes)

# use procGenotypes function to include flags in genotyping data
proc_geno_data <- procGenotypes(geno_data = geno_data2, fulc_data = anno_data)

# use checkGenotypes fucntion to return records with certain flags
issues_geno <- checkGenotypes(geno_data = proc_geno_data, fulc_data = anno_data, return = TRUE)

#=====================================================#
#  Correcting flag_unusual_target_species_name
#=====================================================#
proc_geno_data2 <- proc_geno_data %>%
  dplyr::mutate(species_id = ifelse(species_id == "C. elegans", "Caenorhabditis elegans",
                                    ifelse(species_id == "C. briggsae", "Caenorhabditis briggsae", species_id)))

proc_geno_data3 <- procGenotypes(geno_data = proc_geno_data2, fulc_data = anno_data)

issues_geno2 <- checkGenotypes(geno_data = proc_geno_data3, fulc_data = anno_data, return = T)

#=============================================================#
# Corecting 4 s labels not found in fulcrum
#=============================================================#
# the s-labels "S-13381" "S-13570" "S-13585" "S-13657" are not in fulcrum
# there is no clear reason why these are not found in fulcrum
# they wil be filtered from dataset after geno and fulcrum joined below.

#=============================================================#
# Corecting 6 s labels found in fulcrum but not in genotyping
#=============================================================#
# The s-labels "S-13658" "S-13659" "S-13623" "S-13660" "S-13645" "S-13646"
# all belong to C-5256. Since this collection effectively has no isolation data
# it will be removed from the dataset after geno and fulcrum joined below.

#=============================================================#
# Checking 231 s labels missing its2 genotypes
#=============================================================#
# These will be left as is b/c they are Rockman genotypes that did 
# not get an its2 pcr or blast.

#=============================================================#
#   DO NOT REMOVE Fulcrum Nemtatode Isolation edits 20200814
#=============================================================#
# The S-labels  associated with C-5250 were edited  from 5250.1 - 5250.6 to S-5250.1 - S-5250.6 by TAC 20200814
# S-5421.4a and S-5421.4b S-labels were added to Fulcrum manually by TAC 20200814
# The S-label S-5437.2a and S-5437.2b were added to Fulcrum manually by TAC 20200814. Also, S-5537.2b was renamed S-5437.2b above.
# The S-labels S-5456.1a and S-5456.1b were added to Fulcrum manually by TAC 20200814.
# The S-label 5460.2 was edited to S-5460.2 in Fulcrum by TAC 20200814.
# The S-label S-5508.1a was added to Fulcrum by TAC 20200814.
# The S-labels S-5415.5 and S-5415.6 were added to Fulcrum by TAC 20200824.

# Use joinGenoFulc function to join genotype data to fulcrum data
joingeno_data <- joinGenoFulc(geno = proc_geno_data3, fulc = anno_data)

# Make changes to remove s_labels flagged in genotyping
joingeno_data2 <- joingeno_data %>%
  dplyr::filter(!(s_label %in% issues_geno2$s_label_not_in_fulcrum$s_label)) %>%
  dplyr::filter(c_label != "C-5256")

# test the procPhotos function, output is final dataframe. 
final_data1 <- procPhotos(dir = "data/raw/fulcrum/photos", data = joingeno_data2, max_dim = 500, overwrite = T)

# remove the photos that are not C. elegans for now.
img_rm <- tibble(fs::dir_ls("data/processed/fulcrum/processed_photos", recurse = F, type = "file")) %>%
  dplyr::select(img_path = `fs::dir_ls(...)`) %>%
  dplyr::mutate(strain_name = stringr::str_extract(string = img_path, pattern = "ECA" %R% DGT %R% DGT %R% DGT %R% optional(DGT))) %>%
  dplyr::filter(!(strain_name %in% (joingeno_data2 %>% dplyr::filter(species_id == "Caenorhabditis elegans" & !is.na(strain_name)) %>% dplyr::pull(strain_name))))

fs::file_delete(img_rm$img_path)

img_rm_thumbs <- tibble(fs::dir_ls("data/processed/fulcrum/processed_photos/thumbnails", recurse = F, type = "file")) %>%
  dplyr::select(img_path = `fs::dir_ls(...)`) %>%
  dplyr::mutate(strain_name = stringr::str_extract(string = img_path, pattern = "ECA" %R% DGT %R% DGT %R% DGT %R% optional(DGT))) %>%
  dplyr::filter(!(strain_name %in% (joingeno_data2 %>% dplyr::filter(species_id == "Caenorhabditis elegans" & !is.na(strain_name)) %>% dplyr::pull(strain_name))))

fs::file_delete(img_rm_thumbs$img_path)

# make initial species sheet
species_sheet1 <- makeSpSheet(joingeno_data2)

# correct email addresses to Names
species_sheet2 <- species_sheet1 %>%
  dplyr::mutate(sampled_by = ifelse(sampled_by == "erik.andersen@northwestern.edu", "E. Andersen",
                                    ifelse(sampled_by == "lewis.stevens07@gmail.com", "L. Stevens",
                                           ifelse(sampled_by == "mrockman@nyu.edu", "M. Rockman", "wtf?"))),
                isolated_by = ifelse(isolated_by == "claire.buchanan@northwestern.edu", "C. Buchanan",
                                     ifelse(isolated_by == "nicoleroberto2018@u.northwestern.edu", "N. Roberto",
                                            ifelse(isolated_by == "robyn.tanny@northwestern.edu", "R. Tanny",
                                                   ifelse(isolated_by == "lewis.stevens07@gmail.com", "L. Stevens",
                                                          ifelse(isolated_by == "lorainastinson2024@u.northwestern.edu", "L. Stinson", "wtf?")))))) %>%
  dplyr::select(-flag_sampled_by_is_email_address, -flag_isolated_by_is_email_address, -flag_species_not_in_target_species,
                -flag_unusual_substrate_class, -flag_unusual_landscape_class) %>%
  dplyr::filter(species %in% c("Caenorhabditis elegans",  "Caenorhabditis briggsae", "Caenorhabditis tropicalis")) %>%
  dplyr::left_join(select(rock_collection_info, strain = strain_name, spid_m = species_id_method, loc_d = location)) %>%
  dplyr::mutate(species_id_method = case_when(is.na(spid_m) ~ species_id_method,
                                              !is.na(spid_m) ~ spid_m)) %>% # add Rockman species_id_method
  dplyr::mutate(locality_description = case_when(locality_description == "NA" ~ NA_character_,
                                                 locality_description != "NA" ~ locality_description)) %>% #remove "NA" from locality descriptions and add rockman's
  dplyr::mutate(locality_description = loc_d) %>% #Add rockman's locality descriptions NEED TO ADD TO BIG DATAFRAME TOO IN FUTURE
  dplyr::arrange(strain) %>%
  dplyr::select(-spid_m, -loc_d)

ss_cb <- species_sheet2 %>%
  dplyr::filter(species == "Caenorhabditis briggsae")

ss_ce <- species_sheet2 %>%
  dplyr::filter(species == "Caenorhabditis elegans")

ss_ct <- species_sheet2 %>%
  dplyr::filter(species == "Caenorhabditis tropicalis")

# save species sheet
current_date <- stringr::str_extract(Sys.time(), pattern = regex("^.{10}")) %>%
  stringr::str_replace_all(., pattern = "-", replacement = "")

# export species sheet
rio::export(species_sheet2, glue::glue("reports/{current_date}_{unique(anno_data$project)}_species_sheet.csv"))
rio::export(ss_ce, glue::glue("reports/{current_date}_{unique(anno_data$project)}_species_sheet_elegans.csv"))
rio::export(ss_cb, glue::glue("reports/{current_date}_{unique(anno_data$project)}_species_sheet_briggsae.csv"))

# export final data
saveRDS(final_data1, file = glue::glue("data/processed/fulcrum/{current_date}_{unique(anno_data$project)}.rds"))
