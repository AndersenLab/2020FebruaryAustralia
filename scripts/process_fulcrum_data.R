#!/usr/bin/env Rscript
#Load necessary packages
library(tidyverse)
library(lubridate)
library(geosphere)
library(googlesheets)
library(geonames)
library(rebus)
library(sp)

# Set working directory
setwd(glue::glue("{dirname(rstudioapi::getActiveDocumentContext()$path)}/.."))

###################################################################
### Define functions                                            ###
###################################################################

filter_box <- function(longitude, latitude, coords) {
  between(longitude, coords[1], coords[3]) &
    between(latitude, coords[2], coords[4]) &
    !is.na(longitude)
}

FtoC <- function(F) {
  (F - 32)*(5/9)
}

 ###################################################################
### 1: Read field collection data (C-labels)                    ###
###################################################################

collection <- readr::read_csv("data/fulcrum/nematode_field_sampling.csv") %>%
  dplyr::mutate(c_label = stringr::str_to_upper(c_label)) %>%
  # name created_by to specify who picked up the sample
  dplyr::rename(collection_by = created_by) %>%
  dplyr::select(-updated_at,
                -system_created_at,
                -system_updated_at,
                -date) %>%
  # choose one sample photo only. This takes the first sample photo and warns if additional photos are discarded
  tidyr::separate(col = sample_photo, into = "sample_photo", sep = ",", extra = "warn") %>%
  # this is UTC time (very important if you want to convert to HST time)
  dplyr::mutate(collection_datetime_UTC = lubridate::ymd_hms(created_at, tz = "UTC")) %>% 
  # again this is UTC date (very important if you want to convert to HST date)
  dplyr::mutate(collection_date_UTC = lubridate::date(created_at)) %>% 
  dplyr::select(-created_at) %>%
  # Fix Fahrenheit observations
  dplyr::mutate(substrate_temperature = ifelse(substrate_temperature > 35,
                                               FtoC(substrate_temperature),
                                               substrate_temperature)) %>% 
  # Fix ambient temp F to C
  dplyr::mutate(ambient_temperature = ifelse(ambient_temperature_c > 50,
                                             FtoC(ambient_temperature_c),
                                             ambient_temperature_c)) %>%
  # force ambient temp to numeric
  dplyr::mutate(ambient_temperature = as.numeric(ambient_temperature)) %>%
  # drop ambient temp c
  dplyr::select(-ambient_temperature_c) %>%
  # add flags for runs of temperature data
  dplyr::arrange(collection_datetime_UTC) %>%
  dplyr::mutate(flag_ambient_temperature_run = (ambient_humidity == dplyr::lag(ambient_humidity)) &
                  (ambient_temperature == dplyr::lag(ambient_temperature))
                & (gridsect == "no"))

####################################################################
###      (OPTIONAL) CORRECT DUPLICATE COLLECTIONS (OPTIONAL)     ###
####################################################################

# In 3 instances there were the same c_label was used for differnt collections.
# selecting to retain the collection record based on isolation data and environmental data.
duplicated_collections_to_remove <- c("2e161ffe-01ee-45b3-8d28-28db8e070eb7",
                                      "9bbdd0cf-6eee-4b9f-a127-11e6039bae00",
                                      "1a822595-a9a4-4319-b57a-7045a3c87efc")

collection <- collection %>%
  dplyr::filter(!fulcrum_id %in% duplicated_collections_to_remove)

###################################################################
### 2: Read isolation data (S-labels)                           ###
###################################################################
# Read in S-labels
isolation <- readr::read_csv("data/fulcrum/nematode_isolation.csv") %>%
  dplyr::select(c_label_id = c_label,
                isolation_id = fulcrum_id,
                isolation_datetime_UTC = system_created_at,
                isolation_by = created_by,
                worms_on_sample,
                approximate_number_of_worms,
                isolation_date_UTC = date,
                isolation_local_time = time,
                isolation_latitude = latitude,
                isolation_longitude = longitude)

#############################################################################
### 3: Use exiftool to extract lat, long elevation. ONLY NEED TO RUN ONCE ###
#############################################################################

# # Read in data from photos. Need to install using ‘brew install exiftool’ in terminal.
# comm <- paste0("exiftool -coordFormat '%+.6f' -csv -ext jpg ",
#                 getwd(),
#                 "/photos/*")
# 
# # Exif Data
#  exif <- readr::read_csv(pipe(comm)) %>%
#    dplyr::mutate(SourceFile = stringr::str_replace(basename(SourceFile), ".jpg", "")) %>%
#    dplyr::select(sample_photo = SourceFile,
#                  altitude = GPSAltitude,
#                  latitude = GPSLatitude,
#                  longitude = GPSLongitude,
#                  ExposureTime,
#                  Artist,
#                  Aperture,
#                  BrightnessValue,
#                  FOV) %>%
#    dplyr::mutate(altitude =  as.numeric(stringr::str_replace(altitude, " m", ""))) %>%
#    dplyr::mutate(FOV =  as.numeric(stringr::str_replace(FOV, " deg", ""))) %>%
#    dplyr::group_by(sample_photo) %>%
#    # Only retain data from one sample photo.
#    dplyr::distinct(.keep_all=T)
# save(file = "data/fulcrum/exif.Rda", exif)

# load data from images already processed by Exif
#load("data/fulcrum/exif.Rda")

# load exif data exported directly from Fulcrum
exif <- readr::read_csv("data/fulcrum/nematode_field_sampling_sample_photo.csv") %>%
  dplyr::select(fulcrum_id, exif_gps_latitude, exif_gps_longitude, exif_gps_altitude)

###################################################################
### 4: Join collection, isolation, and location data            ###
###################################################################

#prevent scientific notation
options(scipen=999)

# join collection, isolation, and location data
df1 <- dplyr::full_join(isolation, collection, by = c("c_label_id" = "fulcrum_id")) %>%
  #rename the lat and long from fulcrum to collection_fulcrum_latitude and collection_fulcrum_longitude so that we can specify lat and long from exif tool
  dplyr::rename(collection_fulcrum_latitude = latitude, collection_fulcrum_longitude = longitude) %>%
  dplyr::select(c_label,
                everything(),
                -c_label_id) %>% #-sample_photo_url: we used to remove this but it's helpful
  # Join position data from exif by sample_photo. in some cases there is not position data from the photos
  dplyr::left_join(exif, by = c("sample_photo" = "fulcrum_id")) %>%
  # Create flag to track if lat and long come from record or photo
  dplyr::mutate(collection_lat_long_method = ifelse(is.na(exif_gps_latitude), "fulcrum", "photo")) %>%
  # In cases where lat/lon are not available from photo set to collection_fulcrum_latitude and collection_fulcrum_longitude 
  dplyr::mutate(latitude = ifelse(is.na(exif_gps_latitude), collection_fulcrum_latitude, exif_gps_latitude)) %>%
  dplyr::mutate(longitude = ifelse(is.na(exif_gps_longitude), collection_fulcrum_longitude, exif_gps_longitude)) %>%
  #dplyr::mutate_at(.vars = vars(dplyr::starts_with("gps")),
  #                 .funs = funs(as.numeric)) %>%
  dplyr::rename(fulcrum_altitude = gps_altitude) %>% 
  dplyr::mutate(worms_on_sample = ifelse(is.na(worms_on_sample), "?", worms_on_sample)) %>%
  dplyr::filter(!is.na(c_label)) %>%
  # Calculate the Haversine distance between fulcrum record_latitude and record_longitue and photo latitude and longitude
  dplyr::rowwise() %>%
  dplyr::mutate(collection_lat_long_method_diff = geosphere::distHaversine(c(longitude, latitude),
                                                                           c(collection_fulcrum_longitude, collection_fulcrum_latitude)),
                # adjust collection_lat_long_method_diff to NA if there is only a fulcrum GPS postion
                collection_lat_long_method_diff = ifelse(collection_lat_long_method == "fulcrum", NA, collection_lat_long_method_diff)) %>%
  # rename the latitude and longitude to include 'collection_' prefix
  dplyr::ungroup() %>%
  dplyr::rename(collection_latitude = latitude,
                collection_longitude = longitude,
                collection_local_time = time)

####################################################################
###      (OPTIONAL) CORRECT DUPLICATE ISOLATIONS (OPTIONAL)      ###
####################################################################

# In 3 instances there were two separate isolation records for the same c_label.
# selecting to retain the isolation record based on worm presence. Yes > tracks > no.
# If both isolation records indicate "tracks only" or "no" then we retain the earliest record.
# In no cases did both isolation records indicate worms were present on c_label.
duplicated_isolations_to_remove <- c("3881a3ae-b7a5-409d-a6a9-83eb64403600",
                                     "c69a4084-95d0-42d3-93b1-ff7633c24dc5",
                                     "26214a32-854d-471f-896d-d28447ada398")
df1 <- df1 %>%
  dplyr::filter(!isolation_id %in% duplicated_isolations_to_remove)

###################################################################
### 5: Joining C_lables with S_labels                           ###
###################################################################

df2 <- readr::read_csv("data/fulcrum/nematode_isolation_s_labeled_plates.csv") %>%
  dplyr::select(fulcrum_parent_id, s_label) %>%
  dplyr::full_join(df1, by = c("fulcrum_parent_id" = "isolation_id")) %>% # this used to be a left join
  dplyr::select(-fulcrum_parent_id, -updated_by, -version, -geometry) 

# # OPTIONAL: remove duplicated s_label
# df2 <- df2 %>%
#   # add a count of row number to grouping variable to remove duplicate s_label (S-11690).
#   dplyr::group_by(s_label) %>%
#   dplyr::mutate(n = row_number()) %>%
#   dplyr::mutate(n = ifelse(is.na(s_label), NA, n)) %>%
#   dplyr::filter(is.na(n) | n == "1")

###################################################################
### 6: Handle Substrates                                        ###
###################################################################

# adjust substrate categories if needed

###################################################################
### 7: Get alititudes from ggs locations                        ###
###################################################################

# #only need to run once to get altitudes.
 # options(geonamesUsername="katiesevans")
 # altitudes <- df1 %>%
 #   dplyr::ungroup() %>%
 #   dplyr::select(c_label, collection_latitude, collection_longitude) %>%
 #   dplyr::rowwise() %>%
 #   # Use collection_latitidue and collection_longitude to find altitudes. Note, these lat and longs should be spot checked to ensure proper collection locations.
 #   dplyr::mutate(geonames_altitude = geonames::GNsrtm3(collection_latitude, collection_longitude)$srtm3) %>%
 #   dplyr::ungroup()

# save(altitudes, file = "data/fulcrum/altitude.Rda")

# join geonames altitude data to record altitude data
load("data/fulcrum/altitude.Rda")

df3 <- df2 %>%
  dplyr::ungroup() %>%
  dplyr::inner_join(altitudes %>% dplyr::select(c_label, geonames_altitude)) %>%
  tidyr::unite(., col = altitude_methods_range, c(exif_gps_altitude, geonames_altitude, fulcrum_altitude), sep = ", ", remove = FALSE, na.rm = TRUE) %>%
  # make altitude variable and altitude_method variables to track which altitude is being used
  dplyr::mutate(altitude = ifelse(!is.na(exif_gps_altitude), exif_gps_altitude,
                                  ifelse(is.na(exif_gps_altitude) & !is.na(geonames_altitude), geonames_altitude,
                                         ifelse(is.na(exif_gps_altitude) & is.na(geonames_altitude), fulcrum_altitude, NA))),
                altitude_method = ifelse(!is.na(exif_gps_altitude), "photo",
                                         ifelse(is.na(exif_gps_altitude) & !is.na(geonames_altitude), "geonames",
                                                ifelse(is.na(exif_gps_altitude) & is.na(geonames_altitude), "fulcrum", NA_character_))))
  
###################################################################
### 8: Assign Islands and trails                                ###
###################################################################
# temp df4 assignment for cases where no trails are given
df4 <- df3 %>%
  dplyr::mutate(collection_trail = NA,
                collection_island = NA,
                collection_location = NA)

###################################################################
### 9: Add photo URLs                                           ###
###################################################################

###################################################################
### 10: OPTIONAL Merge automated blast data                     ###
###################################################################

# # Merge in blast data; Take top hit
# blast_results <- readr::read_tsv("data/sanger/blast_results.tsv") %>%
#   dplyr::group_by(s_plate) %>%
#   dplyr::filter(row_number() == 1)

###################################################################
### 11: Merge manual blast results                              ###
###################################################################

# Each collection should have unique gs_key. Your Google Sheet key can be found in your PUBLISHED Google sheets URL.
# Select the string of data found between the slashes after spreadsheets/d in your Google Sheet URL.
# use stringr to find s_labels from s_label column

genotyping_sheet_raw <- googlesheets::gs_key("17qM6cH9TVAhh7DRQeufHUU5iG5rHCP7cwgEtRA50EkE") %>%
  googlesheets::gs_read("genotyping template", na = c("#N/A", "NA", ""),
                        by = c("s_label")) %>%
  dplyr::filter(!is.na(s_label)) %>%
  # remove c_label variable (this column was hand typed and contains at least 2 errors)
  dplyr::select(s_label, species_id, shipment_number, shipment_sent_date, shipment_received_date,
                shipment_color, proliferation_48, proliferation_168, proliferating, pcr_product_its2,
                pcr_product_ssu, notes, manual_blast_notes, ECA_dirty, ECA_clean, possible_new_caeno_sp, make_strain_name, reason_strain_not_named)

# find s_labels in genotyping sheet
slabels <- str_subset(genotyping_sheet_raw$s_label, pattern = "S-")

# filter genotyping sheet by s_labels matching "S-" pattern
genotyping_sheet <- genotyping_sheet_raw %>%
  dplyr::filter(s_label %in% slabels)

# Remove any duplicated s_labels.

# Join genotyping sheet with collection and isolation data
df5 <- df4 %>% 
  dplyr::full_join(genotyping_sheet) %>%
  # Rename variables
  dplyr::rename(project_id = project,
                collection_id = c_label,
                isolation_id = s_label) %>%
  # Fill project_id variable incase there are NAs introduced
  tidyr::fill(project_id) %>%
  # Reorder variables
  dplyr::select(project_id,
                collection_id,
                isolation_id,
                species_id,
                ECA_dirty,
                ECA_clean,
                collection_by,
                collection_datetime_UTC,
                collection_date_UTC,
                collection_local_time,
                collection_island,
                collection_location,
                collection_trail,
                collection_latitude,
                collection_longitude,
                collection_fulcrum_latitude,
                collection_fulcrum_longitude,
                collection_lat_long_method,
                collection_lat_long_method_diff,
                ambient_temperature,
                flag_ambient_temperature_run,
                ambient_humidity,
                substrate_temperature,
                fulcrum_altitude,
                geonames_altitude,
                altitude,
                altitude_method,
                altitude_methods_range,
                landscape,
                sky_view,
                substrate,
                substrate_other,
                substrate_notes,
                sample_photo_url,
                gridsect,
                gridsect_index,
                grid_sect_direction,
                gridsect_radius,
                isolation_by,
                isolation_datetime_UTC,
                isolation_date_UTC,
                isolation_local_time,
                isolation_latitude,
                isolation_longitude,
                worms_on_sample,
                proliferation_48,
                proliferation_168,
                proliferating,
                approximate_number_of_worms,
                shipment_number,
                shipment_sent_date,
                shipment_received_date,
                pcr_product_its2, #ITS2_pcr_product
                pcr_product_ssu, #rhabditid_pcr_product
                manual_blast_notes, 
                notes,
                possible_new_caeno_sp,
                make_strain_name,
                reason_strain_not_named)

###################################################################
### 13: Export processed data                                   ###
###################################################################
# temp assign df5 to fulcrum_dat
`2020FebruaryAustralia_fulcrum_v2` <- df5
# export R dataframe
save(file = "data/fulcrum/2020FebruaryAustralia_fulcrum_v2.Rda", '2020FebruaryAustralia_fulcrum_v2')

###################################################################
### 14: project specific report                                 ###
###################################################################