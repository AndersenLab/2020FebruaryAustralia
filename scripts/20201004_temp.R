library(tidyverse)

# set working directory
setwd(glue::glue("{dirname(rstudioapi::getActiveDocumentContext()$path)}/.."))

# add temporary public url for report map
df <- readr::read_rds("data/processed/fulcrum/20200828_2020FebruaryAustralia.rds") %>%
  dplyr::mutate(pub_url = glue::glue("https://storage.googleapis.com/elegansvariation.org/photos/isolation/fulcrum/2020FebruaryAustralia/sampling_thumbs/{strain_name}.jpg"))

# write temp

rio::export(df, "data/processed/fulcrum/temp_df.rds")
