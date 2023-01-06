#' HIE Patient Extraction
#'
#' This script takes the raw eeg cohort dataset and
#' extracts individuals whose diagnoses indicate HIE.
#' The created list can then be used to manually review charts
#' for potential use of therapeutic hypothermia to create patient flags.
#'
#' Required files: 
#' eeg_hie.csv

library(librarian)
librarian::shelf(tidyverse, dplyr)

getwd()
setwd("Insert working directory here...")

# Extract patient information from raw data
raw_list <- read.csv("./eeg_hie.csv")

# Find all individuals with HIE related diagnoses to check for Therapeutic hypothermia
# and update a running list with any new individuals
filtered_list <- raw_list %>%
  filter(grepl("P91.60|P91.61|P91.62|P91.63", CODE)) %>%
  select(c("PAT_MRN_ID", "ORDERING_DATE", "END_EXAM_DTTM", "CODE")) %>%
  distinct()

write.csv(filtered_list, "EEG_HIE_TH_LIST.csv")
