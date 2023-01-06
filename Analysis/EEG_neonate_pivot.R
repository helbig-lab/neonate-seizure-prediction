#' EEG-HIE data pivot
#'
#' This script takes the raw neonate cohort dataset and
#' adds a HIE and TH flags.
#' The remainder of the script is used to collapse and
#' then pivot data per individual
#'
#' Required files:
#' eeg_hie.csv,
#' EEG_HIE_TH_LIST.xlsx (generated in part using "HIE_patient_extraction.R")

library(librarian)
librarian::shelf(tidyverse, dplyr, readxl)

getwd()
setwd("Insert working directory here...")

## Load initial dataset and 
## add new columns for individuals with HIE who underwent therapeutic hypothermia from hie_th_list
eeg_hie <- read.csv("./eeg_hie.csv")

head(eeg_hie)
str(eeg_hie)

hie_th_list <- read_excel("./EEG_HIE_TH_LIST.xlsx", sheet = "patient_list")

# Ensure all new MRNs already included in master table
test <- hie_th_list %>%
  filter(!PAT_MRN_ID %in% eeg_hie$PAT_MRN_ID) 

## Aggregate multiple entries for same diagnosis at same encounter per individual
## And merge back into original table
temp <- eeg_hie %>%
  select(PAT_MRN_ID, PAT_ENC_CSN_ID, CODE) %>%
  distinct() %>%
  group_by(PAT_MRN_ID, PAT_ENC_CSN_ID) %>%
  summarise(CODE = paste(CODE, collapse = " | "))

eeg_hie2 <- merge(eeg_hie,
                  temp,
                  by = c("PAT_MRN_ID", "PAT_ENC_CSN_ID")
                  ) %>%
  rename(CODE = CODE.y) %>%
  select(-c(DIAGNOSIS, CODE.x))


## Categorize diagnoses and class for "patient who has ever received T.H."
eeg_hie3 <- eeg_hie2 %>%
  mutate(HIE = ifelse(
    grepl("P91.60|P91.61|P91.62|P91.63", CODE), "Y", "N")
  ) %>%
  mutate(THER_HYPO = case_when)
    PAT_MRN_ID %in% hie_th_list$PAT_MRN_ID ~ "Y",
    TRUE ~ "N"
  ))

table(eeg_hie3$HIE, eeg_hie3$THER_HYPO) # Check number of entries

## Similar to above aggregation, but now aggregating element values
# Basically "all other elements being the same, merge duplicate elements with different values"
temp2 <- eeg_hie3 %>%
  select(PAT_MRN_ID, PAT_ENC_CSN_ID, ORDER_PROC_ID, ELEMENT_ID, SDE_CONCEPT_NAME, SMRTDTA_ELEM_VALUE) %>%
  distinct() %>%
  group_by(across(c(-SMRTDTA_ELEM_VALUE))) %>%
  summarise(SDE_VALUE = paste(SMRTDTA_ELEM_VALUE, collapse = " | "))

eeg_hie4 <- merge(eeg_hie3,
                  temp2,
                  by = c("PAT_MRN_ID", "PAT_ENC_CSN_ID", "ORDER_PROC_ID", "ELEMENT_ID", "SDE_CONCEPT_NAME")) %>%
  select(-c(SMRTDTA_ELEM_VALUE)) %>%
  distinct()

## Pivot Smart Data Elements from rows to columns
pivot <- eeg_hie4 %>%
  group_by(PAT_MRN_ID, PAT_ENC_CSN_ID) %>%
  mutate(row_num = str_c("V", row_number())) %>%
  ungroup() %>%
  pivot_wider(names_from = c(ELEMENT_ID, SDE_CONCEPT_NAME), values_from = SDE_VALUE) %>%
  select(-c(row_num)) %>%
  arrange(PAT_MRN_ID, PAT_ENC_CSN_ID, ORDER_PROC_ID, ORDERING_DATE)

## In each column, summarise by the grouping function, using the first data point that is not NA and
## return it as the "summarised" value in that column
pivot2 <- pivot

pivot_grouped <- pivot2 %>%
  group_by(PAT_MRN_ID, PAT_ENC_CSN_ID, ORDER_PROC_ID) %>%
  summarise_each(funs(first(.[!is.na(.)])))

write.csv(pivot_grouped, "eeg_pivot.csv", row.names = FALSE)
