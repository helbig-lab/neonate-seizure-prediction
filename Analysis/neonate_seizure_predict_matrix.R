#' Neonate seizure prediction
#'
#' Explore and create matrix for predictive modeling
#' For all neonates OR HIE only
#'
#' Required Files:
#' eeg_pivot.csv

library(librarian)
librarian::shelf(ggplot2, tidyverse, scales, lubridate, readxl)

getwd()
setwd("Insert your file directory here")

### Load and clean data set for analysis ----
eeg_hie_dc <- read.csv("./eeg_pivot.csv")

# Remove all standard EEG data and update naming convention where necessary
warning("Check formatting in table before running
        as these naming conventions may change depending on software used to collect data.")

eeg_hie_dc <- eeg_hie_dc %>% filter(PROC_NAME != "EEG STANDARD (<60MIN)")
names(eeg_hie_dc) <- gsub(names(eeg_hie_dc), pattern = 'NEURO\\.', replacement = "NEURO#")
names(eeg_hie_dc) <- gsub(names(eeg_hie_dc), pattern = '\\.', replacement = " ")

# Create numeric list of unique MRNs
num_HIE_pts_DC <- eeg_hie_dc %>%
  dplyr::select(PAT_MRN_ID) %>%
  unique() %>%
  unlist() %>%
  as.numeric()

# Add monitoring session labels to each row
eeg_hie_dc <- eeg_hie_dc %>%
  group_by(PAT_MRN_ID) %>%
  arrange(AGE_EEG_DAYS) %>%
  dplyr::mutate(monitor_day = yday(as.Date(END_EXAM_DTTM, format = "%d-%B-%y"))) %>% # Convert dates to days
  dplyr::mutate(visit_group = cumsum(c(0, diff(monitor_day)) >= 0 & PROC_NAME != "EEG WITH VIDEO SUBSEQUENT DAY")) %>%  # Chunk into consecutive monitoring periods and label as such, starting from group label of '1', cumulating logical values (i.e. create new group whenever a 'subsequent monitoring day' has passed)
  ungroup()

### REDACTED #1 - START

source_section("neonate_redacted.R", "REDACTED #1 - START", "REDACTED #1 - END")

### REDACTED #1 - END

## Ensure monitoring sessions line up with procedure names
# (i.e. did any days start a monitoring session but were labeled with a misaligned date? or non-subsequent day? or anything else unusual?)
test <- eeg_hie_dc %>%
  dplyr::select(PAT_MRN_ID, AGE_EEG_DAYS, PROC_BGN_TIME, END_EXAM_DTTM, visit_group, PROC_NAME) %>%
  group_by(PAT_MRN_ID, visit_group) %>%
  arrange(PAT_MRN_ID, visit_group) %>%
  slice(1)

## Switch for HIE+TH only vs. all neonates
hie_switch <- readline(prompt = "Would you like to filter for neonates with HIE who underwent Therapeutic Hypothermia? (y/n): ")
hie_switch <- tolower(as.character(hie_switch))

if (hie_switch == "y") {
  warning("filtered for HIE_TH patients only, filtering out other neonates..")
  eeg_hie_dc <- eeg_hie_dc %>% filter(HIE == "Y" & THER_HYPO == "Y")
} else if (hie_switch == "n") {
  print("Creating matrix of all neonates...")
} else {
  print("This variable was not recognized, continuing analysis with all neonates")
}

### Initial Analysis ----

## Histogram of monitoring duration
dur_mon <- NULL
for (i in seq_along(num_HIE_pts_DC)) {
  tmp <- eeg_hie_dc[eeg_hie_dc$PAT_MRN_ID == num_HIE_pts_DC[i], ]
  dur_mon[i] <- NROW(tmp)
}
hist(dur_mon)
summary(dur_mon)

## Age of first hook-up (in days)
tmp <- eeg_hie_dc %>%
  group_by(PAT_MRN_ID) %>%
  summarise(min_EEG = min(AGE_EEG_DAYS))

summary(tmp$min_EEG)

ggplot(tmp, aes(x = min_EEG)) +
  geom_histogram()

### Features reported on day 1 ----
# This loop with either
# A. look at first encounter day only when dealing with HIE+TH specifically or
# B. look at first day of each "encounter chunk" when looking at all neonates
# (some patients will have multiple rows as a result)
day_1_features <- NULL
for (i in seq_along(num_HIE_pts_DC)) {
  tmp <- eeg_hie_dc[eeg_hie_dc$PAT_MRN_ID == num_HIE_pts_DC[i], ]
  if (hie_switch == "y") { # If filtered for only HIE+TH cohort
    day_1_features <- rbind(day_1_features, tmp[tmp$AGE_EEG_DAYS ==  min(tmp$AGE_EEG_DAYS), ])
  } else {
    # Attach anything that happened on day 1 (including "subsequent" eegs on that day)
    # Or any 1st day in a set (ex. day 5 of days 5-10).
    # Note this may still miss "subsequent day" data in day 1 of a set if not in day 1 (ex. day 5 of days 5-10)
    # If multiple entries appear for say day 1 of a session where age is the same...pick minimum EEG datetime
    tmp2 <- tmp %>%
      arrange(AGE_EEG_DAYS) %>%
      group_by(visit_group) %>%
      slice(1) %>%
      ungroup() # First day of first monitoring session and first day of any new monitoring sessions

    if (NROW(tmp2) > 1)  {
      tmp2 <- tmp2[duplicated(tmp2$AGE_EEG_DAYS) == FALSE, ]
      day_1_features <- rbind(day_1_features, tmp2)
    } else {
      day_1_features <- rbind(day_1_features, tmp2)
    }

  }

}

## Quick checks of day_1_features dataset
test <- eeg_hie_dc[eeg_hie_dc$AGE_EEG_DAYS !=  min(eeg_hie_dc$AGE_EEG_DAYS) & eeg_hie_dc$PROC_NAME == "EEG WITH VIDEO SUBSEQUENT DAY" & eeg_hie_dc$PAT_MRN_ID %in% day_1_features$PAT_MRN_ID, ]
test2 <- day_1_features %>%
  group_by(PAT_MRN_ID) %>%
  summarise(n = n())

# Test that monitoring sessions are grouped properly for select patient tests
test3 <- eeg_hie_dc %>%
  filter(PAT_MRN_ID == "your patient ID")  %>%
  dplyr::mutate(`CHOPNEURO#352_EEG SEIZURES` = case_when(
    grepl("EEG-only seizures", `CHOPNEURO#352_EEG SEIZURES`) ~ 1,
    grepl("electroclinical seizures", `CHOPNEURO#352_EEG SEIZURES`) ~ 1,
    grepl("clinically only seizures", `CHOPNEURO#352_EEG SEIZURES`) ~ 1,
    grepl("none", `CHOPNEURO#352_EEG SEIZURES`) ~ 0)) %>%
  filter(!ORDER_PROC_ID %in% day_1_features$ORDER_PROC_ID[day_1_features$PAT_MRN_ID == "your patient ID"] & !AGE_EEG_DAYS %in% day_1_features$AGE_EEG_DAYS[day_1_features$PAT_MRN_ID == "your patient ID"]) %>%
  arrange(AGE_EEG_DAYS) %>%  # For matching to encounters later
  dplyr::select(c(PAT_MRN_ID, AGE_EEG_DAYS, END_EXAM_DTTM, `CHOPNEURO#352_EEG SEIZURES`)) %>%
  dplyr::mutate(monitor_day = yday(as.Date(END_EXAM_DTTM, format = "%d-%B-%y"))) %>% # Convert dates to days
  dplyr::mutate(visit_group = cumsum(c(0, diff(monitor_day)) >= 0 & PROC_NAME != "EEG WITH VIDEO SUBSEQUENT DAY")) # Chunk into consecutive monitoring periods and label as such, starting from group label of '1', accumulating logical values (i.e. create new group when new monitoring started)


### Seizures on Day 2 or higher ----
num_HIE_encs_DC <- eeg_hie_dc %>%
  dplyr::select(PAT_MRN_ID, visit_group) %>%
  unique() %>%
  arrange(PAT_MRN_ID, visit_group) # By monitoring session AND patient to capture patients with multiple sessions

subday_sz <- NULL
ID_tab <- NULL
for (i in seq_len(NROW(num_HIE_encs_DC))) {

  tmp <- eeg_hie_dc %>%
    filter(PAT_MRN_ID == num_HIE_encs_DC$PAT_MRN_ID[i], visit_group == num_HIE_encs_DC$visit_group[i])  %>%
    dplyr::mutate(`CHOPNEURO#352_EEG SEIZURES` = case_when(
      grepl("EEG-only seizures", `CHOPNEURO#352_EEG SEIZURES`) ~ 1,
      grepl("electroclinical seizures", `CHOPNEURO#352_EEG SEIZURES`) ~ 1,
      grepl("clinically only seizures", `CHOPNEURO#352_EEG SEIZURES`) ~ 1,
      grepl("none", `CHOPNEURO#352_EEG SEIZURES`) ~ 0))

  day1tmp <- day_1_features %>%
    filter(PAT_MRN_ID == num_HIE_encs_DC$PAT_MRN_ID[i], visit_group == num_HIE_encs_DC$visit_group[i])

  tmp <- tmp %>%
    filter(!ORDER_PROC_ID %in% day1tmp$ORDER_PROC_ID & !AGE_EEG_DAYS %in% day1tmp$AGE_EEG_DAYS) %>%  # For when day1 has multiple entries, filter out those rows
        arrange(AGE_EEG_DAYS)

  # If there are no subsequent days of EEGs
  if (NROW(tmp) == 0) {

    subday_sz[i] <- NA
    ID_tab[i] <- paste0(num_HIE_encs_DC$PAT_MRN_ID[i], "-", num_HIE_encs_DC$visit_group[i])

  }
  # If more than one row and all subsequent rows are not NA and are not already rows in the day1 features table...
  # Arrange by monitoring session and attach the 'maximum' subsequent seizure status during that session
  else if (NROW(tmp) >= 1 & any(!is.na(tmp$`CHOPNEURO#352_EEG SEIZURES`[seq_len(NROW(tmp))]))) {

    subday_sz[i] <- tmp[seq_len(NROW(tmp)), ] %>%
      dplyr::select(`CHOPNEURO#352_EEG SEIZURES`) %>%
      max(na.rm = TRUE)  # Record if there was subsequent seizure data after day 1
    ID_tab[i] <- paste0(num_HIE_encs_DC$PAT_MRN_ID[i], "-", num_HIE_encs_DC$visit_group[i])

  } else {
    subday_sz[i] <- 0 # If there are no other values (i.e. NAs only) on subsequent days, then use 0 (i.e. assume subsequent EEGs did not show seizures)
    ID_tab[i] <- paste0(num_HIE_encs_DC$PAT_MRN_ID[i], "-", num_HIE_encs_DC$visit_group[i])

  }
}

### Create subsequent seizure + MRN master table and additional mini-analyses ----
subday_sz_ID <- as.data.frame(cbind(subday_sz, ID_tab))

## Gender breakdown - at the patient level
male <- NROW(unique(day_1_features$PAT_MRN_ID[day_1_features$SEX == "Male"]))
female <- NROW(unique(day_1_features$PAT_MRN_ID[day_1_features$SEX == "Female"]))
gender <- c(male, female)
gender_lbs <- c("Male", "Female")
pie(gender, gender_lbs)

## Race Breakdown - at the patient level
race <- day_1_features %>%
  dplyr::select(PAT_MRN_ID, PATIENT_RACE) %>%
  distinct() %>%
  group_by(PATIENT_RACE) %>%
  summarise(n = n()) %>%
  filter(!is.na(PATIENT_RACE), PATIENT_RACE != "Native Hawaiian or Other Pacific Islander")

race_lbs <- c("Black", "White", "Asian", "Multiple Races", "Other")
pie(race$n, race_lbs)

## Ethnicity Breakdown - at the patient level
ethnicity <- day_1_features %>%
  dplyr::select(PAT_MRN_ID, ETHNICITY) %>%
  distinct() %>%
  group_by(ETHNICITY) %>%
  summarise(n = n()) %>%
  filter(!is.na(ETHNICITY))

ethnicity_lbs <- c("Hispanic or Latino", "Not Hispanic or Latino")
pie(ethnicity$n, ethnicity_lbs)

## Check fill rates for all fields for potential additional features
temp <- day_1_features

min(
NROW(na.omit(temp$SEX)) / NROW(temp),
NROW(na.omit(temp$PATIENT_RACE)) / NROW(temp),
NROW(na.omit(temp$`CHOPNEURO#722_EEG GRAPHOELEMENTS`)) / NROW(temp),
NROW(na.omit(temp$`CHOPNEURO#716_EEG VOLTAGE`)) / NROW(temp),
NROW(na.omit(temp$`CHOPNEURO#695_EEG CONTINUITY NEONATE`)) / NROW(temp),
NROW(na.omit(temp$`CHOPNEURO#729_EEG TRANSIENTS`)) / NROW(temp),
NROW(na.omit(temp$`CHOPNEURO#567_EEG VARIABILITY`)) / NROW(temp),
NROW(na.omit(temp$`CHOPNEURO#568_EEG REACTIVITY`)) / NROW(temp),
NROW(na.omit(temp$`CHOPNEURO#146_EEG IMPRESSION`)) / NROW(temp),
NROW(na.omit(temp$`CHOPNEURO#352_EEG SEIZURES`)) / NROW(temp)
)

# Min frquency is ~95% for full neonate cohort. This will be the limit we set for other fields
feature_list <- list()
for (i in colnames(temp)) {

  if (NROW(na.omit(temp[, i])) / NROW(temp) >= 0.95) {

    print(i)
    print(NROW(na.omit(temp[, i])) / NROW(temp))
    feature_list <- append(i, feature_list)
  }
}

table(temp$`CHOPNEURO#843_EEG LOCATION`)
table(temp$`CHOPNEURO#559_EEG PATIENT STATE`)
table(temp$`CHOPNEURO#409_EEG EVENTS`)
table(temp$HIE)

## Converting specific fields to binary
day_1_features1 <- day_1_features %>% dplyr::mutate(SEX = recode(SEX,
                                                      "Male" = 1,
                                                      "Female" = 0))

day_1_features1 <- day_1_features1 %>%
  dplyr::mutate(`CHOPNEURO#352_EEG SEIZURES` = case_when(
    grepl("EEG-only seizures", `CHOPNEURO#352_EEG SEIZURES`) ~ as.numeric(1),
    grepl("clinically only seizures", `CHOPNEURO#352_EEG SEIZURES`) ~ as.numeric(1), # Confirm with any new patients this counts as seizures
    grepl("electroclinical seizures", `CHOPNEURO#352_EEG SEIZURES`) ~ as.numeric(1),
    grepl("none", `CHOPNEURO#352_EEG SEIZURES`) ~ as.numeric(0)))


## Summarize using count for general features overview
day_1_features %>% count(`CHOPNEURO#352_EEG SEIZURES`)
day_1_features %>% count(`CHOPNEURO#849_EEG ICU CARE`)
day_1_features %>% count(`CHOPNEURO#722_EEG GRAPHOELEMENTS`)
day_1_features %>% count(`CHOPNEURO#716_EEG VOLTAGE`)
day_1_features %>% count(`CHOPNEURO#844_EEG INDICATIONS`)
day_1_features %>% count(`CHOPNEURO#561_EEG PREDOMINANT BACKGROUND FREQUENCIES`)
day_1_features %>% count(`CHOPNEURO#409_EEG EVENTS`)
day_1_features %>% count(`CHOPNEURO#1075_EEG MED CHANGES`)
day_1_features %>% count(`CHOPNEURO#720_EEG DYSMATURITY`)
day_1_features %>% count(`CHOPNEURO#774_EEG IMPRESSION DYSFUNCTION`)
day_1_features %>% count(`CHOPNEURO#697_EEG TYPICAL INTERBURST VOLTAGE`)
day_1_features %>% count(`CHOPNEURO#773_EEG REACTIVITY CHANGE`)
day_1_features %>% count(`CHOPNEURO#695_EEG CONTINUITY NEONATE`)
day_1_features %>% count(`CHOPNEURO#729_EEG TRANSIENTS`)
day_1_features %>% count(`CHOPNEURO#567_EEG VARIABILITY`)
day_1_features %>% count(`CHOPNEURO#568_EEG REACTIVITY`)
day_1_features %>% count(`CHOPNEURO#597_EEG REACTIVITY STIMULUS`)
day_1_features %>% count(`CHOPNEURO#711_EEG TYPICAL INTERBURST DURATION`)
day_1_features %>% count(`CHOPNEURO#146_EEG IMPRESSION`)
day_1_features %>% count(`CHOPNEURO#559_EEG PATIENT STATE`)
day_1_features %>% count(`CHOPNEURO#665_EEG SEIZURE RECOGNITION`)
day_1_features %>% count(`CHOPNEURO#693_EEG PATIENT STATE SLEEP TYPE`)
day_1_features %>% count(`CHOPNEURO#1033_EEG ICU DRUGS`)


##  Simplify matrix to variables of interest
# Note: MRNs to be removed for matrix analysis.  Only kept to merge with sub_seiz_days
day_1_matrix <- day_1_features1[, c("PAT_MRN_ID",
                                 "visit_group",
                                 "END_EXAM_DTTM",
                                 "SEX",
                                 "PATIENT_RACE",
                                 "CHOPNEURO#722_EEG GRAPHOELEMENTS",
                                 "CHOPNEURO#716_EEG VOLTAGE",
                                 "CHOPNEURO#695_EEG CONTINUITY NEONATE",
                                 "CHOPNEURO#729_EEG TRANSIENTS",
                                 "CHOPNEURO#567_EEG VARIABILITY",
                                 "CHOPNEURO#568_EEG REACTIVITY",
                                 "CHOPNEURO#146_EEG IMPRESSION",
                                 "CHOPNEURO#352_EEG SEIZURES",
                                 "HIE", # OPTIONAL, RUN TEST BEFORE USING
                                 "THER_HYPO", # OPTIONAL, RUN TEST BEFORE USING
                                 "CHOPNEURO#409_EEG EVENTS" # OPTIONAL, RUN TEST BEFORE USING
                                 )]

## Combining day1 with subsequent seizures
subday_sz1 <- subday_sz_ID
subday_sz1 <- subday_sz1 %>%
  separate(ID_tab, into = c("PAT_MRN_ID", "visit_group"), sep = "-") %>%
  dplyr::mutate(PAT_MRN_ID = as.integer(PAT_MRN_ID),
                visit_group = as.integer(visit_group))

day_1_matrix <- left_join(day_1_matrix,
                          subday_sz1,
                          by = c("PAT_MRN_ID", "visit_group"))

## Save raw matrix
day_1_matrix_raw <- day_1_matrix

## Similar to above, but for subsequent days
day_1_matrix_sub_sz_raw <- subday_sz

## Convert remaining fields to binary

day_1_matrix <- day_1_matrix %>%
  dplyr::mutate(PATIENT_RACE = recode(PATIENT_RACE,
                                     `White` = 0,
                                     `Black` = 1,
                                     `Asian` = 1,
                                     `Other` = 1,
                                     `Multiple Races` = 1,
                                     `American Indian or Alaska Native` = 1,
                                     `Native Hawaiian or Other Pacific Islander` = 1),
                `CHOPNEURO#722_EEG GRAPHOELEMENTS` = recode(`CHOPNEURO#722_EEG GRAPHOELEMENTS`,
                                      `absent` = 1,
                                      `present` = 0,
                                      `unknown` = 0),
                `CHOPNEURO#716_EEG VOLTAGE` = case_when(
                                      grepl("low voltage suppressed", `CHOPNEURO#716_EEG VOLTAGE`) ~ 1,
                                      grepl("borderline low", `CHOPNEURO#716_EEG VOLTAGE`) ~ 1,
                                      grepl("asymmetry", `CHOPNEURO#716_EEG VOLTAGE`) ~ 1,
                                      grepl("high", `CHOPNEURO#716_EEG VOLTAGE`) ~ 1,
                                      grepl("normal", `CHOPNEURO#716_EEG VOLTAGE`) ~ 0), # Note if any of the above are combined with normal, defaults to not normal (ex. normal/borderline low = 1)
                `CHOPNEURO#695_EEG CONTINUITY NEONATE` = case_when(
                                      grepl("low voltage suppressed", `CHOPNEURO#695_EEG CONTINUITY NEONATE`) ~ 1,
                                      grepl("excessive discontinuity", `CHOPNEURO#695_EEG CONTINUITY NEONATE`) ~ 1,
                                      grepl("burst supression", `CHOPNEURO#695_EEG CONTINUITY NEONATE`) ~ 1,
                                      grepl("asymmetry", `CHOPNEURO#695_EEG CONTINUITY NEONATE`) ~ 1,
                                      grepl("normal continuity", `CHOPNEURO#695_EEG CONTINUITY NEONATE`) ~ 0, # Note if any of the above are combined with normal, defaults to not normal)
                                      grepl("normal discontinuity", `CHOPNEURO#695_EEG CONTINUITY NEONATE`) ~ 0),
                `CHOPNEURO#729_EEG TRANSIENTS` = recode(`CHOPNEURO#729_EEG TRANSIENTS`,
                                      `sharp wave` = 1,
                                      `sharp wave (1 type)` = 1,
                                      `sharp wave (2 types)` = 1,
                                      `sharp wave (3 types)` = 1,
                                      `sharp wave (4 types)` = 1,
                                      `abn sharp waves` = 1,
                                      `abn sharp waves (1 type)` = 1,
                                      `abn sharp waves (2 types)` = 1,
                                      `abn sharp waves (3 types)` = 1,
                                      `absent` = 0),
                `CHOPNEURO#567_EEG VARIABILITY` = recode(`CHOPNEURO#567_EEG VARIABILITY`,
                                      `absent` = 1,
                                      `unknown/unclear/not applicable` = 0,
                                      `present` = 0)
                )

### REDACTED #2 - START

source_section("neonate_redacted.R", "REDACTED #2 - START", "REDACTED #2 - END")

### REDACTED #2 - END

day_1_matrix <- day_1_matrix %>%
  dplyr::mutate(`CHOPNEURO#568_EEG REACTIVITY` = recode(`CHOPNEURO#568_EEG REACTIVITY`,
                                      `absent` = 1,
                                      `unclear` = 0,
                                      `not tested` = 0,
                                      `present` = 0),
                `CHOPNEURO#146_EEG IMPRESSION` = recode(`CHOPNEURO#146_EEG IMPRESSION`,
                                      `abnormal` = 1,
                                      `indeterminate` = 1,
                                      `normal` = 0),
                # OPTIONAL, RUN TEST BEFORE USING
                HIE_TH = case_when(HIE == "Y" & THER_HYPO == "Y" ~ 1,
                                                          TRUE ~ 0),
                #OPTIONAL, RUN TEST BEFORE USING
                `CHOPNEURO#409_EEG EVENTS` = recode(`CHOPNEURO#409_EEG EVENTS`,
                                      `occurred without cerebral electrographic correlate` = 1,
                                      `none` = 0)
                )



## For organizing patients by date of procedure for iterable analyses:
day_1_matrix <- day_1_matrix %>%
  dplyr::mutate(END_EXAM_DTTM = as.Date(END_EXAM_DTTM, format = "%d-%B-%y")) %>%
  arrange(END_EXAM_DTTM) %>%
  dplyr::mutate(month = lubridate::month(END_EXAM_DTTM)) %>%
  dplyr::mutate(year = lubridate::year(END_EXAM_DTTM),
                sex = SEX,
                race = PATIENT_RACE,
                graphoelements = `CHOPNEURO#722_EEG GRAPHOELEMENTS`,
                voltage = `CHOPNEURO#716_EEG VOLTAGE`,
                continuity = `CHOPNEURO#695_EEG CONTINUITY NEONATE`,
                transients = `CHOPNEURO#729_EEG TRANSIENTS`,
                variability = `CHOPNEURO#567_EEG VARIABILITY`,
                reactivity = `CHOPNEURO#568_EEG REACTIVITY`,
                impression = `CHOPNEURO#146_EEG IMPRESSION`,
                seizures = `CHOPNEURO#352_EEG SEIZURES`,
                hie_th = `HIE_TH`, # OPTIONAL, RUN TEST BEFORE USING
                events  = `CHOPNEURO#409_EEG EVENTS` # OPTIONAL, RUN TEST BEFORE USING
  ) %>%
  dplyr::select(-c(PAT_MRN_ID, visit_group, HIE, THER_HYPO)) # Not necessary in matrix form

# Matrix with any seizures as last column [converting NA to 0]
subday_sz <- as.numeric(as.character(day_1_matrix$subday_sz))
subday_sz[is.na(subday_sz)] <- 0

any_sz <- NULL
for (i in seq_len(NROW(day_1_matrix))) {

  tmp <- day_1_matrix %>%
    dplyr::select(c(seizures, subday_sz)) %>%
    dplyr::mutate(subday_sz = as.numeric(as.character(subday_sz))) %>%
    replace(is.na(.), 0)

  if (tmp$seizures[i] == 1 | tmp$subday_sz[i] == 1) {
    any_sz[i] <- TRUE
  } else {
    any_sz[i] <- FALSE
  }

}
day_1_matrix_anysz <- cbind(day_1_matrix, any_sz)

### Heatmap with supervised clustering (subsequent day seizures vs. not) ----
# Note: these packages may no longer be available in R 3.6.0+, so commented out query for now

# if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# BiocManager::install("ComplexHeatmap")
# library(ComplexHeatmap)

# df <- as.matrix(day_1_matrix)

# Heatmap(df)

# df <- subset(df, select = -c(subday_sz, race, sex))
# df <- day_1_matrix[complete.cases(day_1_matrix),] # Remove all NAs
# Heatmap(df, name = "subsequent sz", split = df$subday_sz,
#         column_title = "Variables", #row_title = "Samples",
#         row_names_gp = gpar(fontsize = 7))
# show_row_hclust = TRUE, show_column_hclust = TRUE)
