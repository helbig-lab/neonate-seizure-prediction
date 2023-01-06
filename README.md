# Neonate Seizure Prediction

**Project Summary:** To leverage neonatal EEG data gathered through routine care to generate neonatal seizure prediction models, both for the overall cohort of all neonates and the subgroup with HIE treated with therapeutic hypothermia.

## Before Using:

- In order to allow for additional customizability, scripts `neonate_seizure_predict_matrix.R` and `model_FPR_bootstrap.R` require user input to run.  These lines can be replaced with default values if desired and _should_ be replaced with default values if either script is sourced in relation to `seizure_predict_models.R`.
- There are several redacted sections of the `neonate_seizure_predict_matrix.R` script which previously contained PHI. In order to remove this identifiable information, scripts will not entirely replicate our original analyses.
- When running `seizure_predict_models.R`, `model_FPR_bootstrap.R` and other scripts, be sure to set file paths for any models you wish to save (currently, these paths default to the fake path names such as "your_path_here") and file names for any specific model results that will need to be uploaded.
- Generating H2o models tends to use a lot of processing power and take a while to run.  It is recommended that these models be run on a virtual machine or similar alternative to your local computer. Note that this may require changing the code in source file `neonate_seizure_predict_matrix.R`.

## How to generate data:

### If beginning from scratch with custom data:

**Note: Unless internal to CHOP, provided explicit permission or working with your own dataset, only de-identified data will be provided (see `Data` folder for de-identified modeling CSV files).  In which case, skip to "Modeling and analyses"**

1. Go to the `Extraction` folder and run `EEG_ANALYSIS_HIE.sql` in a compatible database to generate cohort data and save as `eeg_hie.csv`.  All remaining analyses will use the `Analysis` folder.
2. Open the `HIE_patient_extraction.R` script, choose your directory and run the script to generate a reviewable list of individual with HIE who may have been treated using Therapeutic Hypothermia.
3. Manually review the list, saving any patient identifiers (ex. MRNs) to an Excel or CSV file that can be used to distinguish those with HIE who underwent Therapeutic Hypothermia.
4. Open the `EEG_neonate_pivot.R` file, choose your directory and run, using the manually reviewed list and initial raw data to generate a condensed pivot table of cohort data.
5. Once a "pivoted" dataset is created or provided, open the `neonate_seizure_predict_matrix.R`, choose your directory and run the script to do basic analysis of the dataset and to transform into a machine learning model ready dataset. Note that if you simply need the machine learning ready data, this file is already sourced to run in other scripts to generate this type of modeling dataset.

## Modeling and analyses:

- Once data is created (or pre-generated, machine-learning-ready files are uploaded or sourced), several scripts exist to both create, run, test and analyze data. The main source is `seizure_predict_models.R` which generates all singular models and their respective performance metrics using either a premade machine learning ready file, or sourcing data from `neonate_seizure_predict_matrix.R`.
- `model_figures_and_tables.R` can be used to generate basic graphs and tables related to model performance.
- For additional performance metrics, including False Positive Rate and Confidence Intervals for various performance metrics, use the `model_FPR_bootstrap.R` script.

## References:
Below are a list of sources that inspired or guideded the creation of several of these scripts

- https://afit-r.github.io/random_forests
- https://www.datacamp.com/tutorial/bootstrap-r
- https://statistics.berkeley.edu/sites/default/files/tech-reports/666.pdf

---

Use of this software is available to academic and non-profit institutions for research purposes subject to the terms of the 2-Clause BSD License (see copy below). For use or transfers of the software to commercial entities, please inquire with MCKEEJ@chop.edu. © 2022 CHOP.
The FreeBSD Copyright

Copyright and License Information

Copyright (C) 2022 The Children’s Hospital of Philadelphia

Authors: Jillian McKee, Michael Kaufman, Alexander Gonzalez
