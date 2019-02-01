### Readme for processing of RNA seq survival data
### Author: Dr Marion Mateos
### Date: Jan 30 2019

### Summary of RNA seq biology project
### Hypothesis: 
### 1. that the transcriptome using RNAseq data provides additional prognostic information beyond clinicopathological risk features in current clinical practice
### 2. that there are novel transcripts that provide additional prognostic information beyond clinicopathological risk features in current clinical practice

######################################################################################################################################################
### MAIN INPUT
### need input for pData (phenotypic data, derived from Newcastle Medulloblastoma MB cohort "NMB")
### input: /home/nmm199/R/MB_RNAseq/Input_data/database270617_060917.csv

### pData_input_master.R
### test.pData file generated via pData_input_master.R.   If the clinical database file is updated then need to rerun pData_input_master.R

######################################################################################################################################################
### MAIN FUNCTIONS
### use clinical_data_master with clinical_data_functions_master.R file 
### as per notes in clinical_data_master, need input files as listed including test.pData file, GTF file

### clinical database "x.data" (from test.pData)
### 7 molecular group data "meth.data"
### cytogenetic arm data "cytogen.data"
### RNA expression data "RNA.data"


### main function is clinPathAssess
### generates results.master

######################################################################################################################################################
### MAIN OUTPUTS
### use clinical_data_extract_DW.R file with clinical_data_functions_extract_master.R
### input: need results.master
### this will then output individual results objects, dataframes and/or .csv depending on what is hashed out


######################################################################################################################################################
### OTHER FILES:
### FILES TO GENERATE NOVEL TRANSCRIPT DATA AND TO FURTHER CHARACTERISE NOVEL TRANSCRIPT
### novel_trans_v3.R
### novel_gtf.R

### FILES TO VALIDATE TRANSCRIPTS IN CAVALLI DATASET:
### validate_function.R
### validate_transcript.R

### GRAPHICS SCRIPTS
### clinical_data_graphics.R  for the histogram/adjusted p value plots
### Number_risk_KM.R  generates number at risk KM survival curves for individual goi of interest
### KM_number_risk_functions.R  accompanying function file for KM survival curves

