### File introduction
### File name: Biomarker_assess_09_08_2017.R

### Aim of file is to 
# 1. Run basic descriptive statistics on a cohort of children treated for medulloblastoma, whose details are contained within the local clinical database
# 2. Analyse genes of interest in relation to univariate and multivariate risk prediction models for survival (overall, event-free and progression-free)
# 3. This script covers analysis up to and including univariate risk factor analysis
# 4. Multivariate analysis / AUC will be covered by a separate script


### Author: Dr Marion Mateos
### Date: July 3 2017

### Subsequent updates with input from Dr Louise Pease, based on clinical_data_4.R

### R version 3.4.0 (2017-04-21)
### Platform: x86_64-pc-linux-gnu (64-bit)
### Running under: Ubuntu 16.04.2 LTS

# attached base packages:
#   [1] parallel  stats     graphics  grDevices utils     datasets  methods   base     

### Packages and version numbers
## [1] survival_2.41-3     
## RColorBrewer_1.1-2  
## car_2.1-4           
## gplots_3.0.1        
## NMF_0.20.6          
## Biobase_2.36.2      
## BiocGenerics_0.22.0
## cluster_2.0.6       
## rngtools_1.2.4      
## pkgmaker_0.22       
## registry_0.3       


### Libraries required

library(NMF)
#source("http://bioconductor.org/biocLite.R")
#biocLite(c("gplots", "car"))
library('gplots')
library(car)
library(stats)
library('survival')
# library(scales)



### Names of functions used:
#source (file = "/home/nlp71/Marion/biomarker_discovery_functions.R")
#source (file = "/home/nmm199/R/MB_RNAseq/Clinical/clin.script/Current/biomarker_discovery_functions.R")


### "chi.sq"
### "cor.result"
### "lin.reg"
### "km.log.test"
### "km.log.test.OS"
### "cox.result.OS"
### "km.log.test.EFS"


### Will need to run "biomarker_discovery_functions.R" file and "prepare_pData_table.R" files 
### can choose to load an individual gene/transcript (goi) or a list of genes/transcripts. See instructions for unhashing relevant sections below
### ideally in initial biomarker discovery phase, will unhash sink and PDF outputs, due to large volume of data that will be generated

###################################################################################################

#source (file = "/home/nmm199/R/MB_RNAseq/Clinical/clin.script/Current/prepare_pData_table.R")

load("pData_saved.RData")

options(warn=-1)
options(browserNLdisabled = TRUE)
cox_hazard_results_df <- data.frame()
survival_results_df <- data.frame()
Chi_squared_results_df <- data.frame()
logistic_regression_results_df <- data.frame()
Chi_squared_sig_results_df <- data.frame()


#### Can read in a list of genes to test 
#gene_list <- read.csv(file ="genes.csv")

#### check a single gene of interest
gene_list <- as.list(c("ENSG00000000003.14_1", "ENSG00000000005.5_1", "ENSG00000000419.12_1"))

#gene_list <- as.list(c("ENSG00000001036", "ENSG00000001497", "ENSG00000002016"))
#### check a list of genes 
#gene_list <- as.list(c("ENSG00000136997", "ENSG00000197561", "ENSG00000000003", "ENSG00000124693"))
#gene_list <- as.list(c("ENSG00000136997", "ENSG00000197561", "ENSG00000000003", "ENSG00000150991", "ENSG00000101665", "ENSG00000116039", "ENSG00000182979", "ENSG00000133026", "ENSG00000123384", "ENSG00000000005", "ENSG00000000419", "ENSG00000000457", "ENSG00000000460", "ENSG00000000938"))

#### to check every gene in the RNAseq data unhash here (slow!)
#gene_list <- as.list(rownames(mb.vsd))
#gene_list <- gsub("\\..*","",gene_list)
#gene_list <- unique(gene_list)
#groups <- c("G3G4", "others")

for (gene in 1:length(gene_list)){
  #### unhash sink and pdf lines below to output a sink and pdf file per gene in the list
  #sink(paste("pDatalog",gene_list[[gene]],".txt"))
  #pdf(paste("marker.results",gene_list[[gene]],".pdf"), width = 10, height = 10)
  #### create an ordered list of variables measured in the matched data 
  index <- match(colnames(mb.vsd), rownames(test.pData))
  #goi.vsd <- as.numeric(mb.vsd[goi,]) 
  #### Generate a list of expression values for the gene of interest 
  Seq_goi.vsd <- as.numeric(mb.vsd[gene_list[[gene]],])
  names(Seq_goi.vsd) <- gsub("T","",names(mb.vsd))
  #### For a single gene use these, or simply only input one gene into the gene_list
  #goi.vsd <- colnames(mb.vsd[goi,])
  #goi.vsd <- as.numeric(mb.vsd[goi,])
  #index <- match(names(goi.vsd), rownames(test.pData)) 
  #### Since the matched data is based on sample ids in common between RNAseq and pheno data no need to amend this part
  matched.test.pData <- test.pData[index[!is.na(index)],] 
  matched.goi.vsd <- Seq_goi.vsd[!is.na(index)] 
  matched.goi.vsd.cat <- ifelse(matched.goi.vsd>median(Seq_goi.vsd, na.rm = T), "high","low")
  G3G4 <- matched.test.pData[which(matched.test.pData$group3or4fac =="G3G4"),]
  #### now assign these to names prefixed for grouping and append the gene identifier to each of the objects
  assign(paste0("Gene_Seq ",gene_list[[gene]]), Seq_goi.vsd)
  assign(paste0("Seq_gene_matched ",gene_list[[gene]]), matched.goi.vsd)
  assign(paste0("Cat_Seq_gene_matched ",gene_list[[gene]]), matched.goi.vsd.cat)
  assign(paste0("Matched_pDATA ",gene_list[[gene]]), matched.test.pData)
  assign(paste0("Matched_pDATA_G3G4 ",gene_list[[gene]]), G3G4)
  #### Use the prefixes generated above to place the objects in lists that were can itterate through
  CAT_SEQ_LIST <- as.list(mget(ls(pattern="^Cat_Seq_gene_matched")))
  SEQ_GENE_MATCHED <- as.list(mget(ls(pattern="^Seq_gene_matched")))
  GENE_SEQ_LIST <- as.list(mget(ls(pattern="^Gene_Seq")))
  MATCHED_TEST_PDATA <- as.list(mget(ls(pattern="^Matched_pDATA")))
  #### Generate a list of variables that we can use in the loop to run the chi.sq function
  variable_list <- list()
  
  variable_list <- append(variable_list, names(matched.test.pData))
  #### Exlude the nmb column from the list
  #variable_list2 <- variable_list[!variable_list =="NMB"]
  #### Can create a list of data to exclude form the chi.sq test but will need to ammend loop so easier to run all
  #exclude <- c("NMB", "Relapsetodeath", "Event", "EFS", "Followup", "OS.cat", "age.cont", "age.cat.adult.16", "age.cat.adult.21")
  ##### remove the excluded data from the variable list 
  #variable_list2 <- variable_list[!variable_list %in% exclude]
  cat (paste("processing Chi squared test for each variable ",gene_list[[gene]]), sep ="\n")
  #### itterate through the variables list generated and complete a chi squared test for each item in the list 
  for(v in 7:length(variable_list)){
    #### R was generating an error message that was stopping the program, but the objects were being created
    tryCatch({ #### This catches the error and outputs it to screen but allows the program to continue running
      #pdf(paste("Heatmap of chi.sq ",variable_list[v],".pdf"))  #### This outputs a pdf heatmap file per variable in the list
      #### assign the results of the chi squared test to res, [, v] defines a column number 1:
      #res <- chi.sq(matched.test.pData[, v], matched.goi.vsd.cat)  ### works for a single gene
      #### for the gene in the list of categories because the gene expression levels are different per gene a separate category is needed
      res <- chi.sq(na.omit(MATCHED_TEST_PDATA[[gene]][, v]), CAT_SEQ_LIST[[gene]][!is.na(MATCHED_TEST_PDATA[[gene]][, v])]) 
      #### now rename res according to the name of the item in the variable list being assessed
      assign(paste0("chi_res_",variable_list[[v]]), res)
    }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})  #### prints the error message to screen
  }
  #### create a list containing all of the results of the chi.squared test
  chi.sq.results <- as.list(mget(ls(pattern="^chi_res_")))
  #### determine the results in the list that are significant 
  significant_chi_results <- list()
  table_chi_results <- data.frame()
  #### create a matrix of defined size to hold the chisquare test results
  m <- matrix(0, ncol = 35, nrow = 1)
  #### Convert matrix m to a data frame 
  chisqout <- data.frame(m)
  #### create a list of names for results
  chsqnm <- names(chi.sq.results) 
  #### create an empty data frame to hold the results from the significance test 
  Chi_sq_sig_df <- data.frame(m)
  #### set the column names for the dataframe to match those in the test 
  colnames(Chi_sq_sig_df) <- chsqnm
  #### set the column names of the data frame to the names of the results
  #colnames(Chi_squared_sig_results_df) <- chsqnm
  for (c in 1:length(chsqnm)){
    #### extract the pvalue for all of the chi squared tests whether they are significant not
    #sig[[c]] <- chi.sq.results[[c]][[1]]$p.value
    sig <- chi.sq.results[[c]][[1]]$p.value
    tbl <- chi.sq.results[[c]][[2]]
    tbl_df <- as.data.frame(tbl)
    #### add the result from each chi square test to a list 
    significant_chi_results <- append(significant_chi_results, sig)
    table_chi_results <- rbind(table_chi_results, tbl_df)
  }
 
  #### turn the list into a data frame 
  chisqout <- as.data.frame(do.call(cbind, significant_chi_results))
  #### set the appropriate column names 
  colnames(chisqout) <- chsqnm
  #### remove named NA values from the data frame (created from the earlier matrix)
  chisqout <- chisqout[!is.na(names(chisqout))]
  #### assign the gene name to rownames so we know which gene was tested
  rownames(chisqout) <- gene_list[[gene]]
  #### create a large data frame with the results from each gene
  Chi_squared_results_df <- rbind(Chi_squared_results_df, chisqout)
  ### add additional script here to define other outputs including list (proportions) for chi squared
  
  #Chi_squared_results_df <- rbind(Chi_squared_results_df, chisqout)
  cat ("Ordering chi squared tests", sep ="\n")
  child_significant <- Chi_squared_results_df[order(Chi_squared_results_df$chi_res_childfac),]
  #child_significant <- order(Chi_squared_results_df$chi_res_childfac, decreasing=TRUE)
  agegrpfac_significant <- Chi_squared_results_df[order(Chi_squared_results_df$chi_res_agegrpfac),]
  
  cat (paste("Extracting significant results for children from chi squared tests ",gene_list[[gene]]), sep ="\n")
  chsq_significant_in_children <- list()
  for (rw in 1:length(Chi_squared_results_df$chi_res_childfac)){
    if(Chi_squared_results_df$chi_res_childfac[[rw]] < 0.05){
      sig <- Chi_squared_results_df$chi_res_childfac[[rw]]
      names(sig) <- rownames(Chi_squared_results_df[rw, ])
      chsq_significant_in_children <- append(chsq_significant_in_children, sig)
    }
  }
  cat (paste("Extracting significant results for sex from chi squared tests",gene_list[[gene]]), sep ="\n")
  chsq_significant_in_sexes <- list()
  for (rw in 1:length(Chi_squared_results_df$chi_res_sexfac)){
    if(Chi_squared_results_df$chi_res_sexfac[[rw]] < 0.05){
      sig <- Chi_squared_results_df$chi_res_sexfac[[rw]]
      names(sig) <- rownames(Chi_squared_results_df[rw, ])
      chsq_significant_in_sexes <- append(chsq_significant_in_sexes, sig)
    }
  }
  
  ### additional text added 8/8/17
  ### relapse
  
  cat (paste("Extracting significant results for relapse from chi squared tests",gene_list[[gene]]), sep ="\n")
  chsq_significant_in_relapse <- list()
  for (rw in 1:length(Chi_squared_results_df$chi_res_relapse)){
    #print(rownames(Chi_squared_results_df[rw, ]))
    if(Chi_squared_results_df$chi_res_relapse[[rw]] < 0.05){
      sig <- Chi_squared_results_df$chi_res_relapse[[rw]]
      names(sig) <- rownames(Chi_squared_results_df[rw ,])
      chsq_significant_in_relapse <- append(chsq_significant_in_relapse, sig)
    }
  }
  
  ###
  cat (paste("Extracting significant results for 4 subgroup methylation from chi squared tests",gene_list[[gene]]), sep ="\n")
  chsq_significant_in_meth4 <- list()
  for (rw in 1:length(Chi_squared_results_df$chi_res_meth)){
    #print(rownames(Chi_squared_results_df[rw, ]))
    if(Chi_squared_results_df$chi_res_meth[[rw]] < 0.05){
      sig <- Chi_squared_results_df$chi_res_meth[[rw]]
      names(sig) <- rownames(Chi_squared_results_df[rw ,])
      chsq_significant_in_meth4 <- append(chsq_significant_in_meth4, sig)
    }
  }
  
  ###
  
  cat (paste("Extracting significant results for 7 subgroup methylation from chi squared tests",gene_list[[gene]]), sep ="\n")
  chsq_significant_in_meth7 <- list()
  for (rw in 1:length(Chi_squared_results_df$chi_res_meth7)){
    #print(rownames(Chi_squared_results_df[rw, ]))
    if(Chi_squared_results_df$chi_res_meth7[[rw]] < 0.05){
      sig <- Chi_squared_results_df$chi_res_meth7[[rw]]
      names(sig) <- rownames(Chi_squared_results_df[rw ,])
      chsq_significant_in_meth7 <- append(chsq_significant_in_meth7, sig)
    }
  }
  
  cat (paste("Extracting significant results for MYC from chi squared tests",gene_list[[gene]]), sep ="\n")
  chsq_significant_in_MYC <- list()
  for (rw in 1:length(Chi_squared_results_df$chi_res_MYC.cat)){
    #print(rownames(Chi_squared_results_df[rw, ]))
    if(Chi_squared_results_df$chi_res_MYC.cat[[rw]] < 0.05){
      sig <- Chi_squared_results_df$chi_res_MYC.cat[[rw]]
      names(sig) <- rownames(Chi_squared_results_df[rw ,])
      chsq_significant_in_MYC <- append(chsq_significant_in_MYC, sig)
    }
  }
  
  cat (paste("Extracting significant results for MYCN from chi squared tests",gene_list[[gene]]), sep ="\n")
  chsq_significant_in_MYCN <- list()
  for (rw in 1:length(Chi_squared_results_df$chi_res_MYCN.cat)){
    #print(rownames(Chi_squared_results_df[rw, ]))
    if(Chi_squared_results_df$chi_res_MYCN.cat[[rw]] < 0.05){
      sig <- Chi_squared_results_df$chi_res_MYCN.cat[[rw]]
      names(sig) <- rownames(Chi_squared_results_df[rw ,])
      chsq_significant_in_MYCN <- append(chsq_significant_in_MYCN, sig)
    }
  }
  
  ### MYCMYCN
  
  cat (paste("Extracting significant results for MYCMYCN from chi squared tests",gene_list[[gene]]), sep ="\n")
  chsq_significant_in_MYCMYCN <- list()
  for (rw in 1:length(Chi_squared_results_df$chi_res_MYCMYCN.cat)){
    #print(rownames(Chi_squared_results_df[rw, ]))
    if(Chi_squared_results_df$chi_res_MYCMYCN.cat[[rw]] < 0.05){
      sig <- Chi_squared_results_df$chi_res_MYCMYCN.cat[[rw]]
      names(sig) <- rownames(Chi_squared_results_df[rw ,])
      chsq_significant_in_MYCMYCN <- append(chsq_significant_in_MYCMYCN, sig)
    }
  }
  
  
  ### LCA
  
  cat (paste("Extracting significant results for LCA from chi squared tests",gene_list[[gene]]), sep ="\n")
  chsq_significant_in_LCA <- list()
  for (rw in 1:length(Chi_squared_results_df$chi_res_MYCMYCN.cat)){
    #print(rownames(Chi_squared_results_df[rw, ]))
    if(Chi_squared_results_df$chi_res_LCA[[rw]] < 0.05){
      sig <- Chi_squared_results_df$chi_res_LCA[[rw]]
      names(sig) <- rownames(Chi_squared_results_df[rw ,])
      chsq_significant_in_LCA <- append(chsq_significant_in_LCA, sig)
    }
  }
  
  ### mstatus

  cat (paste("Extracting significant results for metastatic status from chi squared tests",gene_list[[gene]]), sep ="\n")
  chsq_significant_in_mstatus <- list()
  for (rw in 1:length(Chi_squared_results_df$chi_res_mstatus)){
    #print(rownames(Chi_squared_results_df[rw, ]))
    if(Chi_squared_results_df$chi_res_mstatus[[rw]] < 0.05){
      sig <- Chi_squared_results_df$chi_res_mstatus[[rw]]
      names(sig) <- rownames(Chi_squared_results_df[rw ,])
      chsq_significant_in_mstatus <- append(chsq_significant_in_mstatus, sig)
    }
  }
  
  ### resection
  
  
  #cat (paste("Extracting significant results for resection status from chi squared tests",gene_list[[gene]]), sep ="\n")
  #chsq_significant_in_resection <- list()
  #for (rw in 1:length(Chi_squared_results_df$chi_res_resection)){
    #print(rownames(Chi_squared_results_df[rw, ]))
  #  if(Chi_squared_results_df$chi_res_resection[[rw]] < 0.05){
  #    sig <- Chi_squared_results_df$chi_res_resection[[rw]]
  #    names(sig) <- rownames(Chi_squared_results_df[rw ,])
  #    chsq_significant_in_resection <- append(chsq_significant_in_resection, sig)
   # }
  #}
  
  
  
  ### TP53
  
  cat (paste("Extracting significant results for TP53 status from chi squared tests",gene_list[[gene]]), sep ="\n")
  chsq_significant_in_TP53 <- list()
  for (rw in 1:length(Chi_squared_results_df$chi_res_TP53.cat)){
    #print(rownames(Chi_squared_results_df[rw, ]))
    if(Chi_squared_results_df$chi_res_TP53.cat[[rw]] < 0.05){
      sig <- Chi_squared_results_df$chi_res_TP53.cat[[rw]]
      names(sig) <- rownames(Chi_squared_results_df[rw ,])
      chsq_significant_in_TP53 <- append(chsq_significant_in_TP53, sig)
    }
  }
  
  ### TERT
  cat (paste("Extracting significant results for TP53 status from chi squared tests",gene_list[[gene]]), sep ="\n")
  chsq_significant_in_TERT <- list()
  for (rw in 1:length(Chi_squared_results_df$chi_res_TERT.cat)){
    #print(rownames(Chi_squared_results_df[rw, ]))
    if(Chi_squared_results_df$chi_res_TERT.cat[[rw]] < 0.05){
      sig <- Chi_squared_results_df$chi_res_TERT.cat[[rw]]
      names(sig) <- rownames(Chi_squared_results_df[rw ,])
      chsq_significant_in_TERT <- append(chsq_significant_in_TERT, sig)
    }
  }
  
  
  
  ###
  cat (paste("Extracting significant results for q13 loss from chi squared tests",gene_list[[gene]]), sep ="\n")
  chsq_significant_in_q13loss <- list()
  for (rw in 1:length(Chi_squared_results_df$chi_res_q13loss)){
    #print(rownames(Chi_squared_results_df[rw, ]))
    if(Chi_squared_results_df$chi_res_q13loss[[rw]] < 0.05){
      sig <- Chi_squared_results_df$chi_res_q13loss[[rw]]
      names(sig) <- rownames(Chi_squared_results_df[rw ,])
      chsq_significant_in_q13loss <- append(chsq_significant_in_q13loss, sig)
    }
  }
  
  ### is this biomarker overrepresented in group that received RTX or CSI, or those classified as curative
  
  cat (paste("Extracting significant results for RTX from chi squared tests",gene_list[[gene]]), sep ="\n")
  chsq_significant_in_RTX <- list()
  for (rw in 1:length(Chi_squared_results_df$chi_res_RTX)){
    #print(rownames(Chi_squared_results_df[rw, ]))
    if(Chi_squared_results_df$chi_res_RTX[[rw]] < 0.05){
      sig <- Chi_squared_results_df$chi_res_RTX[[rw]]
      names(sig) <- rownames(Chi_squared_results_df[rw ,])
      chsq_significant_in_RTX <- append(chsq_significant_in_RTX, sig)
    }
  }
  
  
  cat (paste("Extracting significant results for CSI from chi squared tests",gene_list[[gene]]), sep ="\n")
  chsq_significant_in_CSI <- list()
  for (rw in 1:length(Chi_squared_results_df$chi_res_CSI)){
    #print(rownames(Chi_squared_results_df[rw, ]))
    if(Chi_squared_results_df$chi_res_CSI[[rw]] < 0.05){
      sig <- Chi_squared_results_df$chi_res_CSI[[rw]]
      names(sig) <- rownames(Chi_squared_results_df[rw ,])
      chsq_significant_in_CSI <- append(chsq_significant_in_CSI, sig)
    }
  }
  
  
  cat (paste("Extracting significant results for curative from chi squared tests",gene_list[[gene]]), sep ="\n")
  chsq_significant_in_curative <- list()
  for (rw in 1:length(Chi_squared_results_df$chi_res_curative)){
    #print(rownames(Chi_squared_results_df[rw, ]))
    if(Chi_squared_results_df$chi_res_curative[[rw]] < 0.05){
      sig <- Chi_squared_results_df$chi_res_curative[[rw]]
      names(sig) <- rownames(Chi_squared_results_df[rw ,])
      chsq_significant_in_curative <- append(chsq_significant_in_curative, sig)
    }
  }
  
  ###creating a list that outputs all the chisquare results
  
  ###creating a list that outputs all the significant chisquare results
  significant_chisquare_results <- as.list(mget(ls(pattern="chsq_significant_in_")))
  #### hashed out as unused, can be unhashed if needed
  #### Because we have a list of significant results each item within contains the table data, we can now run fishers test on each table in the list 
  #### First create a list to hold the results of the analyses 
 #fishers_results <- list()
  #cat (paste("Extracting fishers test results for histopathology ",gene_list[[gene]]), sep ="\n")
 # for (s in 1:nrow(Chi_squared_results_df)){
 # table <- table(as.factor(MATCHED_TEST_PDATA[[gene]]$histopath), as.factor(CAT_SEQ_LIST[[gene]]))
 # result <- fisher.test(table)
 # assign(paste0("Fishers_test_histopath ",gene_list[[gene]]), result)
 # }
 # fishers_results <- as.list(mget(ls(pattern="Fishers_test_histopath")))

  
  cat (paste("Generating Boxplots ",gene_list[[gene]]), sep ="\n")
  histopath.boxplot <- boxplot(SEQ_GENE_MATCHED[[gene]] ~ MATCHED_TEST_PDATA[[gene]]$histopath, col = c("red", "blue"), xlab = "Histopathology subtype", ylab = "Biomarker expression", main = "Correlation between biomarker and histopathology")
  age.boxplot <- boxplot(SEQ_GENE_MATCHED[[gene]] ~ MATCHED_TEST_PDATA[[gene]]$age.cat.infant, col = c("red", "blue"), xlab = "Infant", ylab = "Biomarker expression", main = "Correlation between biomarker and age (infant vs non infant)")
  sex.boxplot <- boxplot (SEQ_GENE_MATCHED[[gene]] ~ MATCHED_TEST_PDATA[[gene]]$sex, col = c("red", "blue"), xlab = "Gender", ylab = "Expression of biomarker", main = "Biomarker expression and gender")
  mstatus.boxplot <- boxplot(SEQ_GENE_MATCHED[[gene]] ~ MATCHED_TEST_PDATA[[gene]]$mstatus, col = c("red", "blue"), xlab = "M status", ylab = "Biomarker expression", main = "Correlation between biomarker and metastatic status")
  relapse.boxplot <- boxplot(SEQ_GENE_MATCHED[[gene]] ~ MATCHED_TEST_PDATA[[gene]]$relapse, col = c("red", "blue"), xlab = "Relapse status", ylab = "Biomarker expression",  main = "Correlation between biomarker and relapse")
  resection.boxplot <- boxplot(SEQ_GENE_MATCHED[[gene]] ~  MATCHED_TEST_PDATA[[gene]]$resection, col = c("red", "blue"), xlab = "Resection status", ylab = "Biomarker expression", main = "Correlation between biomarker and resection status")
  
  
  #### This stoped working when the output of chi squared tests was changed, new loop above only for histopath
  #### Now itterate through each of the items in the significant  results list 
  #for (s in 1:length(significant_chi_results)){
  #for (s in 1:nrow(Chi_squared_results_df)){
  #### create the individual tables
  #tab <- Chi_squared_results_df
  #### run the fishers test on the extracted table
  #result <- fisher.test(tab)
  #### generate a list of fishers test results 
  #fishers_results <- append(fishers_results, result)
  #}
  #### Errors will come up here if the biomarker of interest is not found in the RNAseq data
  ####  This will result in the SEQ_GENE_MATCHED[[gene]] matrix returning as empty (NAs) causing the program to crash
  cat (paste("processing Correlation coefficients for age variable ",gene_list[[gene]]), sep ="\n")
  ### Correlation coefficients
  #### run the tests you want to run on the continuous age data 
  x <- matched.test.pData$age.cont
  #### for a single gene
  #y <- matched.goi.vsd
  #### for a list of genes
  y <- SEQ_GENE_MATCHED[[gene]]
  age_results.cor <- cor.result(x,y)
  assign(paste0("age_correlation_",gene_list[[gene]]), age_results.cor)
  age_results.lin.reg <- lin.reg(x,y)
  assign(paste0("age_linear_regression_",gene_list[[gene]]), age_results.lin.reg)
  age_results.wilcox <- wilcox.test(x,y)
  assign(paste0("age_wilcox_HR_",gene_list[[gene]]), age_results.wilcox)
  
  age_correlation_results <- as.list(mget(ls(pattern="age_correlation_")))
  age_linear_regression_results <- as.list(mget(ls(pattern="age_linear_regression_")))
  age_wilcox_results <- as.list(mget(ls(pattern="age_wilcox_HR_")))
  
  ##################################
  ### logistic regression
  
  cat (paste("processing logistic regression for each variable ",gene_list[[gene]]), sep ="\n")
  #### create a vector of test factors
  test_factors <- as.vector(names(matched.test.pData))
  #### generate an empty list to hold the factors 
  fac <- list()
  #### itterate through the list of test factorss and generate the heading names add them into the fac list
  for(i in 1:length(test_factors)){
    fac <- append(fac, as.name(paste0("matched.test.pData$",test_factors[[i]], sep="")))
  }
  fac2 <- fac[c(8,13,26,14,19,20,22,23,24,25,34)]
  #### itterate through each of the factors in the list s and run a lofgistic regression on each 
  for(t in 3:length(fac2)){
    #### For a single gene 
    #regression <- logisticRegression(matched.test.pData[, t], matched.goi.vsd, matched.test.pData)
    #### For a list of genes 
    regression <- logisticRegression(matched.test.pData[, t], SEQ_GENE_MATCHED[[gene]], matched.test.pData)
    #### rename the output variable from each run through the loop so that the factor tested is attached to the object 
    assign(paste0("log_reg_",test_factors[[t]]), regression)
  }
  
  
  #### errors produced for logistic regression solved
  #### warnings() exlcluded using options (prevents the program crashing due to warning messages)
  #### create a list of all the logistic regression outputs based on the prefix assigned to the objects
  reg.log.list <- as.list(mget(ls(pattern="log_reg_")))
  #print(reg.log.list)
  cat (paste("processing pairwise t test for each variable and creating a list of results ",gene_list[[gene]]), sep ="\n")
  for(t in 3:length(fac)){
    #### for a single gene
    #pairwise <- pairwise.t.test(matched.goi.vsd, matched.test.pData[, t])
    #### for a list of genes 
    pairwise <- pairwise.t.test(SEQ_GENE_MATCHED[[gene]], matched.test.pData[, t])
    #### renmae the output variable from each run through the loop so that the factor tested is attached to the object 
    assign(paste0("pairwise_t_test",test_factors[[t]]), regression)
  }
  
  #### Create a list pairwise t test p value results
  pairwise_t_tests <- as.list(mget(ls(pattern="pairwise_t_test")))
  
  cat (paste("creating combined dataframe to assess biomarker in G3 G4 combined group, for survival cohort, aged 3-16 years, curative intent ",gene_list[[gene]]), sep = "\n")
  Curative_Treated_all_groups <- test.pData[which(test.pData$curative == "curative" & test.pData$childfac == "Child.M" | test.pData$childfac == "Child.F"),]
  
  Curative_Treated_G3G4 <- test.pData[which(test.pData$curative == "curative" & test.pData$childfac == "Child.M" | test.pData$childfac == "Child.F" & test.pData$subgroup4fac == "G3" | test.pData$subgroup4fac == "G4"),]

  
  #### Creating matched data frames containing RNAseq expression data and curative data for samples in test.pData
  treat_grps <- as.list(mget(ls(pattern="Curative_Treated_")))
  #names_treat_grps <- names(treat_grps)
  names_treat_grps <- c("all_groups", "G3_G4")
  cat (paste("creating curative treatment groups ",gene_list[[gene]]), sep = "\n")
  for (i in 1:length(treat_grps)){
    #### for a single gene
    #index.incl <- match(names(Seq_goi.vsd), rownames(treat_grps[[i]]))
    #### for a list of genes
    index.incl <- match(names(GENE_SEQ_LIST[[gene]]), rownames(treat_grps[[i]]))
    matched.test.incl.pData <- treat_grps[[i]][index.incl[!is.na(index.incl)],]
    assign(paste0("Matched_Curative_",names_treat_grps[[i]]),matched.test.incl.pData)
    #### for a single gene
    #matched.goi.vsd.incl <- Seq_goi.vsd[!is.na(index.incl)] 
    #### for a list of genes
    matched.goi.vsd.incl <- GENE_SEQ_LIST[[gene]][!is.na(index.incl)] 
    assign(paste0("Matched_GOI_Curative_",names_treat_grps[i]),matched.goi.vsd.incl)
    matched.goi.vsd.cat.incl <- ifelse(matched.goi.vsd.incl>median(GENE_SEQ_LIST[[gene]], na.rm = T), "high","low")
    assign(paste0("Matched_GOI_category_Curative_",names_treat_grps[i]),matched.goi.vsd.cat.incl)
  }
  
  #### create the lists for downstream processing 
  cat (paste("creating the lists for downstream processing ",gene_list[[gene]]), sep = "\n")
  curatives <- as.list(mget(ls(pattern="^Matched_Curative_")))
  names_curatives <- names(curatives)
  genesofinterest <- as.list(mget(ls(pattern="Matched_GOI_category_Curative")))
  
  #### change the name list to use otherwise re-run will cause the list to be regenrated containing binary data 
  cat (paste("creating binary relapse variables labelled 0,1 for event analysis",gene_list[[gene]]), sep = "\n")
  ### creating binary relapse variables labelled 0,1 for event analysis
  for (c in 1:length(curatives)){
    relapse_binary <- ifelse(curatives[[c]]$relapse == "relapse", 1, 0)
    assign(paste0("Binary_data_relapse_",names_curatives[[c]]),relapse_binary)
    Overall_survival_binary <- ifelse(curatives[[c]]$OS.cat == "Dead", 1, 0)
    assign(paste0("Binary_data_OS_",names_curatives[[c]]), Overall_survival_binary)
    EventFreeSurvival_binary <- ifelse(curatives[[c]]$Event == "Event", 1, 0)
    assign(paste0("Binary_data_EFS_",names_curatives[[c]]), EventFreeSurvival_binary)
  }
  
  cat (paste("creating binary lists for downstream processing ",gene_list[[gene]]), sep = "\n")
  #### create the lists for downstream processing 
  relapse_binaries <- as.list(mget(ls(pattern="Binary_data_relapse")))
  OS_binaries <- as.list(mget(ls(pattern="Binary_data_OS_")))
  EFS_binaries <- as.list(mget(ls(pattern="Binary_data_EFS_")))
  named_relapse_binaries <- names(relapse_binaries)
  named_OS_binaries <- names(OS_binaries)
  named_EFS_binaries <- names(EFS_binaries)
  groups <- c("all_groups", "G3_G4")
  binaries <- as.list(mget(ls(pattern="Binary_data_")))
  named_binaries <- names(binaries)
  
  #### Run Kaplan-Meier estimates for curative data using the binary survival data for the gene of interest 
  #### Plot the survival curves 
  #### Progression free survival
  cat (paste("Running Kaplan-Meier estimates for curative data for Progression free survival ",gene_list[[gene]]), sep = "\n")
  for (cur in 1:length(curatives)){
    km_results_log.test<- km.log.test(time = curatives[[cur]]$PFS, event = relapse_binaries[[cur]], marker = genesofinterest[[cur]])
    assign(paste0("Kaplan_Meier_",named_relapse_binaries[[cur]]), km_results_log.test)
    km.PFS.incl <- survfit(Surv(curatives[[cur]]$PFS, relapse_binaries[[cur]])~genesofinterest[[cur]], type = "kaplan-meier", conf.type = "log")
    #km.PFS.incl <- survfit(Surv(curatives[[cur]]$PFS, binaries[[cur]])~genesofinterest[[cur]], type = "kaplan-meier", conf.type = "log")
    plot(km.PFS.incl, yaxt="n", col = c("red", "blue"), xlab = "time to progression/relapse (years)", ylab = "PFS (%)", xlim = c(0,10), main = paste("Expression of ", gene_list[[gene]]," and progression-free survival (PFS)",groups[[cur]]), lty = 1:2)
    PFS.names <- c("biomarker - high", "biomarker - low")
    legend (x="topright", PFS.names,  lty= 1:2, col = c("red","blue"))
    axis(2, at=pretty(relapse_binaries[[cur]]), lab=pretty(binaries[[cur]]) * 100, las=TRUE)
    PFS.incl.logrank <- survdiff(Surv(curatives[[cur]]$PFS, binaries[[cur]]) ~ genesofinterest[[cur]])
    1 - pchisq(PFS.incl.logrank$chisq, length(PFS.incl.logrank$obs)-1) -> surv.p.val
    text(4,0.1,paste("p =",round(surv.p.val, 3)), pos = 4, cex = 1)
    assign(paste0("Survival_Fit_",named_relapse_binaries[[cur]],gene_list[[gene]]), km.PFS.incl)
    assign(paste0("Survival_pval_",named_relapse_binaries[[cur]],gene_list[[gene]]), surv.p.val)
  }
  cat (paste("Running Kaplan-Meier estimates for curative data for Overall survival ",gene_list[[gene]]), sep = "\n")
  #### Overall survival
  for (cur in 1:length(curatives)){
    km_results_log.test<- km.log.test(time = curatives[[cur]]$Followup, event = OS_binaries[[cur]], marker = genesofinterest[[cur]])
    assign(paste0("Kaplan_Meier_",named_OS_binaries[[cur]]), km_results_log.test)
    #km.PFS.incl <- survfit(Surv(curatives[[cur]]$Followup, binaries[[cur]])~genesofinterest[[cur]], type = "kaplan-meier", conf.type = "log")
    km.PFS.incl <- survfit(Surv(curatives[[cur]]$Followup, OS_binaries[[cur]])~genesofinterest[[cur]], type = "kaplan-meier", conf.type = "log")
    plot(km.PFS.incl, yaxt="n", col = c("red", "blue"), xlab = "overall survival (years)", ylab = "OS (%)", xlim = c(0,10), main = paste("Expression of ",gene_list[[gene]]," and overall survival (OS)",groups[[cur]]), lty = 1:2)
    PFS.names <- c("biomarker - high", "biomarker - low")
    legend (x="topright", PFS.names,  lty= 1:2, col = c("red","blue"))
    axis(2, at=pretty(binaries[[cur]]), lab=pretty(binaries[[cur]]) * 100, las=TRUE)
    PFS.incl.logrank <- survdiff(Surv(curatives[[cur]]$Followup, binaries[[cur]]) ~ genesofinterest[[cur]])
    1 - pchisq(PFS.incl.logrank$chisq, length(PFS.incl.logrank$obs)-1) -> surv.p.val
    text(4,0.1,paste("p =",round(surv.p.val, 3)), pos = 4, cex = 1)
    assign(paste0("Survival_Fit_",named_OS_binaries[[cur]],gene_list[[gene]]), km.PFS.incl)
    assign(paste0("Survival_pval_",named_OS_binaries[[cur]],gene_list[[gene]]), surv.p.val)
  }
  cat (paste("Running Kaplan-Meier estimates for curative data for Event free survival ",gene_list[[gene]]), sep = "\n")
  #### Event free survival
  for (cur in 1:length(curatives)){
    km_results_log.test<- km.log.test(time = curatives[[cur]]$EFS, event = EFS_binaries[[cur]], marker = genesofinterest[[cur]])
    assign(paste0("Kaplan_Meier_",named_EFS_binaries[[cur]]), km_results_log.test)
    #km.PFS.incl <- survfit(Surv(curatives[[cur]]$EFS, binaries[[cur]])~ genesofinterest[[cur]], type = "kaplan-meier", conf.type = "log")
    km.PFS.incl <- survfit(Surv(curatives[[cur]]$EFS, EFS_binaries[[cur]])~ genesofinterest[[cur]], type = "kaplan-meier", conf.type = "log")
    #pdf(file=paste("marker_results_EFS",groups[[cur]],".pdf"))
    plot(km.PFS.incl, yaxt="n", col = c("red", "blue"), xlab = "event-free survival (years)", ylab = "EFS (%)", xlim = c(0,10), main = paste("Expression of ",gene_list[[gene]]," and event-free survival (EFS)",groups[[cur]]), lty = 1:2)
    PFS.names <- c("biomarker - high", "biomarker - low")
    legend (x="topright", PFS.names,  lty= 1:2, col = c("red","blue"))
    axis(2, at=pretty(binaries[[cur]]), lab=pretty(binaries[[cur]]) * 100, las=TRUE)
    PFS.incl.logrank <- survdiff(Surv(curatives[[cur]]$EFS, binaries[[cur]]) ~ genesofinterest[[cur]])
    1 - pchisq(PFS.incl.logrank$chisq, length(PFS.incl.logrank$obs)-1) -> surv.p.val
    text(4,0.1,paste("p =",round(surv.p.val, 3)), pos = 4, cex = 1)
    assign(paste0("Survival_Fit_",named_EFS_binaries[[cur]],gene_list[[gene]]), km.PFS.incl)
    assign(paste0("Survival_pval_",named_EFS_binaries[[cur]],gene_list[[gene]]), surv.p.val)
  }
  
  cat (paste("extracting the 5 year survival statistics for the biomarker ",gene_list[[gene]]), sep = "\n")
  #### extract the 5 year survival statistics for the biomarker
  summaries_survival_fits <- as.list(mget(ls(pattern="Survival_Fit_")))
  naming_survival_fits <- paste(names(summaries_survival_fits),gene_list[[gene]])
  for (sum in 1:length(summaries_survival_fits)){
    res <- summary(summaries_survival_fits[[sum]])
    cols <- lapply(c(2:6, 8:11) , function(x) res[x])
    tbl <- do.call(data.frame, cols)
    #write.table(tbl, file=paste("kaplan meier ",naming_survival_fits[[sum]]), sep="\t", quote=FALSE)
    greater_than_5_yrs <- tbl[which(tbl$time > 5),]
    assign(paste0("greater_than_5_years_",gene_list[[gene]]),greater_than_5_yrs)
  }
  cat (paste("Running Cox haxards ratio test for event free survival ",gene_list[[gene]]), sep = "\n")
  #### Cox haxards ratio test for event free survival
  for (cur in 1:length(curatives)){
    cox.relapse.incl <- coxph (Surv(curatives[[cur]]$EFS, EFS_binaries[[cur]]) ~ genesofinterest[[cur]])
    assign(paste0("cox_relapse_",named_EFS_binaries[[cur]]), cox.relapse.incl)
  }
  #### Cox haxards ratio test for overlall survival
  for (cur in 1:length(curatives)){
    cox.relapse.incl <- coxph (Surv(curatives[[cur]]$Followup, OS_binaries[[cur]]) ~ genesofinterest[[cur]])
    assign(paste0("cox_relapse_",named_OS_binaries[[cur]]), cox.relapse.incl)
  }
  #### Cox haxards ratio test for progression free survival
  for (cur in 1:length(curatives)){
    cox.relapse.incl <- coxph (Surv(curatives[[cur]]$PFS, relapse_binaries[[cur]]) ~ genesofinterest[[cur]])
    assign(paste0("cox_relapse_",named_relapse_binaries[[cur]]), cox.relapse.incl)
  }
  cat (paste("generating cox relapse list summaries for downstream filtering",gene_list[[gene]]), sep = "\n")
  Output_Cox <- as.list(mget(ls(pattern="cox_relapse_")))
  for (cox in 1:length(Output_Cox)){
    summary(Output_Cox[[cox]])$logtest
    print(summary(Output_Cox[[cox]]))
  }
  cat (paste("Create data frames to hold the desired output results ",gene_list[[gene]]), sep = "\n")
  #### Create a data frame to hold all of the desired output results 
  cat (paste("Cox hazards data frame being compiled ",gene_list[[gene]]), sep = "\n")
  cox_names <- toupper(names(Output_Cox))
  
  ##### Extraction of relevant values 
  for (cox in 1:length(Output_Cox)){
    cat (paste("Extracting relevant Cox hazard test values",gene_list[[gene]]), sep = "\n")
    cox.pvalue <- summary(Output_Cox[[cox]])$coefficients[5]
    cox.HazR <- summary(Output_Cox[[cox]])$coefficients[2]
    cox.n <- summary(Output_Cox[[cox]])$n
    cox.nevent <- summary(Output_Cox[[cox]])$nevent
    cox.UCI.95 <- summary(Output_Cox[[cox]])$conf.int[4]
    cox.LCI.95 <- summary(Output_Cox[[cox]])$conf.int[3]
  }
  COX_OUT_DF <- data.frame(cox.pvalue, cox.HazR, cox.n, cox.nevent, cox.UCI.95, cox.LCI.95)
  rownames(COX_OUT_DF) <- gene_list[[gene]]
  assign(paste0("COX_data_frame_",gene_list[[gene]], cox_names[[cox]]), COX_OUT_DF)
  ALL_COX_Results <- as.list(mget(ls(pattern="^COX_data_frame")))
  cat (paste("Extracting significant COX results ",gene_list[[gene]]), sep = "\n")
  for (cr in 1:length(ALL_COX_Results)){
    most_significant_Cox <- ALL_COX_Results[[cr]][which(ALL_COX_Results[[cr]]$cox.pvalue < 0.05)]
  }
  
  
  cat (paste("Creating Logistic Regression data frame ",gene_list[[gene]]), sep = "\n")
  m <- matrix()
  logistic_regression_out <- data.frame(m)
  names_reg_log <- toupper(names(reg.log.list))
  for (res in 1:length(reg.log.list)){
    lr.pval <- reg.log.list[[res]][2, 7]
    lr.OR <- reg.log.list[[res]][2, 1]
    lr.CI.97.5 <- reg.log.list[[res]][2, 3]
    lr.CI.2.5 <- reg.log.list[[res]][2, 2]
  }
  
  LR_DF <- data.frame(lr.pval, lr.OR, lr.CI.97.5, lr.CI.2.5)
  rownames(LR_DF) <- gene_list[[gene]]
  assign(paste0("LR_data_frame_", gene_list[[gene]],names_reg_log[[res]]), LR_DF)
  ALL_Log_Reg_Results <- as.list(mget(ls(pattern="^LR_data_frame")))
  cat (paste("logistic regression data frame created ",gene_list[[gene]]), sep = "\n")
  cat (paste("Extracting significant Logistic Regression results ",gene_list[[gene]]), sep = "\n")
  for (lr in 1:length(ALL_Log_Reg_Results)){
    most_significant_Log_Reg <- ALL_Log_Reg_Results[[lr]][which(ALL_Log_Reg_Results[[lr]]$lr.pvalue < 0.05)]
  }
}

significant_chisquare_results_df <- data.frame()
significant_chisquare_results_df = as.data.frame(do.call(rbind, significant_chisquare_results)) 
for (rw in 1:nrow(significant_chisquare_results_df)){
  significant_chisquare_results_df[, rw] <- as.numeric(significant_chisquare_results_df[, rw])
}
#sink()
#dev.off()
biomarkers_greater_than_5_years <- as.list(mget(ls(pattern = "greater_than_5_years_")))
biomarkers_greater_than_5_years_df <- data.frame(biomarkers_greater_than_5_years)
  #### write a table of significant results to output 
  #write.table(most_significant_Cox, file=paste("most_significant_Cox",names_sig_cox), sep="\t", quote=FALSE)
  write.table(most_significant_Log_Reg, file="most_significant_logistic_regression_results.txt", sep="\t", quote=FALSE, col.names=TRUE, row.names=TRUE)
  write.table(most_significant_Cox, file="most_significant_Cox_hazard_results.txt", sep="\t", quote=FALSE, col.names=TRUE, row.names=TRUE)               
  write.table(biomarkers_greater_than_5_years_df , file="most_significant_survival_results.txt", sep="\t", quote=FALSE, col.names=TRUE, row.names=TRUE)
  #write.table(significant_in_children, file="most_significant_Chi_squared_results_for_children.txt", sep="\t", quote=FALSE, col.names=TRUE, row.names=TRUE)
  write.table(significant_chisquare_results_df, file="most_significant_Chi_squared_results_for_all_test_categories.txt", sep="\t", quote=FALSE, col.names=TRUE, row.names=TRUE)
  
  ###
 ### want to have G3G4 output for the above as well
  #write(significant_chisquare_results_df, file="most_significant_Chi_squared_results_for_all_test_categories.txt", sep="\n")
  