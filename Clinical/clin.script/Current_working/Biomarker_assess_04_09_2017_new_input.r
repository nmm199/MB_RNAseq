#!/usr/bin/Rscript
library(NMF)
#source("http://bioconductor.org/biocLite.R")
#biocLite(c("gplots", "car"))
library('gplots')
library(car)
library(stats)
library('survival')
library(scales)
library(foreach)
library(doParallel)
registerDoParallel(10) ##### changed to 10 from 32

setwd("/home/nmm199/R/MB_RNAseq/") ### to change this depending on who is working on the script
load("pData_saved.RData")
options(warn=-1)
options(browserNLdisabled = TRUE)
cox_hazard_results_df <- data.frame()
survival_results_df <- data.frame()
Chi_squared_results_df <- data.frame()
logistic_regression_results_df <- data.frame()
Chi_squared_sig_results_df <- data.frame()

#### runs the first ten gene ids in the list
gene_list <- as.list(head(rownames(mb.vsd), 28))
#gene_list <- as.list(rownames(mb.vsd))
for (gene in 1:length(gene_list)){
#foreach (gene = 1:length(gene_list))%dopar%{
  #### unhash sink and pdf lines below to output a sink and pdf file per gene in the list
  #sink(paste("pDatalog",gene_list[[gene]],".txt"))
  #pdf(paste("marker.results",gene_list[[gene]],".pdf"), width = 10, height = 10)
  #### create an ordered list of variables measured in the matched data 
  index <- match(colnames(mb.vsd), rownames(test.pData))
  #### Generate a list of expression values for the gene of interest 
  Seq_goi.vsd <- as.numeric(mb.vsd[gene_list[[gene]],])
  names(Seq_goi.vsd) <- gsub("T","",names(mb.vsd))
  ############################################## SET THE GROUPS TO BE TESTED IN THE ANALYSIS ############################
  ############## matched.test.pData is all groups #################
  matched.test.pData <- test.pData[index[!is.na(index)],] 
  matched.goi.vsd <- Seq_goi.vsd[!is.na(index)] 
  matched.goi.vsd.cat <- ifelse(matched.goi.vsd>median(Seq_goi.vsd, na.rm = T), "high","low")
  #### to look at group 3 and 4 unhash below
  ############# This will run only G3G4, if this is hashed out all run, if this is unhased please hash out matched.test.pData
  #G3G4 <- matched.test.pData[which(matched.test.pData$group3or4fac =="G3G4"),]
  #### now assign these to names prefixed for grouping and append the gene identifier to each of the objects
  assign(paste0("Gene_Seq ",gene_list[[gene]]), Seq_goi.vsd)
  assign(paste0("Seq_gene_matched ",gene_list[[gene]]), matched.goi.vsd)
  assign(paste0("Cat_Seq_gene_matched ",gene_list[[gene]]), matched.goi.vsd.cat)
  assign(paste0("Matched_pDATA ",gene_list[[gene]]), matched.test.pData)
  #### to look at group 3 and 4 unhash below
  #assign(paste0("Matched_pDATA_G3G4 ",gene_list[[gene]]), G3G4)
  #### Use the prefixes generated above to place the objects in lists that were can itterate through
  CAT_SEQ_LIST <- as.list(mget(ls(pattern="^Cat_Seq_gene_matched")))
  SEQ_GENE_MATCHED <- as.list(mget(ls(pattern="^Seq_gene_matched")))
  GENE_SEQ_LIST <- as.list(mget(ls(pattern="^Gene_Seq")))
  MATCHED_TEST_PDATA <- as.list(mget(ls(pattern="^Matched_pDATA")))
  cat (paste("processing Chi squared test for each variable ",gene_list[[gene]]), sep ="\n")
  
  chi.sq.results <- foreach (v = 1:length(MATCHED_TEST_PDATA[[gene]]), .combine=append, .packages=c("data.table"))%dopar%{
    tryCatch({ #### This catches the error and outputs it to screen but allows the program to continue running
      res <- chi.sq(as.matrix(na.omit(MATCHED_TEST_PDATA[[gene]][, v])), as.matrix(t(CAT_SEQ_LIST[[gene]][!is.na(MATCHED_TEST_PDATA[[gene]][, v])])))  
      chi_results <- list(res[[1]], res[[2]], res[[3]], res[[4]])
      names(chi_results) <- paste0(colnames(MATCHED_TEST_PDATA[[gene]][v]), "_", gene_list[[gene]])
      return(chi_results)
    }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})  #### prints the error message to screen
  }
  
  
  all_chi.sq <- as.list(chi.sq.results)

  cat (paste("Extraction significant chi.sq.results ",gene_list[[gene]]), sep ="\n")
  #vts <- rapply(chi.sq.results, function(x) head(x, 1, x))
  #vts <- rapply(all_chi.sq, function(x) head(x, 1, x))
  vts <- rapply(all_chi.sq, function(x) head(x[1], 1, x))
  #chi.sq.p.values <- vts[seq(3, length(vts), 19)]
  #chi.sq.p.values <- vts[seq(3, length(vts), 12)]
  chi.sq.p.values <- subset(vts, grepl(glob2rx("*p.value"), names(vts)))
  #chi.sq.p.list <- as.list(chi.sq.p.values)
  chi.sq.p.list <- list()
  chi.sq.p.list <- append(chi.sq.p.list, chi.sq.p.values)
  chisquare_results_df = as.data.frame(do.call(rbind, chi.sq.p.list)) 
  chisquare_results_df$V1 <- as.numeric(as.character(chisquare_results_df$V1))
  colnames(chisquare_results_df) <- "P.value"
  Chi_squared_results_df <- rbind(Chi_squared_results_df, chisquare_results_df)
  
  #significant_chisq_res <- as.data.frame(chisquare_results_df[which(chisquare_results_df$P.value < 0.05), drop=FALSE,])
  significant_chisq_res <- as.data.frame(Chi_squared_results_df[which(Chi_squared_results_df$P.value < 0.05), drop=FALSE,])
  
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
  #x <- matched.test.pData$age.cont
  x <- MATCHED_TEST_PDATA[[gene]]$age.cont
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
  
  #logisticRegression <- function(x,y, data) {
  #  temp.glm <- glm (x ~ y, family = binomial (link = 'logit'), data= data)
  #  return(cbind(OR=exp(coef(temp.glm)), exp(confint(temp.glm)), summary(temp.glm)$coefficients))
  #}
  cat (paste("processing logistic regression for each variable ",gene_list[[gene]]), sep ="\n")
  #subset_df <- MATCHED_TEST_PDATA[[gene]][c(7,8,10,11,13,15,16,17,18,26,14,19,20,21,22,23,24,25,27,28,34,35)]  ## rm 12, 20
  subset_df <- MATCHED_TEST_PDATA[[gene]][c(7,8,10,11,13,15,16,17,18,26,14,19,20,21,22,23,25,27,28,34,35)]  ## rm 12, 20, 24
  #subset_df$t <- factor(subset_df$t)
  #subset_df <- as.factor(subset_df[, 20])
  ## rm 12
  #### create a vector of test factors
  test_factors <- as.vector(names(subset_df))
  #### itterate through each of the factors in the list s and run a lofgistic regression on each 
  for(t in 1:length(test_factors)){
    subset_df[, t] <- factor(subset_df[, t])
    #### For a single gene 
    #regression <- logisticRegression(matched.test.pData[, t], matched.goi.vsd, matched.test.pData)
    #### For a list of genes 
    #regression <- logisticRegression(na.omit(subset_df[, t], SEQ_GENE_MATCHED[[gene]][!is.na(subset_df[, t])], subset_df))
    regression <- logisticRegression(subset_df[, t], SEQ_GENE_MATCHED[[gene]], subset_df)
    #### rename the output variable from each run through the loop so that the factor tested is attached to the object 
    assign(paste0("log_reg_",test_factors[[t]]), regression)
  }
  
  #### errors produced for logistic regression solved
  #### warnings() exlcluded using options (prevents the program crashing due to warning messages)
  #### create a list of all the logistic regression outputs based on the prefix assigned to the objects
  reg.log.list <- as.list(mget(ls(pattern="log_reg_")))
  #print(reg.log.list)
  cat (paste("processing pairwise t test for each variable and creating a list of results ",gene_list[[gene]]), sep ="\n")
  #for(t in 3:length(fac)){
  for(t in 1:length(test_factors)){
    #### for a single gene
    #pairwise <- pairwise.t.test(matched.goi.vsd, matched.test.pData[, t])
    #### for a list of genes 
    pairwise <- pairwise.t.test(SEQ_GENE_MATCHED[[gene]], MATCHED_TEST_PDATA[[gene]][, t])
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
	curatives[[cur]]$PFS <- as.numeric(as.character(curatives[[cur]]$PFS))
    #curatives[[cur]]$PFS <- as.numeric(unique(levels(curatives[[cur]]$PFS)))
	#curatives[[cur]]$PFS <- as.numeric(curatives[[cur]]$PFS)
    km_results_log.test<- km.log.test(time = as.numeric(curatives[[cur]]$PFS), event = relapse_binaries[[cur]], marker = genesofinterest[[cur]])
    assign(paste0("Kaplan_Meier_",named_relapse_binaries[[cur]]), km_results_log.test)
    km.PFS.incl <- survfit(Surv(as.numeric(curatives[[cur]]$PFS), relapse_binaries[[cur]])~genesofinterest[[cur]], type = "kaplan-meier", conf.type = "log")
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
	curatives[[cur]]$OS <- as.numeric(as.character(curatives[[cur]]$Followup))
    #curatives[[cur]]$OS <- as.numeric(unique(levels(curatives[[cur]]$Followup)))
	#curatives[[cur]]$OS <- as.numeric(curatives[[cur]]$Followup)
    km_results_log.test<- km.log.test(time = as.numeric(curatives[[cur]]$Followup), event = OS_binaries[[cur]], marker = genesofinterest[[cur]])
    assign(paste0("Kaplan_Meier_",named_OS_binaries[[cur]]), km_results_log.test)
    #km.PFS.incl <- survfit(Surv(curatives[[cur]]$Followup, binaries[[cur]])~genesofinterest[[cur]], type = "kaplan-meier", conf.type = "log")
    km.PFS.incl <- survfit(Surv(as.numeric(curatives[[cur]]$Followup), OS_binaries[[cur]])~genesofinterest[[cur]], type = "kaplan-meier", conf.type = "log")
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
	curatives[[cur]]$EFS <- as.numeric(as.character(curatives[[cur]]$EFS))
    #curatives[[cur]]$EFS <- as.numeric(unique(levels(curatives[[cur]]$EFS)))
	#curatives[[cur]]$EFS <- as.numeric(curatives[[cur]]$EFS)
    km_results_log.test<- km.log.test(time = as.numeric(curatives[[cur]]$EFS), event = EFS_binaries[[cur]], marker = genesofinterest[[cur]])
    assign(paste0("Kaplan_Meier_",named_EFS_binaries[[cur]]), km_results_log.test)
    #km.PFS.incl <- survfit(Surv(curatives[[cur]]$EFS, binaries[[cur]])~ genesofinterest[[cur]], type = "kaplan-meier", conf.type = "log")
    km.PFS.incl <- survfit(Surv(as.numeric(curatives[[cur]]$EFS), EFS_binaries[[cur]])~ genesofinterest[[cur]], type = "kaplan-meier", conf.type = "log")
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
  
  cat (paste("extracting the significant survival statistics for the biomarker ",gene_list[[gene]]), sep = "\n")
  #### extract significant survival analysis results 
  survival_pvals <- as.list(mget(ls(pattern="Survival_pval_")))
  survival_pvals_df <- as.data.frame(do.call("rbind",survival_pvals)) 
  colnames(survival_pvals_df) <- ("p.value")
  significant_survival_pvals_df <- as.data.frame(survival_pvals_df[which(survival_pvals_df$p.value < 0.05), drop=FALSE,])
  
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
    COX_OUT_DF <- data.frame(cox.pvalue, cox.HazR, cox.n, cox.nevent, cox.UCI.95, cox.LCI.95)
    assign(paste0("COX_data_frame_",gene_list[[gene]], cox_names[[cox]]), COX_OUT_DF)
  }
  ALL_COX_Results <- as.list(mget(ls(pattern="^COX_data_frame")))
  COX_DF <- as.data.frame(do.call("rbind",ALL_COX_Results)) 
  cat (paste("Extracting significant COX results ",gene_list[[gene]]), sep = "\n")
  
  most_significant_Cox <- COX_DF[which(COX_DF$cox.pvalue < 0.05),]
  
  
  cat (paste("Creating Logistic Regression data frame ",gene_list[[gene]]), sep = "\n")
  m <- matrix()
  logistic_regression_out <- data.frame(m)
  names_reg_log <- toupper(names(reg.log.list))
  for (res in 1:length(reg.log.list)){
    lr.pval <- reg.log.list[[res]][2, 7]
    lr.OR <- reg.log.list[[res]][2, 1]
    lr.CI.97.5 <- reg.log.list[[res]][2, 3]
    lr.CI.2.5 <- reg.log.list[[res]][2, 2]
    LR_DF <- data.frame(lr.pval, lr.OR, lr.CI.97.5, lr.CI.2.5)
    assign(paste0("LR_DF_",names_reg_log[[res]],"_",gene_list[[gene]]), LR_DF)
  }
  LR_list <- as.list(mget(ls(pattern="LR_DF_")))
  DF_LR <- as.data.frame(do.call("rbind",LR_list)) 
  cat (paste("logistic regression data frame created ",gene_list[[gene]]), sep = "\n")
  cat (paste("Extracting significant Logistic Regression results ",gene_list[[gene]]), sep = "\n")
  most_significant_Log_Reg <- DF_LR[which(DF_LR$lr.pval < 0.05),]
  
  
  
  #sink()
  #dev.off()
  cat (paste("wrting significance tables ",gene_list[[gene]]), sep = "\n")
  write.table(most_significant_Log_Reg, file="most_significant_logistic_regression_results_all_genes.txt", sep="\t", quote=FALSE, col.names=TRUE, row.names=TRUE)
  write.table(most_significant_Cox, file="most_significant_Cox_hazard_results_all_genes.txt", sep="\t", quote=FALSE, col.names=TRUE, row.names=TRUE)               
  write.table(significant_chisq_res, file="most_significant_Chi_squared_results_for_all_test_categories_all_genes.txt", sep="\t", quote=FALSE, col.names=TRUE, row.names=TRUE)
}

  
  