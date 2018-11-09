
### Script for analysing the output from univariate analysis of RNA expression data compared to survival outcomes
### Date: September 25 2017
### Author: Dr Marion Mateos

### choose the file input

# results.master <- readRDS (file = "/home/nmm199/R/MB_RNAseq/Clinical/clin.results/master/results.filt.complete.20180220.rds")   ### samples filtered only

 results.master <- readRDS (file = "/home/nmm199/R/MB_RNAseq/Clinical/clin.results/master/Complete_transcripts/results.filt.genefilter.20181031.rds") ### samples and genes filtered

# results.master <- readRDS (file = "/home/nmm199/R/MB_RNAseq/Clinical/clin.results/master/Novel_transcripts/results.filt.genefilter.novel.20181031.rds")

# results.master <- readRDS (file = "/home/nmm199/R/MB_RNAseq/Clinical/clin.results/master/results.filt.genefilter.random.20180327.rds")


### read in functions file

source(file = "/home/nmm199/R/MB_RNAseq/Clinical/clin.script/master_script/clinical_data_functions_master.R")
source(file = "/home/nmm199/R/MB_RNAseq/Clinical/clin.script/master_script/clinical_data_functions_extract_master.R")

### file output 

### dataframes with p value, adjusted p value
### all.survival.p.bothgroups: includes OS, EFS and PFS for overall group and G3G4 
### significant adjusted p values for OS, EFS, PFS separately in dataframes for overall group and for G3G4
### graphical output for adjusted p values compared to p values
### univariate cox regression p values, Z scores, HR, 95 CI
### multivariate cox regression p values Z scores, HR, 95 CI
### working on logistic regression p values, OR, 95 CI
### working on chi squared p values, adjusted p values
### graphical output to include survival curves, relationship to subgroup


##########################################################################

### Extracted dataframes

### KM Survival overall
### see extraction functions for KM

######################
### OS p values
km.OS.all.results <- extract.km.OS.pval(results.master, name = "surv.km.OS.all")               
km.OS.G3G4.results <- extract.km.OS.pval (results.master, name = "surv.km.OS.G3G4")  
km.OS.G4.results <- extract.km.OS.pval (results.master, name = "surv.km.OS.G4")
km.OS.SHH.results <- extract.km.OS.pval(results.master, name = "surv.km.OS.SHH")             
km.OS.SHH.old.results <- extract.km.OS.pval(results.master, name = "surv.km.OS.SHH.old")            

### EFS p values
km.EFS.all.results <- extract.km.EFS.pval(results.master, name = "surv.km.EFS.all")  
km.EFS.G3G4.results <- extract.km.EFS.pval(results.master, name = "surv.km.EFS.G3G4") 

### PFS p values
km.PFS.all.results <- extract.km.PFS.pval(results.master, name = "surv.km.PFS.all")  
km.PFS.G3G4.results <- extract.km.PFS.pval(results.master, name = "surv.km.PFS.G3G4") 
km.PFS.G4.results <- extract.km.PFS.pval(results.master, name = "surv.km.PFS.G4") ### added 25/10/18 G4 specific
km.PFS.SHH.results <- extract.km.PFS.pval(results.master, name = "surv.km.PFS.SHH")  
km.PFS.SHH.old.results <- extract.km.PFS.pval(results.master, name = "surv.km.PFS.SHH.old")  
 

### significant dataframes for adjusted p values

significant.km.EFS.all <- km.EFS.all.results[which(km.EFS.all.results[, 2]<0.05),]
significant.km.EFS.G3G4 <- km.EFS.G3G4.results [which(km.EFS.G3G4.results[, 2]<0.05),]

significant.km.OS.all <- km.OS.all.results [which (km.OS.all.results[,2] <0.05),]
significant.km.OS.G3G4 <- km.OS.G3G4.results[which(km.OS.G3G4.results [, 2]<0.05), ]
significant.km.OS.G4 <- km.OS.G4.results[which(km.OS.G4.results [,2]<0.05), ] ### added 25/10/18
significant.km.OS.SHH <- km.OS.SHH.results [which(km.OS.SHH.results [, 2]<0.05), ]
significant.km.OS.SHH.old <- km.OS.SHH.old.results[which(km.OS.SHH.old.results[, 2]<0.05),]

significant.km.PFS.all <- km.PFS.all.results[which(km.PFS.all.results[, 2]<0.05),]              
significant.km.PFS.G3G4 <- km.PFS.G3G4.results[which(km.PFS.G3G4.results[, 2]<0.05),]           
significant.km.PFS.G4 <- km.PFS.G4.results[which(km.PFS.G4.results[,2]<0.05),]                                                                      
significant.km.PFS.SHH <- km.PFS.SHH.results[which(km.PFS.SHH.results[, 2]<0.05),]              
significant.km.PFS.SHH.old <- km.PFS.SHH.old.results[which(km.PFS.SHH.old.results[, 2]<0.05),]  


##########################################################################################################################
##########################################################################################################################

### cox PFS for continuous variable, overall category
### this section has been updated below , specify results.master then subset.index as a name (27/2/18)

cox.PFS.cat.all.df <- extract.cox (results.master, "surv.cox.relapse.incl.cat" ) ### was 3

cox.PFS.cont.all.df <- extract.cox(results.master, "surv.cox.relapse.incl.contin" ) ### was 4, and so forth for below

cox.PFS.cat.G3G4.df <- extract.cox(results.master, "surv.cox.relapse.incl.G3G4.cat")

cox.PFS.cont.G3G4.df <- extract.cox(results.master, "surv.cox.relapse.incl.G3G4.contin")

cox.PFS.cat.G4.df <- extract.cox(results.master, "surv.cox.relapse.incl.G4.cat") ### added 25/10/18

cox.PFS.cont.G4.df <- extract.cox(results.master, "surv.cox.relapse.incl.G4.contin") ### added 25/10/18

cox.PFS.cat.SHH.df <- extract.cox (results.master, "surv.cox.relapse.incl.SHH.cat")

cox.PFS.cont.SHH.df <- extract.cox (results.master, "surv.cox.relapse.incl.SHH.contin")

cox.PFS.cat.SHH.old.df <- extract.cox (results.master, "surv.cox.relapse.incl.SHH.old.cat") 
                                       
cox.PFS.cont.SHH.old.df <- extract.cox (results.master, "surv.cox.relapse.incl.SHH.old.contin") 

### generate dataset for GSEA/IPA analysis with cox univariate regression using ranked.file function 
### note sometimes need to rerun ranked.file function (in extract functions) or annotate.HTseq.IDs function(in master function) or rerun part of command 
### only useful if gene is a known transcript

annot.cox.PFS.cat.all.df <- ranked.file (cox.PFS.cat.all.df)
annot.cox.PFS.cont.all.df <- ranked.file (cox.PFS.cont.all.df)
annot.cox.PFS.cat.G3G4.df <- ranked.file (cox.PFS.cat.G3G4.df)
annot.cox.PFS.cont.G3G4.df <- ranked.file (cox.PFS.cont.G3G4.df)
annot.cox.PFS.cat.G4.df <- ranked.file (cox.PFS.cat.G4.df) ### added 25/10/18
annot.cox.PFS.cont.G4.df <- ranked.file (cox.PFS.cont.G4.df) ### added 25/10/18
annot.cox.PFS.cat.SHH.df <- ranked.file (cox.PFS.cat.SHH.df)
annot.cox.PFS.cont.SHH.df <- ranked.file (cox.PFS.cont.SHH.df)
annot.cox.PFS.cat.SHH.old.df <- ranked.file (cox.PFS.cat.SHH.old.df)
annot.cox.PFS.cont.SHH.old.df <- ranked.file (cox.PFS.cont.SHH.old.df)

write.csv(annot.cox.PFS.cat.all.df, file = "/home/nmm199/R/MB_RNAseq/Clinical/clin.results/Oct_25_2018/Complete_transcripts_filtered/annot.cox.PFS.cat.all.csv")
write.csv (annot.cox.PFS.cont.all.df, file = "/home/nmm199/R/MB_RNAseq/Clinical/clin.results/Oct_25_2018/Complete_transcripts_filtered/annot.cox.PFS.cont.all.csv") 
write.csv (annot.cox.PFS.cat.G3G4.df, file = "/home/nmm199/R/MB_RNAseq/Clinical/clin.results/Oct_25_2018/Complete_transcripts_filtered/annot.cox.PFS.cat.G3G4.csv" )
write.csv (annot.cox.PFS.cont.G3G4.df, file = "/home/nmm199/R/MB_RNAseq/Clinical/clin.results/Oct_25_2018/Complete_transcripts_filtered/annot.cox.PFS.cont.G3G4.csv" )
write.csv (annot.cox.PFS.cat.G4.df, file = "/home/nmm199/R/MB_RNAseq/Clinical/clin.results/Oct_25_2018/Complete_transcripts_filtered/annot.cox.PFS.cat.G4.csv" )
write.csv (annot.cox.PFS.cont.G4.df, file = "/home/nmm199/R/MB_RNAseq/Clinical/clin.results/Oct_25_2018/Complete_transcripts_filtered/annot.cox.PFS.cont.G4.csv" )
write.csv (annot.cox.PFS.cat.SHH.df, file = "/home/nmm199/R/MB_RNAseq/Clinical/clin.results/Oct_25_2018/Complete_transcripts_filtered/annot.cox.PFS.cat.SHH.csv" )
write.csv (annot.cox.PFS.cont.SHH.df, file = "/home/nmm199/R/MB_RNAseq/Clinical/clin.results/Oct_25_2018/Complete_transcripts_filtered/annot.cox.PFS.cont.SHH.csv" )
write.csv (annot.cox.PFS.cat.SHH.old.df, file = "/home/nmm199/R/MB_RNAseq/Clinical/clin.results/Oct_25_2018/Complete_transcripts_filtered/annot.cox.PFS.cat.SHH.old.csv" )
write.csv (annot.cox.PFS.cont.SHH.old.df , file = "/home/nmm199/R/MB_RNAseq/Clinical/clin.results/Oct_25_2018/Complete_transcripts_filtered/annot.cox.PFS.cont.SHH.old.csv" )

##########################################################################################################################
sig.cox.PFS.cat.all <- cox.PFS.cat.all.df [which(cox.PFS.cat.all.df[, 2]<0.05),]  ### nrow = 81 for complete transcripts

sig.cox.PFS.cont.all <- cox.PFS.cont.all.df[which(cox.PFS.cont.all.df[, 2]<0.05),] ### nrow = 698 for complete transcripts

# write.csv (sig.cox.PFS.cont.all, file = "/home/nmm199/R/MB_RNAseq/Clinical/clin.results/Oct_25_2018/Complete_transcripts_filtered/sig.cox.PFS.cont.all.csv")

sig.cox.PFS.cat.G3G4 <- cox.PFS.cat.G3G4.df[which(cox.PFS.cat.G3G4.df[, 2]<0.05),]

sig.cox.PFS.cont.G3G4 <- cox.PFS.cont.G3G4.df[which(cox.PFS.cont.G3G4.df[, 2]<0.05),]

sig.cox.PFS.cat.G4 <- cox.PFS.cat.G4.df[which (cox.PFS.cat.G4.df [,2]<0.05),]

sig.cox.PFS.cont.G4 <- cox.PFS.cont.G4.df [which (cox.PFS.cont.G4.df [,2]<0.05),]

sig.cox.PFS.cat.SHH <- cox.PFS.cat.SHH.df[which(cox.PFS.cat.SHH.df[, 2]<0.05),]

sig.cox.PFS.cont.SHH <- cox.PFS.cont.SHH.df[which(cox.PFS.cont.SHH.df[, 2]<0.05)] ### there is an error here I believe 24/4/18

sig.cox.PFS.cat.SHH.old <- cox.PFS.cat.SHH.old.df[which(cox.PFS.cat.SHH.old.df[, 2]<0.05),] ### no significant results 11/10/17

sig.cox.PFS.cont.SHH.old <- cox.PFS.cont.SHH.old.df[which(cox.PFS.cont.SHH.old.df[, 2]<0.05),] ### makes more sense as continuous variable 

###########################################################################################

### annotate those with ensembl gene IDs, removing those with NA ### need to run annotate.HTseq.IDs function again prior to this working
### remember useful for known transcripts only
### need to run annotate.HTseq.IDs function again prior to this working

# try(annot.cox.PFS.cont.all <- annotate.HTseq.IDs(rownames(cox.PFS.cont.all.df)), silent = T) 
try(annot.sig.cox.PFS.cont.all <- annotate.HTseq.IDs(rownames(sig.cox.PFS.cont.all)),silent = T) ### worked after reload load a few times of function 24/4/18
try(annot.sig.cox.PFS.cont.all.df <- cbind (sig.cox.PFS.cont.all, annot.sig.cox.PFS.cont.all), silent = T) ### sometimes works better without "try" statement ###file contains significant results and annotated gene name

try(annot.sig.cox.PFS.cat.all <- annotate.HTseq.IDs(rownames(sig.cox.PFS.cat.all)), silent = T)
try(annot.sig.cox.PFS.cat.all.df <- cbind(sig.cox.PFS.cat.all, annot.sig.cox.PFS.cat.all), silent = T)


try(annot.cox.PFS.cont.G3G4 <- annotate.HTseq.IDs(rownames(sig.cox.PFS.cont.G3G4)), silent = T)
annot.sig.cox.PFS.cont.G3G4 <- cbind(annot.cox.PFS.cont.G3G4, sig.cox.PFS.cont.G3G4)
# annot.sig.cox.PFS.cont.G3G4.clean <- annot.cox.PFS.cont.G3G4[complete.cases(annot.sig.cox.PFS.cont.G3G4),]


try(annot.cox.PFS.cat.G3G4 <- annotate.HTseq.IDs(rownames(sig.cox.PFS.cat.G3G4)), silent = T) ### error if not data
annot.sig.cox.PFS.cat.G3G4 <- cbind(annot.cox.PFS.cat.G3G4, sig.cox.PFS.cat.G3G4)

try(annot.cox.PFS.cat.G4 <- annotate.HTseq.IDs(rownames(sig.cox.PFS.cat.G4)), silent = T) ### may not create any transcripts, added 25/10/18

try(annot.cox.PFS.cont.G4 <- annotate.HTseq.IDs(rownames(sig.cox.PFS.cont.G4)), silent = T) ### may not create any transcripts, added 25/10/18

# try(annot.cox.PFS.cont.SHH <- annotate.HTseq.IDs(rownames(sig.cox.PFS.cont.SHH)), silent = T) ### error as does not contain data
# annot.sig.cox.PFS.cont.SHH <- cbind(annot.cox.PFS.cont.SHH, sig.cox.PFS.cont.SHH)
# try(annot.cox.PFS.cat.SHH <- annotate.HTseq.IDs(rownames(sig.cox.PFS.cat.SHH)), silent = T)  ### error as does not contain data

try(annot.cox.PFS.cont.SHH.old <- annotate.HTseq.IDs(rownames(sig.cox.PFS.cont.SHH.old)), silent = T) ### error, not reading annotate.HTseq.IDs properly
annot.sig.cox.PFS.cont.SHH.old <- cbind (annot.cox.PFS.cont.SHH.old, sig.cox.PFS.cont.SHH.old)

try(annot.cox.PFS.cat.SHH.old <- annotate.HTseq.IDs(rownames(sig.cox.PFS.cat.SHH.old)), silent = T)
annot.sig.cox.PFS.cat.SHH.old <- cbind(annot.cox.PFS.cat.SHH.old, sig.cox.PFS.cat.SHH.old)

### save annotated files
 
### unhash SHH cat and SHH contin if useful
# clean.annot.sig.cox.PFS.cont.all <- annot.sig.cox.PFS.cont.all[complete.cases(annot.sig.cox.PFS.cont.all),] ### complete.cases removes NAs

try(write.csv(annot.sig.cox.PFS.cont.all.df, file = "/home/nmm199/R/MB_RNAseq/Clinical/clin.results/Oct_25_2018/Complete_transcripts_filtered/annot.sig.cox.PFS.cont.all.csv"), silent = T)
### note that this file above contains both the significant results and the annotated gene name

write.csv(annot.sig.cox.PFS.cat.all.df, file = "/home/nmm199/R/MB_RNAseq/Clinical/clin.results/Oct_25_2018/Complete_transcripts_filtered/annot.sig.cox.PFS.cat.all.csv")

write.csv(annot.sig.cox.PFS.cont.G3G4, file = "/home/nmm199/R/MB_RNAseq/Clinical/clin.results/Oct_25_2018/Complete_transcripts_filtered/annot.sig.cox.PFS.cont.G3G4.csv")
write.csv(annot.sig.cox.PFS.cat.G3G4, file = "/home/nmm199/R/MB_RNAseq/Clinical/clin.results/Oct_25_2018/Complete_transcripts_filtered/annot.sig.cox.PFS.cat.G3G4.csv")


# write.csv(annot.sig.cox.PFS.cont.SHH, file = "/home/nmm199/R/MB_RNAseq/Clinical/clin.results/Oct_25_2018/Complete_transcripts_filtered/significant.cox.PFS.cont.SHH.csv")
# write.csv (annot.sig.cox.PFS.cat.SHH, file = ""/home/nmm199/R/MB_RNAseq/Clinical/clin.results/Oct_25_2018/Complete_transcripts_filtered/significant.cox.PFS.cat.SHH.csv")

try(write.csv (annot.sig.cox.PFS.cont.SHH.old, file = "/home/nmm199/R/MB_RNAseq/Clinical/clin.results/Oct_25_2018/Complete_transcripts_filtered/annot.sig.cox.PFS.cont.SHH.old.csv"), silent = T)
try(write.csv(annot.sig.cox.PFS.cat.SHH.old, file = "/home/nmm199/R/MB_RNAseq/Clinical/clin.results/Oct_25_2018/Complete_transcripts_filtered/annot.sig.cox.PFS.cat.SHH.old.csv"), silent = T)

### for files where annotate did not work or not feasible

write.csv(sig.cox.PFS.cont.SHH.old, file = "/home/nmm199/R/MB_RNAseq/Clinical/clin.results/Oct_25_2018/Complete_transcripts_filtered/sig.cox.PFS.cont.SHHold.csv")
write.csv(sig.cox.PFS.cont.SHH, file = "/home/nmm199/R/MB_RNAseq/Clinical/clin.results/Oct_25_2018/Complete_transcripts_filtered/sig.cox.PFS.cont.SHH.csv") 

########################################################################################
########################################################################################

### Cox OS overall for categorical variable 

cox.OS.cat.all.df <- extract.cox.OS (results.master, "surv.cox.result.OS.all.cat") ### was 11; and so forth below, note different function for extract.cox.OS

cox.OS.cont.all.df <- extract.cox.OS (results.master, "surv.cox.result.OS.all.contin")  

cox.OS.cat.G3G4.df <- extract.cox.OS (results.master, "surv.cox.result.OS.G3G4.cat" )  

cox.OS.cont.G3G4.df <- extract.cox.OS (results.master, "surv.cox.result.OS.G3G4.contin")  

cox.OS.cat.G4.df <- extract.cox.OS (results.master, "surv.cox.results.OS.G4.cat") ### 25/10/18

cox.OS.cont.G4.df <- extract.cox.OS (results.master, "surv.cox.results.OS.G4.contin") ### 25/10/18

cox.OS.cat.SHH.df <- extract.cox (results.master, "surv.cox.result.OS.SHH.cat")

cox.OS.cont.SHH.df <- extract.cox (results.master, "surv.cox.result.OS.SHH.contin")

cox.OS.cat.SHH.old.df <- extract.cox (results.master, "surv.cox.result.OS.SHH.old.cat")

cox.OS.cont.SHH.old.df <- extract.cox (results.master, "surv.cox.result.OS.SHH.old.contin")

##################
### generate complete datasets for ranking in GSEA or IPA
### easiest if generated annotated dataframes for export with gene names if known transcripts

annot.cox.OS.cat.all.df <- ranked.file (cox.OS.cat.all.df)
annot.cox.OS.cont.all.df <- ranked.file (cox.OS.cont.all.df)
annot.cox.OS.cat.G3G4.df <- ranked.file (cox.OS.cat.G3G4.df)
annot.cox.OS.cont.G3G4.df <- ranked.file (cox.OS.cont.G3G4.df)
annot.cox.OS.cat.G4.df <- ranked.file (cox.OS.cat.G4.df) ### 25/10/18
annot.cox.OS.cont.G4.df <- ranked.file (cox.OS.cont.G4.df) ### 25/10/18
annot.cox.OS.cat.SHH.df <- ranked.file (cox.OS.cat.SHH.df)
annot.cox.OS.cont.SHH.df <- ranked.file (cox.OS.cont.SHH.df)
annot.cox.OS.cat.SHH.old.df <- ranked.file (cox.OS.cat.SHH.old.df)
annot.cox.OS.cont.SHH.old.df <- ranked.file (cox.OS.cont.SHH.old.df)

write.csv(annot.cox.OS.cat.all.df, file = "/home/nmm199/R/MB_RNAseq/Clinical/clin.results/Oct_25_2018/Complete_transcripts_filtered/annot.cox.OS.cat.all.csv")
write.csv (annot.cox.OS.cont.all.df, file = "/home/nmm199/R/MB_RNAseq/Clinical/clin.results/Oct_25_2018/Complete_transcripts_filtered/annot.cox.OS.cont.all.csv") 
write.csv (annot.cox.OS.cat.G3G4.df, file = "/home/nmm199/R/MB_RNAseq/Clinical/clin.results/Oct_25_2018/Complete_transcripts_filtered/annot.cox.OS.cat.G3G4.csv" )
write.csv (annot.cox.OS.cont.G3G4.df, file = "/home/nmm199/R/MB_RNAseq/Clinical/clin.results/Oct_25_2018/Complete_transcripts_filtered/annot.cox.OS.cont.G3G4.csv" )
write.csv (annot.cox.OS.cat.G4.df, file = "/home/nmm199/R/MB_RNAseq/Clinical/clin.results/Oct_25_2018/Complete_transcripts_filtered/annot.cox.OS.cat.G4.csv" )
write.csv (annot.cox.OS.cont.G4.df, file = "/home/nmm199/R/MB_RNAseq/Clinical/clin.results/Oct_25_2018/Complete_transcripts_filtered/annot.cox.OS.cont.G4.csv" )
write.csv (annot.cox.OS.cat.SHH.df, file = "/home/nmm199/R/MB_RNAseq/Clinical/clin.results/Oct_25_2018/Complete_transcripts_filtered/annot.cox.OS.cat.SHH.csv" )
write.csv (annot.cox.OS.cont.SHH.df, file = "/home/nmm199/R/MB_RNAseq/Clinical/clin.results/Oct_25_2018/Complete_transcripts_filtered/annot.cox.OS.cont.SHH.csv" )
write.csv (annot.cox.OS.cat.SHH.old.df, file = "/home/nmm199/R/MB_RNAseq/Clinical/clin.results/Oct_25_2018/Complete_transcripts_filtered/annot.cox.OS.cat.SHH.old.csv" )
write.csv (annot.cox.OS.cont.SHH.old.df , file = "/home/nmm199/R/MB_RNAseq/Clinical/clin.results/Oct_25_2018/Complete_transcripts_filtered/annot.cox.OS.cont.SHH.old.csv" )

### significant dataframes & annotation

sig.cox.OS.cat.all <- cox.OS.cat.all.df[which(cox.OS.cat.all.df[, 2]<0.05),]
try(annot.sig.cox.OS.cat.all <- annotate.HTseq.IDs(rownames(sig.cox.OS.cat.all)), silent = T)
write.csv(sig.cox.OS.cat.all,  file = "/home/nmm199/R/MB_RNAseq/Clinical/clin.results/Oct_25_2018/Complete_transcripts_filtered/sig.cox.OS.cat.all.csv")

###

sig.cox.OS.cont.all <- cox.OS.cont.all.df[which(cox.OS.cont.all.df[, 2]<0.05),]
try(annot.sig.cox.OS.cont.all <- annotate.HTseq.IDs(rownames(sig.cox.OS.cont.all)), silent = T)
write.csv(sig.cox.OS.cont.all,  file = "/home/nmm199/R/MB_RNAseq/Clinical/clin.results/Oct_25_2018/Complete_transcripts_filtered/sig.cox.OS.cont.all.csv")
write.csv(annot.sig.cox.OS.cont.all, file = "/home/nmm199/R/MB_RNAseq/Clinical/clin.results/Oct_25_2018/Complete_transcripts_filtered/annot.sig.cox.OS.cont.all.csv" )
annot.sig.cox.OS.cont.all.df <- cbind(sig.cox.OS.cont.all, annot.sig.cox.OS.cont.all)
write.csv (annot.sig.cox.OS.cont.all.df, file = "/home/nmm199/R/MB_RNAseq/Clinical/clin.results/Oct_25_2018/Complete_transcripts_filtered/annot.sig.cox.OS.cont.all.complete.csv" )
###

sig.cox.OS.cat.G3G4 <- cox.OS.cat.G3G4.df[which(cox.OS.cat.G3G4.df[, 2]<0.05),]  ### no results
try(annot.sig.cox.OS.cat.G3G4 <- annotate.HTseq.IDs(rownames(sig.cox.OS.cat.G3G4)), silent = T)
write.csv(sig.cox.OS.cat.G3G4,  file = "/home/nmm199/R/MB_RNAseq/Clinical/clin.results/Oct_25_2018/Complete_transcripts_filtered/sig.cox.OS.cat.G3G4.csv")
annot.sig.cox.OS.cat.G3G4.df <- cbind (sig.cox.OS.cat.G3G4, annot.sig.cox.OS.cat.G3G4)
write.csv(annot.sig.cox.OS.cat.G3G4.df, file = "/home/nmm199/R/MB_RNAseq/Clinical/clin.results/Oct_25_2018/Complete_transcripts_filtered/annot.sig.cox.OS.cat.G3G4.complete.csv")
###

sig.cox.OS.cont.G3G4 <- cox.OS.cont.G3G4.df[which(cox.OS.cont.G3G4.df[, 2]<0.05),]
try(annot.sig.cox.OS.cont.G3G4 <- annotate.HTseq.IDs(rownames(sig.cox.OS.cont.G3G4)), silent = T)
annot.sig.cox.OS.cont.G3G4.df <- cbind (sig.cox.OS.cont.G3G4, annot.sig.cox.OS.cont.G3G4)
write.csv(sig.cox.OS.cont.G3G4,  file = "/home/nmm199/R/MB_RNAseq/Clinical/clin.results/Oct_25_2018/Complete_transcripts_filtered/sig.cox.OS.cont.G3G4.csv")
write.csv(annot.sig.cox.OS.cont.G3G4.df, file = "/home/nmm199/R/MB_RNAseq/Clinical/clin.results/Oct_25_2018/Complete_transcripts_filtered/annot.sig.cox.OS.cont.G3G4.csv")

###

sig.cox.OS.cat.SHH <- cox.OS.cat.SHH.df[which(cox.OS.cat.SHH.df[, 2]<0.05),] 

sig.cox.OS.cont.SHH.df <- cox.OS.cont.SHH.df[which(cox.OS.cont.SHH.df[, 2]<0.05), ]

write.csv(sig.cox.OS.cont.SHH.df, file = "/home/nmm199/R/MB_RNAseq/Clinical/clin.results/Oct_25_2018/Complete_transcripts_filtered/sig.cox.OS.cont.SHH.csv")
try(annot.sig.cox.OS.cont.SHH.df <- annotate.HTseq.IDs(rownames(sig.cox.OS.cont.SHH.df)),silent = T) ### this worked after reran annotate.HTseq.IDs function within clinical_data_functions_master.R
try(write.csv(annot.sig.cox.OS.cont.SHH.df, file = "/home/nmm199/R/MB_RNAseq/Clinical/clin.results/Oct_25_2018/Complete_transcripts_filtered/annot.sig.cox.OS.cont.SHH.csv"), silent = TRUE)

sig.cox.OS.cat.SHH.old <- cox.OS.cat.SHH.old.df[which(cox.OS.cat.SHH.old.df[, 2]<0.05),]
try(write.csv(sig.cox.OS.cat.SHH.old, file = "/home/nmm199/R/MB_RNAseq/Clinical/clin.results/Oct_25_2018/Complete_transcripts_filtered/sig.cox.OS.cat.SHHold.csv"), silent = T)

sig.cox.OS.cont.SHH.old.df <- cox.OS.cont.SHH.old.df [which(cox.OS.cont.SHH.old.df[,2]<0.05),]
try(write.csv(sig.cox.OS.cont.SHH.old.df, file = "/home/nmm199/R/MB_RNAseq/Clinical/clin.results/Oct_25_2018/Complete_transcripts_filtered/sig.cox.OS.cont.SHHold.csv"), silent = T)
try(annot.sig.cox.OS.cont.SHH.old <- annotate.HTseq.IDs(rownames(sig.cox.OS.cont.SHH.old.df)), silent = T)
try(write.csv(annot.sig.cox.OS.cont.SHH.old, file = "/home/nmm199/R/MB_RNAseq/Clinical/clin.results/Oct_25_2018/Complete_transcripts_filtered/annot.sig.cox.OS.cont.SHHold.csv"), silent = T)
########################################################################

### Cox EFS for all - these are all categorical expression data

# cox.EFS.cat.all.df <- extract.cox.OS (results.master, 1) ### previous script which is now replaced by extract.cox function and specific names
cox.EFS.cat.all.df <- extract.cox(results.master, "surv.cox.EFS.incl") 

cox.EFS.cat.G3G4.df <- extract.cox(results.master, "surv.cox.EFS.incl.G3G4")


sig.cox.EFS.cat.all <- cox.EFS.cat.all.df[which(cox.EFS.cat.all.df[, 2]<0.05),]

sig.cox.EFS.cat.G3G4 <- cox.EFS.cat.G3G4.df[which(cox.EFS.cat.G3G4.df[,2]<0.05), ]


try(write.csv(sig.cox.EFS.cat.all,  file = "/home/nmm199/R/MB_RNAseq/Clinical/clin.results/Oct_25_2018/Complete_transcripts_filtered/sig.cox.EFS.all.csv"), silent =T)

try(write.csv(sig.cox.EFS.cat.G3G4,  file = "/home/nmm199/R/MB_RNAseq/Clinical/clin.results/Oct_25_2018/Complete_transcripts_filtered/sig.cox.EFS.G3G4.csv"), silent = T)

#################################################################################################
#################################################################################################

### extract logistic regression p value 

### decision made to hard code as difficulty generating function, in hindsight could replace with "reg.log.list" like in chi square function (extract.chi.all, name = "....")

logreg.LCA.pval <- lapply(results.master, function(x){return(x$reg.log.list$log.reg.LCA$LR.pval)})  ### x$reg.log.list$log.reg.LCA[[1]]
logreg.LCA.OR <- lapply(results.master, function(x){return(x$reg.log.list$log.reg.LCA$LR.OR.val)})
logreg.LCA.L95CI <- lapply(results.master, function(x){return(x$reg.log.list$log.reg.LCA$lower.95CI)})
logreg.LCA.U95CI <- lapply(results.master, function(x){return(x$reg.log.list$log.reg.LCA$upper.95CI)})
extract.logreg.LCA.df <- log.reg.dataframe(pval = logreg.LCA.pval, OR = logreg.LCA.OR, L95CI = logreg.LCA.L95CI, U95CI= logreg.LCA.U95CI)

###

logreg.relapse.pval <- lapply(results.master, function(x){return(x$reg.log.list$log.reg.relapse$LR.pval)})
logreg.relapse.OR <- lapply(results.master, function(x){return(x$reg.log.list$log.reg.relapse$LR.OR.val)})
logreg.relapse.L95CI <- lapply(results.master, function(x){return(x$reg.log.list$log.reg.relapse$lower.95CI)})
logreg.relapse.U95CI <- lapply(results.master, function(x){return(x$reg.log.list$log.reg.relapse$upper.95CI)})
extract.logreg.relapse.df <- log.reg.dataframe(pval = logreg.relapse.pval, OR = logreg.relapse.OR, L95CI = logreg.relapse.L95CI, U95CI = logreg.relapse.U95CI)  

###
logreg.mstatus.pval <-lapply(results.master, function(x){return(x$reg.log.list$log.reg.mstatus$LR.pval)}) 
logreg.mstatus.OR <- lapply(results.master, function(x){return(x$reg.log.list$log.reg.mstatus$LR.OR.val)}) 
logreg.mstatus.L95CI <- lapply(results.master, function(x){return(x$reg.log.list$log.reg.mstatus$lower.95CI)}) 
logreg.mstatus.U95CI <- lapply(results.master, function(x){return(x$reg.log.list$log.reg.mstatus$upper.95CI)})
extract.logreg.mstatus.df <- log.reg.dataframe(pval = logreg.mstatus.pval, OR = logreg.mstatus.OR, L95CI = logreg.mstatus.L95CI, U95CI = logreg.mstatus.U95CI)

###
logreg.age.cat.pval <-lapply(results.master, function(x){return(x$reg.log.list$log.reg.age.cat$LR.pval)}) 
logreg.age.cat.OR <- lapply(results.master, function(x){return(x$reg.log.list$log.reg.age.cat$LR.OR.val)}) 
logreg.age.cat.L95CI <- lapply(results.master, function(x){return(x$reg.log.list$log.reg.age.cat$lower.95CI)}) 
logreg.age.cat.U95CI <- lapply(results.master, function(x){return(x$reg.log.list$log.reg.age.cat$upper.95CI)})
extract.logreg.age.cat.df <- log.reg.dataframe(pval = logreg.age.cat.pval , OR = logreg.age.cat.OR , L95CI = logreg.age.cat.L95CI, U95CI = logreg.age.cat.U95CI)

###
logreg.meth.pval <-lapply(results.master, function(x){return(x$reg.log.list$log.reg.meth$LR.pval)}) 
logreg.meth.OR <- lapply(results.master, function(x){return(x$reg.log.list$log.reg.meth$LR.OR.val)}) 
logreg.meth.L95CI <- lapply(results.master, function(x){return(x$reg.log.list$log.reg.meth$lower.95CI)}) 
logreg.meth.U95CI <- lapply(results.master, function(x){return(x$reg.log.list$log.reg.meth$upper.95CI)})
extract.logreg.meth.df <- log.reg.dataframe(pval = logreg.meth.pval , OR = logreg.meth.OR , L95CI = logreg.meth.L95CI, U95CI = logreg.meth.U95CI)

###

logreg.meth7.pval <-lapply(results.master, function(x){return(x$reg.log.list$log.reg.meth7$LR.pval)}) 
logreg.meth7.OR <- lapply(results.master, function(x){return(x$reg.log.list$log.reg.meth7$LR.OR.val)}) 
logreg.meth7.L95CI <- lapply(results.master, function(x){return(x$reg.log.list$log.reg.meth7$lower.95CI)}) 
logreg.meth7.U95CI <- lapply(results.master, function(x){return(x$reg.log.list$log.reg.meth7$upper.95CI)})
extract.logreg.meth7.df <- log.reg.dataframe(pval = logreg.meth7.pval , OR = logreg.meth7.OR , L95CI = logreg.meth7.L95CI, U95CI = logreg.meth7.U95CI)


###

logreg.MYC.pval <-lapply(results.master, function(x){return(x$reg.log.list$log.reg.MYC$LR.pval)}) 
logreg.MYC.OR <- lapply(results.master, function(x){return(x$reg.log.list$log.reg.MYC$LR.OR.val)}) 
logreg.MYC.L95CI <- lapply(results.master, function(x){return(x$reg.log.list$log.reg.MYC$lower.95CI)}) 
logreg.MYC.U95CI <- lapply(results.master, function(x){return(x$reg.log.list$log.reg.MYC$upper.95CI)})
extract.logreg.MYC.df <- log.reg.dataframe(pval = logreg.MYC.pval , OR = logreg.MYC.OR , L95CI = logreg.MYC.L95CI, U95CI = logreg.MYC.U95CI)


###

logreg.MYCN.pval <-lapply(results.master, function(x){return(x$reg.log.list$log.reg.MYCN$LR.pval)}) 
logreg.MYCN.OR <- lapply(results.master, function(x){return(x$reg.log.list$log.reg.MYCN$LR.OR.val)}) 
logreg.MYCN.L95CI <- lapply(results.master, function(x){return(x$reg.log.list$log.reg.MYCN$lower.95CI)}) 
logreg.MYCN.U95CI <- lapply(results.master, function(x){return(x$reg.log.list$log.reg.MYCN$upper.95CI)})
extract.logreg.MYCN.df <- log.reg.dataframe(pval = logreg.MYCN.pval , OR = logreg.MYCN.OR , L95CI = logreg.MYCN.L95CI, U95CI = logreg.MYCN.U95CI)

###
logreg.MYCMYCN.pval <-lapply(results.master, function(x){return(x$reg.log.list$log.reg.MYCMYCN$LR.pval)}) 
logreg.MYCMYCN.OR <- lapply(results.master, function(x){return(x$reg.log.list$log.reg.MYCMYCN$LR.OR.val)}) 
logreg.MYCMYCN.L95CI <- lapply(results.master, function(x){return(x$reg.log.list$log.reg.MYCMYCN$lower.95CI)}) 
logreg.MYCMYCN.U95CI <- lapply(results.master, function(x){return(x$reg.log.list$log.reg.MYCMYCN$upper.95CI)})
extract.logreg.MYCMYCN.df <- log.reg.dataframe(pval = logreg.MYCMYCN.pval , OR = logreg.MYCMYCN.OR , L95CI = logreg.MYCMYCN.L95CI, U95CI = logreg.MYCMYCN.U95CI)

###

logreg.resection.pval <- lapply(results.master, function(x){return(x$reg.log.list$log.reg.resection$LR.pval)})
logreg.resection.OR <- lapply(results.master, function(x){return(x$reg.log.list$log.reg.resection$LR.OR.val)})
logreg.resection.L95CI <- lapply(results.master, function(x){return(x$reg.log.list$log.reg.resection$lower.95CI)})
logreg.resection.U95CI <- lapply(results.master, function(x){return(x$reg.log.list$log.reg.resection$upper.95CI)})
extract.logreg.resection.df <- log.reg.dataframe(pval = logreg.resection.pval, OR = logreg.resection.OR, L95CI = logreg.resection.L95CI, U95CI = logreg.resection.U95CI) 

###
logreg.sex.pval <- lapply(results.master, function(x){return(x$reg.log.list$log.reg.sex$LR.pval)})
logreg.sex.OR <- lapply(results.master, function(x){return(x$reg.log.list$log.reg.sex$LR.OR.val)})
logreg.sex.L95CI <- lapply(results.master, function(x){return(x$reg.log.list$log.reg.sex$lower.95CI)})
logreg.sex.U95CI <- lapply(results.master, function(x){return(x$reg.log.list$log.reg.sex$upper.95CI)})
extract.logreg.sex.df <- log.reg.dataframe(pval = logreg.sex.pval, OR = logreg.sex.OR, L95CI = logreg.sex.L95CI, U95CI = logreg.sex.U95CI) 

###

logreg.TERT.pval <- lapply(results.master, function(x){return(x$reg.log.list$log.reg.TERT$LR.pval)})
logreg.TERT.OR <- lapply(results.master, function(x){return(x$reg.log.list$log.reg.TERT$LR.OR.val)})
logreg.TERT.L95CI <- lapply(results.master, function(x){return(x$reg.log.list$log.reg.TERT$lower.95CI)})
logreg.TERT.U95CI <- lapply(results.master, function(x){return(x$reg.log.list$log.reg.TERT$upper.95CI)})
extract.logreg.TERT.df <- log.reg.dataframe(pval = logreg.TERT.pval, OR = logreg.TERT.OR, L95CI = logreg.TERT.L95CI, U95CI = logreg.TERT.U95CI) 


###

logreg.TP53.pval <- lapply(results.master, function(x){return(x$reg.log.list$log.reg.TP53$LR.pval)})
logreg.TP53.OR <- lapply(results.master, function(x){return(x$reg.log.list$log.reg.TP53$LR.OR.val)})
logreg.TP53.L95CI <- lapply(results.master, function(x){return(x$reg.log.list$log.reg.TP53$lower.95CI)})
logreg.TP53.U95CI <- lapply(results.master, function(x){return(x$reg.log.list$log.reg.TP53$upper.95CI)})
extract.logreg.TP53.df <- log.reg.dataframe(pval = logreg.TP53.pval, OR = logreg.TP53.OR, L95CI = logreg.TP53.L95CI, U95CI = logreg.TP53.U95CI) 


######################################################################
### if wish to make a list of all the logistic regression results 
logistic.reg.results <- as.list(mget(ls(pattern="extract.logreg"))) 

######################################################################

### creating log reg list for all significant variables
# significant.logreg.df.all <- logistic.reg.results[which(logistic.reg.results[, 2]<0.05),] ### this does not work because logistic.reg.results object is a list
# logistic.reg.df <- as.data.frame(logistic.reg.results)

#######################################################################

### extract chi square p value

chi.age.cat.infant.result <- extract.chi.all(results.master, name = "list.age.cat.infant")
chi.CSI.result <- extract.chi.all(results.master, name = "list.CSI")
chi.LCA.result <- extract.chi.all(results.master, name = "list.LCA")
chi.meth4.result <-extract.chi.all(results.master, name = "list.meth.4")
chi.meth7.result <- extract.chi.all(results.master, name = "list.meth.7")
chi.mstatus.result <- extract.chi.all(results.master, name = "list.mstatus")
chi.MYC.result <- extract.chi.all(results.master, name = "list.MYC")
chi.MYCMYCN.result <- extract.chi.all(results.master, name = "list.MYCMYCN")
chi.MYCN.result <- extract.chi.all(results.master, name = "list.MYCN")
chi.q13loss.result <- extract.chi.all(results.master, name = "list.q13loss")
chi.relapse.result <- extract.chi.all(results.master, name = "list.relapse")
chi.resection.result <- extract.chi.all(results.master, name = "list.resection")
chi.RTX.result <- extract.chi.all(results.master, name = "list.RTX")
chi.sex.result <- extract.chi.all(results.master, name = "list.sex")
chi.TERT.result <- extract.chi.all(results.master, name = "list.TERT")
chi.TP53.result <- extract.chi.all(results.master, name = "list.TP53")

### extract adj p <0.05 for relapse, mstatus, MYC, MYCN, MYCMYCN

significant.chi.relapse <- chi.relapse.result[which(chi.relapse.result[,2]<0.05), ]  ### n=4388 4/12/17 for mb.vsd
significant.chi.mstatus <- chi.mstatus.result[which(chi.mstatus.result[,2]<0.05), ]  ### n=3875, 4/12/17 for mb.vsd
significant.chi.MYC <- chi.MYC.result [which(chi.MYC.result[,2]<0.05), ]            ### n=4640, 4/12/17 for mb.vsd
significant.chi.MYCN <- chi.MYCN.result[which(chi.MYCN.result[,2]<0.05),]
significant.chi.MYCMYCN <- chi.MYCMYCN.result[which(chi.MYCMYCN.result[,2]<0.05), ]  ### n=214 4/12/17 for mb.vsd


########################################################################

###  multivariate cox, looking for transcripts that are significant beyond either the current PNET5, the Lancet oncology paper (Schwalbe et al 2017) or a combined model taking both models together
### not adjusted for age (25/10/18 new analysis)

multivar.cox.OS.combined.cat.df <- extract.multivar.cox(results.master,  "multivar.cox.OS.combined.cat" )  ### section updated 6/2/18 ### prev updated so that p value is for biomarker not overall model p val 21/11/17

multivar.cox.OS.combined.cont.df <- extract.multivar.cox(results.master,"multivar.cox.OS.combined.contin")  

multivar.cox.OS.PNET.G3G4.cat.df <- extract.multivar.cox(results.master, "multivar.cox.OS.PNET5.G3G4.cat") ### added 25/10/18

multivar.cox.OS.PNET.G3G4.cont.df <- extract.multivar.cox(results.master, "multivar.cox.OS.PNET5.G3G4.contin") ### added 25/10/18

multivar.cox.OS.lancetG3G4.cat.df <- extract.multivar.cox(results.master, "multivar.cox.OS.Lancet.G3G4.cat") 

multivar.cox.OS.lancetG3G4.cont.df <- extract.multivar.cox(results.master, "multivar.cox.OS.Lancet.G3G4.contin" )

multivar.cox.OS.PNET5.cat.df <- extract.multivar.cox(results.master, "multivar.cox.OS.PNET5.cat")

multivar.cox.OS.PNET5.cont.df <- extract.multivar.cox(results.master, "multivar.cox.OS.PNET5.contin") 

multivar.cox.OS.SHHold.cat.df <- extract.multivar.cox(results.master, "multivar.cox.OS.SHH.old.cat") ### something to sort out 20/02/18 as generating 0/inf for 95CI and same p val/adj p val

multivar.cox.OS.SHHold.cont.df <- extract.multivar.cox(results.master, "multivar.cox.OS.SHH.old.contin") ### as above 20/02/18


### PFS
multivar.cox.PFS.combined.cat.df <- extract.multivar.cox (results.master, "multivar.cox.PFS.combined.cat")   

multivar.cox.PFS.combined.cont.df <- extract.multivar.cox (results.master, "multivar.cox.PFS.combined.contin") 

multivar.cox.PFS.PNET5.G3G4.cat.df <- extract.multivar.cox (results.master, "multivar.cox.PFS.PNET5.G3G4.cat") ### added 25/10/18

multivar.cox.PFS.PNET5.G3G4.cont.df <- extract.multivar.cox (results.master, "multivar.cox.PFS.PNET5.G3G4.contin") ### added 25/10/18

multivar.cox.PFS.lancetG3G4.cat.df <- extract.multivar.cox (results.master, "multivar.cox.PFS.Lancet.G3G4.cat")   

multivar.cox.PFS.lancetG3G4.cont.df <- extract.multivar.cox (results.master, "multivar.cox.PFS.Lancet.G3G4.contin") 

multivar.cox.PFS.PNET5.cat.df <- extract.multivar.cox (results.master,  "multivar.cox.PFS.PNET5.cat" ) 

multivar.cox.PFS.PNET5.cont.df <- extract.multivar.cox (results.master, "multivar.cox.PFS.PNET5.contin") 

multivar.cox.PFS.SHHold.cat.df <- extract.multivar.cox (results.master, "multivar.cox.PFS.SHH.old.cat")  

multivar.cox.PFS.SHHold.cont.df <- extract.multivar.cox (results.master, "multivar.cox.PFS.SHH.old.contin") 


### need to generate annotated files with trimmed ENSG ID names and z scores 

### OS files

annot.multivar.cox.OS.combined.cat <- annotate.HTseq.IDs(rownames(multivar.cox.OS.combined.cat.df))
annot.multivar.cox.OS.combined.cat.bind <- cbind (annot.multivar.cox.OS.combined.cat, multivar.cox.OS.combined.cat.df)



annot.multivar.cox.OS.lancetG3G4.cont <- annotate.HTseq.IDs(rownames(multivar.cox.OS.lancetG3G4.cont.df))
annot.multivar.cox.OS.lancetG3G4.cont.bind <- cbind (annot.multivar.cox.OS.lancetG3G4.cont, multivar.cox.OS.lancetG3G4.cont.df)


### PFS files

annot.multivar.cox.PFS.combined.cat <- annotate.HTseq.IDs(rownames(multivar.cox.PFS.combined.cat.df)) 
annot.multivar.cox.PFS.combined.cat.bind <- cbind (annot.multivar.cox.PFS.combined.cat, multivar.cox.PFS.combined.cat.df)

annot.multivar.cox.PFS.combined.cont <- annotate.HTseq.IDs (rownames(multivar.cox.PFS.combined.cont.df))
annot.multivar.cox.PFS.combined.cont.bind <- cbind (annot.multivar.cox.PFS.combined.cont, multivar.cox.PFS.combined.cont.df)

annot.multivar.cox.PFS.lancetG3G4.cont <- annotate.HTseq.IDs(rownames(multivar.cox.PFS.lancetG3G4.cont.df))
annot.multivar.cox.PFS.lancetG3G4.cont.bind <- cbind (annot.multivar.cox.PFS.lancetG3G4.cont, multivar.cox.PFS.lancetG3G4.cont.df)

annot.multivar.cox.PFS.PNET5.cont <- annotate.HTseq.IDs(rownames(multivar.cox.PFS.PNET5.cont.df))
annot.multivar.cox.PFS.PNET5.cont.bind <- cbind (annot.multivar.cox.PFS.PNET5.cont, multivar.cox.PFS.PNET5.cont.df)
  
annot.multivar.cox.PFS.SHHold.cont <- annotate.HTseq.IDs(rownames(multivar.cox.PFS.SHHold.cont.df))
annot.multivar.cox.PFS.SHHold.cont.bind <- cbind(annot.multivar.cox.PFS.SHHold.cont, multivar.cox.PFS.SHHold.cont.df)

### write out csv to interrogate Z scores for specific targets
### OS files

write.csv (annot.multivar.cox.OS.combined.cat.bind, file = "/home/nmm199/R/MB_RNAseq/Clinical/clin.results/Oct_25_2018/Complete_transcripts_filtered/annot.multivar.cox.OS.cat.bind.csv")

write.csv(annot.multivar.cox.OS.lancetG3G4.cont.bind, file = "/home/nmm199/R/MB_RNAseq/Clinical/clin.results/Oct_25_2018/Complete_transcripts_filtered/annot.multivar.cox.OS.lancetG3G4.ocnt.bind.csv")

### PFS files

write.csv(annot.multivar.cox.PFS.combined.cat.bind, file = "/home/nmm199/R/MB_RNAseq/Clinical/clin.results/Oct_25_2018/Complete_transcripts_filtered/annot.multivar.cox.PFS.combined.cat.bind.csv")

write.csv (annot.multivar.cox.PFS.combined.cont.bind, file = "/home/nmm199/R/MB_RNAseq/Clinical/clin.results/Oct_25_2018/Complete_transcripts_filtered/annot.multivar.cox.PFS.combined.cont.bind.csv")

write.csv(annot.multivar.cox.PFS.lancetG3G4.cont.bind, file = "/home/nmm199/R/MB_RNAseq/Clinical/clin.results/Oct_25_2018/Complete_transcripts_filtered/annot.multivar.cox.PFS.lancetG3G4.cont.bind.csv")

write.csv (annot.multivar.cox.PFS.PNET5.cont.bind, file = "/home/nmm199/R/MB_RNAseq/Clinical/clin.results/Oct_25_2018/Complete_transcripts_filtered/annot.multivar.cox.PFS.PNET5.cont.bind.csv")

write.csv (annot.multivar.cox.PFS.SHHold.cont.bind, file = "/home/nmm199/R/MB_RNAseq/Clinical/clin.results/Oct_25_2018/Complete_transcripts_filtered/annot.multivar.cox.PFS.SHHold.cont.bind.csv")


# write.table(annot.multivar.cox.PFS.combined.cat.bind, file = "/home/nmm199/R/MB_RNAseq/Clinical/clin.results/Oct_25_2018/Complete_transcripts_filtered/annot.multi.cox.PFS.combined.cat.bind.txt", sep = "\t",
#      row.names = TRUE, col.names = TRUE)



### files to annotate, prioritise G3G4

multivar.cox.PFS.PNET5.G3G4.cat.df 

multivar.cox.PFS.PNET5.G3G4.cont.df 

multivar.cox.PFS.lancetG3G4.cat.df 


multivar.cox.PFS.PNET5.cat.df <- 

 

multivar.cox.PFS.SHHold.cat.df 

multivar.cox.PFS.SHHold.cont.df 

### write csv






write.csv (multivar.cox.PFS.combined.cont.df, file = "/home/nmm199/R/MB_RNAseq/Clinical/clin.results/Oct_25_2018/Complete_transcripts_filtered/multivar.cox.PFS.combined.cont.csv")

write.csv (multivar.cox.PFS.combined.cat.df, file = "/home/nmm199/R/MB_RNAseq/Clinical/clin.results/Oct_25_2018/Complete_transcripts_filtered/multivar.cox.PFS.combined.cat.csv")

write.csv (multivar.cox.PFS.lancetG3G4.cat.df, file = "/home/nmm199/R/MB_RNAseq/Clinical/clin.results/Oct_25_2018/Complete_transcripts_filtered/multivar.cox.PFS.lancetG3G4.cat.csv")

write.csv (multivar.cox.PFS.lancetG3G4.cont.df, file = "/home/nmm199/R/MB_RNAseq/Clinical/clin.results/Oct_25_2018/Complete_transcripts_filtered/multivar.cox.PFS.lancetG3G4.cont.csv")


### generating significant dataframes for the multivariate cox modelling, ie transcripts that perform about and beyond current clinical risk models

significant.multivar.cox.OS.combined.cat <- multivar.cox.OS.combined.cat.df [which(multivar.cox.OS.combined.cat.df[,2]<0.05), ]

significant.multivar.cox.OS.combined.cont <- multivar.cox.OS.combined.cont.df [which(multivar.cox.OS.combined.cont.df[,2]<0.05), ] 

significant.multivar.cox.OS.PNET5.G3G4.cat <- multivar.cox.OS.PNET.G3G4.cat.df [which(multivar.cox.OS.PNET.G3G4.cat.df[,2]<0.05), ] ### added 25/10/18

significant.multivar.cox.OS.PNET5.G3G4.cont <- multivar.cox.OS.PNET.G3G4.cont.df [which(multivar.cox.OS.PNET.G3G4.cont.df[,2]<0.05), ] ### added 25/10/18

significant.multivar.cox.OS.lancetG3G4.cat <- multivar.cox.OS.lancetG3G4.cat.df[which(multivar.cox.OS.lancetG3G4.cat.df[,2]<0.05),] ### n=0 6/2/18 for all transcripts (allgenes.20180104.rds)

significant.multivar.cox.OS.lancetG3G4.cont <- multivar.cox.OS.lancetG3G4.cont.df[which(multivar.cox.OS.lancetG3G4.cont.df[,2]<0.05),] ### n= 2 6/2/18

significant.multivar.cox.OS.PNET5.cat <- multivar.cox.OS.PNET5.cat.df[which(multivar.cox.OS.PNET5.cat.df[,2]<0.05),]  ### n=43 6/2/18 and Jan 2018

significant.multivar.cox.OS.PNET5.cont <- multivar.cox.OS.PNET5.cont.df[which(multivar.cox.OS.PNET5.cont.df[,2]<0.05),] ### n=30 Jan 2018, n=48 6/2/18, need to filter real results

significant.multivar.cox.OS.SHHold.cat <- multivar.cox.OS.SHHold.cat.df[which(multivar.cox.OS.SHHold.cat.df[,2]<0.05), ] ### n=0 6/2/18

significant.multivar.cox.OS.SHHold.cont <- multivar.cox.OS.SHHold.cont.df[which(multivar.cox.OS.SHHold.cont.df[,2]<0.05),] ### n=0 6/2/18


### multivariate PFS modelling

significant.multivar.cox.PFS.combined.cat <- multivar.cox.PFS.combined.cat.df [which (multivar.cox.PFS.combined.cat.df[,2]<0.05), ] ### n=0 6/2/18 allgenes.20180104.rds

significant.multivar.cox.PFS.combined.cont <- multivar.cox.PFS.combined.cont.df [which (multivar.cox.PFS.combined.cont.df[,2]<0.05), ]  ### n=0 6/2/18 allgenes.20180104.rds

significant.multivar.cox.PFS.PNET.G3G4.cat <- multivar.cox.PFS.PNET5.G3G4.cat.df[which (multivar.cox.PFS.PNET5.G3G4.cat.df [,2]<0.05), ]
  
significant.multivar.cox.PFS.PNET.G3G4.cont <- multivar.cox.PFS.PNET5.G3G4.cont.df [which(multivar.cox.PFS.PNET5.G3G4.cont.df[,2]<0.05),]
  
significant.multivar.cox.PFS.lancetG3G4.cat <- multivar.cox.PFS.lancetG3G4.cat.df [which (multivar.cox.PFS.lancetG3G4.cat.df[,2]<0.05), ]  ### n=4 6/2/18 allgenes.20180104.rds, n=4 filt.mb.vsd (results.filt.genefilter.20180220.rds)

significant.multivar.cox.PFS.lancetG3G4.cont <- multivar.cox.PFS.lancetG3G4.cont.df[which (multivar.cox.PFS.lancetG3G4.cont.df[,2]<0.05),]  ### n=0 6/2/18 allgenes.20180104.rds

significant.multivar.cox.PFS.PNET5.cat <- multivar.cox.PFS.PNET5.cat.df[which(multivar.cox.PFS.PNET5.cat.df[,2]<0.05),]  ### n=0 6/2/18 allgenes.20180104.rds

significant.multivar.cox.PFS.PNET5.cont <- multivar.cox.PFS.PNET5.cont.df[which(multivar.cox.PFS.PNET5.cont.df[,2]<0.05),] ### n=0 6/2/18 allgenes.20180104.rds

significant.multivar.cox.PFS.SHHold.cat <- multivar.cox.PFS.SHHold.cat.df[which(multivar.cox.PFS.SHHold.cat.df[,2]<0.05),] ### n=0 6/2/18 allgenes.20180104.rds

significant.multivar.cox.PFS.SHHold.cont <- multivar.cox.PFS.SHHold.cont.df[which(multivar.cox.PFS.SHHold.cont.df[,2]<0.05),]  ### n=2, however not realistic 6/2/18 allgenes.20180104.rds

########################################################################
### writing out relevant csv files for significant transcripts on multivariate analyses, with annotation attached
### examples below from gp.filt.mb.vsd file
### annotation does not work if there is no object to annotate; sometimes need to run the script twice and OUTSIDE of the try statement separately, as is temperamental.
### need to run the annotate.HTseq.IDs script twice in order to overcome "error" about "the given dataset: hsapiens_gene_ensembl , is not valid.  


# write.csv(significant.multivar.cox.OS.combined.cat, file = "/home/nmm199/R/MB_RNAseq/Clinical/clin.results/Oct_25_2018/Complete_transcripts_filtered/results.sig.multivar.cox.OS.combined.cat.csv")

try(write.csv(significant.multivar.cox.OS.lancetG3G4.cont, file = "/home/nmm199/R/MB_RNAseq/Clinical/clin.results/Oct_25_2018/Complete_transcripts_filtered/results.sig.multi.cox.OS.lancetG3G4.cont.csv"), silent = T)

write.csv(significant.multivar.cox.PFS.lancetG3G4.cont, file = "/home/nmm199/R/MB_RNAseq/Clinical/clin.results/Oct_25_2018/Complete_transcripts_filtered/results.sig.multivar.cox.PFS.lancetG3G4.cont.csv")

try(write.csv(significant.multivar.cox.PFS.lancetG3G4.cat, file =  "/home/nmm199/R/MB_RNAseq/Clinical/clin.results/Oct_25_2018/Complete_transcripts_filtered/results.sig.multivar.cox.PFS.lancetG3G4.cat.csv"  ), silent = T)

try(write.csv(significant.multivar.cox.PFS.lancetG3G4.cont, file =  "/home/nmm199/R/MB_RNAseq/Clinical/clin.results/Oct_25_2018/Complete_transcripts_filtered/results.sig.multivar.cox.PFS.lancetG3G4.cont.csv"  ), silent = T)


write.csv(significant.multivar.cox.PFS.combined.cont, file = "/home/nmm199/R/MB_RNAseq/Clinical/clin.results/Oct_25_2018/Complete_transcripts_filtered/results.sig.multivar.cox.PFS.combined.cont.csv")

write.csv(significant.multivar.cox.PFS.combined.cat, file = "/home/nmm199/R/MB_RNAseq/Clinical/clin.results/Oct_25_2018/Complete_transcripts_filtered/results.sig.multivar.cox.PFS.combined.cat.csv")

write.csv(significant.multivar.cox.PFS.PNET.G3G4.cont,file = "/home/nmm199/R/MB_RNAseq/Clinical/clin.results/Oct_25_2018/Complete_transcripts_filtered/results.sig.multivar.cox.PFS.PNET.G3G4.cont.csv")

write.csv (significant.multivar.cox.OS.lancetG3G4.cont, file =  "/home/nmm199/R/MB_RNAseq/Clinical/clin.results/Oct_25_2018/Complete_transcripts_filtered/results.sig.multivar.cox.OS.lancetG3G4.cont.csv"  )  

write.csv(significant.multivar.cox.PFS.PNET5.cat, file = "/home/nmm199/R/MB_RNAseq/Clinical/clin.results/Oct_25_2018/Complete_transcripts_filtered/results.sig.multivar.cox.PFS.PNET5.cat.csv")



### annotated dataframes for significant multivar cox  OS and PFS

try(annot.sig.multi.cox.OS.PNET5.cont <- annotate.HTseq.IDs(rownames(significant.multivar.cox.OS.PNET5.cont)),silent = T) ### worked when updated annotate.HTseq.IDs function 6/11/18

try(annot.sig.multi.cox.OS.combined.cat <- annotate.HTseq.IDs(rownames(significant.multivar.cox.OS.combined.cat)),silent = T)

try(annot.sig.multi.cox.OS.combined.cont <- annotate.HTseq.IDs(rownames(significant.multivar.cox.OS.combined.cont)), silent =T) 

try(annot.sig.multi.cox.OS.lancetG3G4.cont <- annotate.HTseq.IDs(rownames(significant.multivar.cox.OS.lancetG3G4.cont)), silent = T) 

try(annot.sig.multivar.cox.OS.SHHold.cont <- annotate.HTseq.IDs(rownames(significant.multivar.cox.OS.SHHold.cont)), silent = T) 

try(annot.sig.multi.cox.OS.PNET5.cont <- annotate.HTseq.IDs(rownames(significant.multivar.cox.OS.PNET5.cont)),silent = T)

try(annot.sig.multivar.cox.PFS.combined.cat <- annotate.HTseq.IDs(rownames(significant.multivar.cox.PFS.combined.cat)), silent = T) ### reran several times before worked 6/11/18 despite new function changes

annot.sig.multivar.cox.PFS.SHHold.cat <- annotate.HTseq.IDs(rownames(significant.multivar.cox.PFS.SHHold.cat))




### write csv


write.csv (annot.sig.multi.cox.OS.combined.cat, file = "/home/nmm199/R/MB_RNAseq/Clinical/clin.results/Oct_25_2018/Complete_transcripts_filtered/annot.sig.multi.cox.OS.combined.cat.csv") 

write.csv (annot.sig.multivar.cox.PFS.combined.cat, file = "/home/nmm199/R/MB_RNAseq/Clinical/clin.results/Oct_25_2018/Complete_transcripts_filtered/annot.sig.multi.cox.PFS.combined.cat.csv") 

# write.csv(annot.sig.multivar.cox.PFS.SHHold.cat, file = "/home/nmm199/R/MB_RNAseq/Clinical/clin.results/Oct_25_2018/Complete_transcripts_filtered/annot.sig.multivar.cox.PFS.SHHold.cat.csv")

# annot.sig.multivar.cox.PFS.lancetG3G4.cat <- annotate.HTseq.IDs(rownames(significant.multivar.cox.PFS.lancetG3G4.cat)) 

# annot.sig.multivar.cox.PFS.lancetG3G4.cont <- annotate.HTseq.IDs(rownames(significant.multivar.cox.PFS.lancetG3G4.cont)) 


############################################################################
############################################################################

### schema when integrating new subsetting into a function
### 1. determine relevant formulae e. extract.cox.OS, cox.dataframe
### 2. can use hardcoding first to determine the cut off subset required e.g x[[5]]<15, see example for extract.multivar.cox.PFS.comb.cat.pval below
### 3. then create function using x, subset.index; and keep this function within the current file to tweak until all subset indices are working, or generate new functions as required
### 4. then move new function into main "clinical_data_functions_extract_master.R" file
### 5. OR can name the index, such as $p.val rather than [[subset.index]] described as a number
  

### now have replaced subset x[[5]]<15 for all, used to have following tips:
### there are fewer candidates found when have more liberal NA rule e.g for x[[4]]<17 (liberal, fewer candidates, less NAs) compared to x[[4]]<6 (strict)


#########################################################
#########################################################

