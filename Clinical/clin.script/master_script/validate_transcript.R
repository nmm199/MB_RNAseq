
### Script to validate transcripts in Cavalli dataset
### Yura Grabovska and Marion Mateos July 17 2018

# source("https://bioconductor.org/biocLite.R") ### use http if https does not work
# biocLite() ### to install program


library(car)
library(stats)
library(survival)
library(BiocInstaller)

#biocLite(pkgs = "affy")
#biocLite("survival")
#biocLite(pkgs = "dplyr")

library(affy)
library(biomaRt)
library (ggplot2)
library(dplyr)

# install.packages("affy")

# sessionInfo()

require(affy)

### unsure if need to run parallel but may help 

# library(foreach)
# library(tictoc)
library(parallel)
library(doParallel)
registerDoParallel(16)

##expression array data
# eset <- readRDS("/home/yuri/eSet_for_Dan_GSE85218.rds")
eset <- readRDS ("/home/yuri/eSet_GSE85218_26Sep2018.rds") 

##obtain a dataframe of the phenotype factors from the paper # (Yura omitted the copy number data because it wasn’t straightforward to tabulate from the paper supplements) 

pd <- pData(eset)

eset.exprs <- exprs(eset) ## obtain the expression values as a matrix

rownames(eset.exprs) <- gsub("_at", "", rownames(eset.exprs)) ## gets rid of the "_at" after the ENGS IDs ### length = 763

head(eset.exprs)

# “/home/dan/mygit/consensus_scripts/taylor_data”

# taylor <- View("/home/dan/mygit/consensus_scripts/taylor_data")

## the newer versions of the R biomaRt package run a series of smaller queries which are less likely to timeout internally without error install.packages("/home/yuri/biomaRt_2.37.0.tar.gz", repos = NULL, type="source") ## devel level version

require(biomaRt)
useMart("ENSEMBL_MART_ENSEMBL") ### used to be useMart()
## load the 'mart' and tell the program we want to use Ensembl and specifically the human annotations 

listMarts ### sometimes temperamental and have to run the useMart or listMarts more than once, then rerun ensembl.mart line. It worked after I ran the listMarts command then listMarts() then reran ensembl.mart (18/9/18)
listMarts() ### added 18/9/18

# listMarts(mart = NULL, host="www.ensembl.org", path="/biomart/martservice",
         # port=80, includeHosts = FALSE, archive = FALSE, ssl.verifypeer = TRUE, 
         # ensemblRedirect = TRUE, verbose = FALSE)

# ?listMarts
ensembl.mart <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")

# ensembl.mart <- useMart("ENSEMBL_MART_ENSEMBL", dataset = "hsapiens_gene_ensembl")


## set up the query: attributes are all the things we want to pull and filters are the IDs we want to filter by, takes 5-10 minutes without doParallel
hugene.anno <- getBM(attributes = c("ensembl_gene_id",
"hgnc_symbol",
"external_gene_name",
"description",
"chromosome_name",
"start_position",
"end_position",
"strand",
"band",
"gene_biotype",
"go_id",
"name_1006",
"definition_1006",
"go_linkage_type"),
filters = "ensembl_gene_id",
values = rownames(eset.exprs),
mart = ensembl.mart,
verbose = FALSE) ## setting verbose to TRUE will throw out all the raw API interface calls

## collapse data to unique gene-level, with all information concatenated by a separator; it's a messy call but it works and is relatively fast
require(plyr)
hugene.anno2 <- ddply(hugene.anno, .(ensembl_gene_id), summarize,
                     external_gene_name_u = paste(unique(external_gene_name), collapse=","),
                     hgnc_symbol_u = paste(unique(hgnc_symbol), collapse=","),
                     description_u = paste(unique(description), collapse=","),
                     chromosome_name_u = paste(unique(chromosome_name), collapse=","),
                     start_position_u = paste(unique(start_position), collapse=","),
                     end_position_u = paste(unique(end_position), collapse=","),
                     strand_u = paste(unique(strand), collapse=","),
                     band_u = paste(unique(band), collapse=","),
                     gene_biotype_u = paste(unique(gene_biotype), collapse=","),
                     go_id_u = paste(unique(go_id), collapse=","),
                     name_1006_u = paste(unique(name_1006), collapse=","),
                     definition_1006_u = paste(unique(definition_1006), collapse=","),
                     go_linkage_type_u = paste(unique(go_linkage_type), collapse=","))

# write.csv(hugene.anno, file = "/home/nmm199/R/MB_RNAseq/Clinical/clin.results/April_13_2018/hugene.anno.20180717.csv")

### now to validate the transcripts in survival analysis
View(pData(eset))

eset$OS <- eset$OS_.years.

### goi of interest
### create categorical variable
### then run survival analysis

# goi <- "ENSG00000165304"
# goi.df <- eset.exprs[rowname = goi,] 
# goi.cat <- ifelse(goi >median(goi, na.rm = T), "high","low")

### can add in MYC and MYCN data here

eset_match <- read.csv(file = "/home/nmm199/R/MB_RNAseq/Clinical/clin.results/April_13_2018/Complete_transcripts_filtered/Validation/eset_master_20181015.csv", header = TRUE, sep = ",", quote = "\"", row.names = 1)
head(eset_match)

#rownames(eset_match)<- eset_match[,1]

names(eset_match)
eset_match$Sample_Name
head(eset_match)
# View(pData(eset))

# colnames(eset_match)
### eset_match$Sample_Name this is the MB_SubtypeStudy_number

### next step is to move MYC and MYCN columns from eset_match into the eset that uses pData
### then can create this as main dataframe from which to subset expression datasets/ match in expression data. 

eset$MYC <- eset_match$MYC  ### this is OK as the sample numbers are in same order
eset$MYCN <- eset_match$MYCN
eset$meth7 <- eset_match$meth7
eset$q13loss <- eset_match$q13loss_YN
identical(rownames(pData(eset)),rownames(eset_match)) ### TRUE indicates that rownames are identical therefore valid to use above column additions

View(pData(eset)) 


##############################################################################################################
##############################################################################################################
### look to validate the goi in the Cavalli dataset

source(file = "/home/nmm199/R/MB_RNAseq/Clinical/clin.script/master_script/Validate_function.R")

##########################################################################################################   

### running for individual goi
  
list.goi.MELK <- goi.validate(goi = "ENSG00000165304", data = eset )
list.goi.NQO2 <- goi.validate(goi =  "ENSG00000124588", data = eset )
# list.goi.LOC256880 <- goi.validate (goi = "ENSG00000245322", data = eset) ### error subscript out of bounds
list.goi.MRPS12 <- goi.validate (goi = "ENSG00000128626", data = eset)
list.goi.MLYCD <- goi.validate (goi = "ENSG00000103150", data = eset)
list.goi.AKTIP <- goi.validate(goi = "ENSG00000166971", data = eset)
list.goi.SLC25A34 <- goi.validate (goi ="ENSG00000162461", data = eset )
list.goi.POU6F1 <- goi.validate (goi = "ENSG00000184271", data = eset)
list.goi.CXXC4 <- goi.validate (goi ="ENSG00000168772", data = eset )
list.goi.CHEK1 <- goi.validate (goi = "ENSG00000149554", data = eset)
list.goi.CDK4 <- goi.validate (goi = "ENSG00000135446", data = eset)
list.goi.ROGDI <- goi.validate (goi ="ENSG00000067836", data = eset )
list.goi.ASB16 <- goi.validate (goi = "ENSG00000161664", data = eset)
list.goi.KIAA0895L <- goi.validate (goi = "ENSG00000196123", data = eset)
list.goi.SATB2 <- goi.validate (goi = "ENSG00000119042", data = eset)
# list.goi.LINC01544 <- goi.validate (goi = "ENSG00000260440", data = eset) ### subscript out of bounds
# list.goi.ENSG261534 <- goi.validate (goi = "ENSG00000261534", data = eset) ### subscript out of bounds

write.csv(list.goi.NQO2$cox_summary$summary_nogender, file = "/home/nmm199/R/MB_RNAseq/Clinical/clin.results/April_13_2018/Complete_transcripts_filtered/Validation/NQO2_cox_nogender.csv")
write.csv(list.goi.NQO2$cox_summary$summary_gender, file = "/home/nmm199/R/MB_RNAseq/Clinical/clin.results/April_13_2018/Complete_transcripts_filtered/Validation/NQO2_cox_gender.csv")

write.csv(list.goi.MELK[[4]], file = "/home/nmm199/R/MB_RNAseq/Clinical/clin.results/Oct_25_2018/Complete_transcripts_filtered/2018_10_31/Validation/MELK_cox_nogender.csv")
write.csv (list.goi.NQO2[[4]], file = "/home/nmm199/R/MB_RNAseq/Clinical/clin.results/Oct_25_2018/Complete_transcripts_filtered/2018_10_31/Validation/NQO2_cox_nogender.csv")
write.csv (list.goi.MRPS12[[4]], file = "/home/nmm199/R/MB_RNAseq/Clinical/clin.results/Oct_25_2018/Complete_transcripts_filtered/2018_10_31/Validation/MRPS12_cox_nogender.csv")
write.csv (list.goi.MLYCD[[4]], file = "/home/nmm199/R/MB_RNAseq/Clinical/clin.results/Oct_25_2018/Complete_transcripts_filtered/2018_10_31/Validation/MLYCD_cox_nogender.csv")
write.csv (list.goi.AKTIP[[4]], file = "/home/nmm199/R/MB_RNAseq/Clinical/clin.results/Oct_25_2018/Complete_transcripts_filtered/2018_10_31/Validation/AKTIP_cox_nogender.csv")
write.csv (list.goi.SLC25A34[[4]], file = "/home/nmm199/R/MB_RNAseq/Clinical/clin.results/Oct_25_2018/Complete_transcripts_filtered/2018_10_31/Validation/SLC25A34__cox_nogender.csv")
write.csv (list.goi.POU6F1[[4]], file = "/home/nmm199/R/MB_RNAseq/Clinical/clin.results/Oct_25_2018/Complete_transcripts_filtered/2018_10_31/Validation/POU6F1_cox_nogender.csv")
write.csv (list.goi.CXXC4[[4]], file = "/home/nmm199/R/MB_RNAseq/Clinical/clin.results/Oct_25_2018/Complete_transcripts_filtered/2018_10_31/Validation/CXXC4_cox_nogender.csv") 
write.csv (list.goi.CHEK1[[4]], file =  "/home/nmm199/R/MB_RNAseq/Clinical/clin.results/Oct_25_2018/Complete_transcripts_filtered/2018_10_31/Validation/CHEK1_cox_nogender.csv")
write.csv (list.goi.CDK4[[4]], file =  "/home/nmm199/R/MB_RNAseq/Clinical/clin.results/Oct_25_2018/Complete_transcripts_filtered/2018_10_31/Validation/CDK4_cox_nogender.csv")
write.csv (list.goi.ROGDI[[4]], file =  "/home/nmm199/R/MB_RNAseq/Clinical/clin.results/Oct_25_2018/Complete_transcripts_filtered/2018_10_31/Validation/ROGD1_cox_nogender.csv")
write.csv (list.goi.ASB16 [[4]], file = "/home/nmm199/R/MB_RNAseq/Clinical/clin.results/Oct_25_2018/Complete_transcripts_filtered/2018_10_31/Validation/ASB16_cox_nogender.csv")
write.csv (list.goi.KIAA0895L[[4]], file = "/home/nmm199/R/MB_RNAseq/Clinical/clin.results/Oct_25_2018/Complete_transcripts_filtered/2018_10_31/Validation/KIAA0895L_cox_nogender.csv")
write.csv (list.goi.SATB2[[4]], file = "/home/nmm199/R/MB_RNAseq/Clinical/clin.results/Oct_25_2018/Complete_transcripts_filtered/Validation/2018_10_31/SATB2_cox_nogender.csv")

### cox overall combined (PNET5 and lancet) results

write.csv(list.goi.MELK[[1]], file = "/home/nmm199/R/MB_RNAseq/Clinical/clin.results/Oct_25_2018/Complete_transcripts_filtered/2018_10_31/Validation/MELK_cox_overall.csv")
write.csv (list.goi.NQO2[[1]], file = "/home/nmm199/R/MB_RNAseq/Clinical/clin.results/Oct_25_2018/Complete_transcripts_filtered/2018_10_31/Validation/NQO2_cox_overall.csv")
write.csv (list.goi.MRPS12[[1]], file = "/home/nmm199/R/MB_RNAseq/Clinical/clin.results/Oct_25_2018/Complete_transcripts_filtered/2018_10_31/Validation/MRPS12_cox_overall.csv")
write.csv (list.goi.MLYCD[[1]], file = "/home/nmm199/R/MB_RNAseq/Clinical/clin.results/Oct_25_2018/Complete_transcripts_filtered/2018_10_31/Validation/MLYCD_cox_overall.csv")
write.csv (list.goi.AKTIP[[1]], file = "/home/nmm199/R/MB_RNAseq/Clinical/clin.results/Oct_25_2018/Complete_transcripts_filtered/2018_10_31/Validation/AKTIP_cox_overall.csv")
write.csv (list.goi.SLC25A34[[1]], file = "/home/nmm199/R/MB_RNAseq/Clinical/clin.results/Oct_25_2018/Complete_transcripts_filtered/2018_10_31/Validation/SLC25A34_cox_overall.csv")
write.csv (list.goi.POU6F1[[1]], file = "/home/nmm199/R/MB_RNAseq/Clinical/clin.results/Oct_25_2018/Complete_transcripts_filtered/2018_10_31/Validation/POU6F1_cox_overall.csv")
write.csv (list.goi.CXXC4[[1]], file = "/home/nmm199/R/MB_RNAseq/Clinical/clin.results/Oct_25_2018/Complete_transcripts_filtered/2018_10_31/Validation/CXXC4_cox_overall.csv") 
write.csv (list.goi.CHEK1[[1]], file =  "/home/nmm199/R/MB_RNAseq/Clinical/clin.results/Oct_25_2018/Complete_transcripts_filtered/2018_10_31/Validation/CHEK1_cox_overall.csv")
write.csv (list.goi.CDK4[[1]], file =  "/home/nmm199/R/MB_RNAseq/Clinical/clin.results/Oct_25_2018/Complete_transcripts_filtered/2018_10_31/Validation/CDK4_cox_overall.csv")
write.csv (list.goi.ROGDI[[1]], file =  "/home/nmm199/R/MB_RNAseq/Clinical/clin.results/Oct_25_2018/Complete_transcripts_filtered/2018_10_31/Validation/ROGD1_cox_overall.csv")
write.csv (list.goi.ASB16 [[1]], file = "/home/nmm199/R/MB_RNAseq/Clinical/clin.results/Oct_25_2018/Complete_transcripts_filtered/2018_10_31/Validation/ASB16_cox_overall.csv")
write.csv (list.goi.KIAA0895L[[1]], file = "/home/nmm199/R/MB_RNAseq/Clinical/clin.results/Oct_25_2018/Complete_transcripts_filtered/2018_10_31/Validation/KIAA0895L_cox_overall.csv")
write.csv (list.goi.SATB2[[1]], file = "/home/nmm199/R/MB_RNAseq/Clinical/clin.results/Oct_25_2018/Complete_transcripts_filtered/2018_10_31/Validation/SATB2_cox_overall.csv")

### saveRDS 

# saveRDS(list.goi.NQO2, file = "/home/nmm199/R/MB_RNAseq/Clinical/clin.results/April_13_2018/Complete_transcripts_filtered/Validation/NQ02.rds")
# saveRDS(list.goi.MELK, file = "/home/nmm199/R/MB_RNAseq/Clinical/clin.results/April_13_2018/Complete_transcripts_filtered/Validation/MELK.rds")

### write out specific .csv based on the files in the rds object

# NQO2 <- readRDS(file = "/home/nmm199/R/MB_RNAseq/Clinical/clin.results/April_13_2018/Complete_transcripts_filtered/Validation/NQ02.rds")
# MELK <- readRDS(file = "/home/nmm199/R/MB_RNAseq/Clinical/clin.results/April_13_2018/Complete_transcripts_filtered/Validation/MELK.rds")
# NQO2$pval


### using continuous data

list.goi.MELK.cont <- goi.validate.contin  (goi = "ENSG00000165304", data = eset )
list.goi.NQO2.cont <- goi.validate.contin(goi =  "ENSG00000124588", data = eset )
list.goi.MRPS12.cont <- goi.validate.contin (goi = "ENSG00000128626", data = eset)
list.goi.MLYCD.cont <- goi.validate.contin (goi = "ENSG00000103150", data = eset)
list.goi.AKTIP.cont <- goi.validate.contin (goi = "ENSG00000166971", data = eset)
list.goi.SLC25A34.cont <- goi.validate.contin (goi ="ENSG00000162461", data = eset )
list.goi.POU6F1.cont <- goi.validate.contin (goi = "ENSG00000184271", data = eset)
list.goi.CXXC4.cont <- goi.validate.contin (goi ="ENSG00000168772", data = eset )
list.goi.CHEK1.cont <- goi.validate.contin (goi = "ENSG00000149554", data = eset)
list.goi.CDK4.cont <- goi.validate.contin (goi = "ENSG00000135446", data = eset)
list.goi.ROGDI.cont <- goi.validate.contin (goi ="ENSG00000067836", data = eset )
list.goi.ASB16.cont <- goi.validate.contin (goi = "ENSG00000161664", data = eset)
list.goi.KIAA0895L.cont <- goi.validate.contin (goi = "ENSG00000196123", data = eset)
list.goi.SATB2.cont <- goi.validate.contin (goi = "ENSG00000119042", data = eset)


list.goi.results <- as.list(mget(ls(pattern = "list.goi")))

#########################################################
#########################################################

