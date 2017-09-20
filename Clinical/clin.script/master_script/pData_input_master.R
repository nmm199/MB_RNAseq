#####################################################################################
### update your pData object



### Functions used

source(file = "/home/nmm199/R/MB_RNAseq/Clinical/clin.script/master_script/clinical_data_functions_master.R")

x.data <- "/home/nmm199/R/MB_RNAseq/Input_data/database270617_060917.csv" ### changed as previous file database270617_310817 had PFS_R with non-numerical values
cat ("reading in clinical database", sep ="\n")
### add in row names to original database file
pData <- read.csv(x.data, row.names = 1)

meth.data <- "/home/nmm199/R/MB_RNAseq/Input_data/all7subgroupCalls.csv"
meth7 <- read.csv(meth.data, header=TRUE, sep=",", quote="\"", dec=".", row.names=1)

cytogen.data <- "/home/nmm199/R/MB_RNAseq/Input_data/arm_calls_clean280617.txt"
cytogen <- read.table (cytogen.data, header=T, sep="\t", row.names = 1)

test.pData <- updatepData(pData, meth7, cytogen, pdf.file = "/home/nmm199/R/MB_RNAseq/Clinical/clin.results/temp.pdf", log.file = "/home/nmm199/R/MB_RNAseq/Clinical/clin.results/temp.log.txt")
save(test.pData, file = "/home/nmm199/R/MB_RNAseq/Clinical/test.pData")





################################################################################