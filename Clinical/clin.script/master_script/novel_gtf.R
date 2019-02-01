### Identification of novel transcripts in original gtf file 
### Author: Dr Marion Mateos
### Date: October 12 2018

library(foreach)
library(doParallel)
biocLite(pckgs = c("rtracklayer", "GenomicFeatures", "zoo", "BSgenome.Hsapiens.UCSC.hg19", "scatterplot3d"), suppressUpdates = FALSE, suppressAutoUpdate = FALSE)

library(rtracklayer)
library(GenomicFeatures)
library(zoo)
library(BSgenome.Hsapiens.UCSC.hg19)
library(scatterplot3d)
registerDoParallel(cores = 10)

### gtf files used for the data
novel.gtf <- import("/home/dan/novel.mb.merged.gtf", form = "gtf")

### find required novel transcripts to confirm chromosomal location
### XLOC_050167 for PFS (categorical and continuous expression)in G3G4 and OS (continuous) in G3G4
### XLOC_003973 for PFS (categorical) in G3G4
### XLOC_056325 for PFS (continuous) in G3G4
### XLOC_043858 for OS (continuous) in G3G4

# class(novel.gtf)
# which(novel.gtf@elementMetadata@listData$gene_id =="XLOC_050167")

str(novel.gtf)

XLOC_050167 <- grep ("XLOC_050167", novel.gtf@elementMetadata@listData$gene_id)
XLOC_050167_grep <- novel.gtf[29363:29388, ]
XLOC_050167_grep_df <- as.data.frame(XLOC_050167_grep)
write.csv(XLOC_050167_grep_df,  file  = "/home/nmm199/R/MB_RNAseq/Clinical/clin.results/April_13_2018/Novel_transcripts/XLOC_050167.csv")

XLOC_003973 <- grep ("XLOC_003973", novel.gtf@elementMetadata@listData$gene_id )
XLOC_003973_grep_df <- as.data.frame(novel.gtf[875:877, ])
write.csv(XLOC_003973_grep_df, file = "/home/nmm199/R/MB_RNAseq/Clinical/clin.results/April_13_2018/Novel_transcripts/XLOC_003973.csv")

XLOC_056325 <- grep ("XLOC_056325", novel.gtf@elementMetadata@listData$gene_id )
XLOC_056325_grep_df <- as.data.frame(novel.gtf[32997, ])
write.csv(XLOC_056325_grep_df, file = "/home/nmm199/R/MB_RNAseq/Clinical/clin.results/April_13_2018/Novel_transcripts/XLOC_056325.csv")

XLOC_043858 <- grep ("XLOC_043858", novel.gtf@elementMetadata@listData$gene_id )
XLOC_043858_grep_df <- as.data.frame(novel.gtf[25623:25624, ])
write.csv(XLOC_043858_grep_df, file = "/home/nmm199/R/MB_RNAseq/Clinical/clin.results/April_13_2018/Novel_transcripts/XLOC_043858.csv")
