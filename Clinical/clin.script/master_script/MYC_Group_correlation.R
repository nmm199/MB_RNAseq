### MYC_Group_correlation
### March 27th 2018

RNA.data <- "/home/dan/mygit/rna_seq_mb/paper/MB.vsd.txt" 
mb.vsd <- read.delim(RNA.data)
### up to here

source(file = "/home/nmm199/R/MB_RNAseq/Clinical/clin.script/master_script/clinical_data_functions_master.R")
load("/home/nmm199/R/MB_RNAseq/Clinical/test.pData")

###need to define goi and goi.vsd

# goi <- "ENSG00000136997" ### this will only generate MYC expression for each sample, you can use your list of 1000 if you are interested in this
goi <- "ENSG00000249859" ###  this is PVT1
goi.vsd <- as.numeric(mb.vsd[goi,])

### need to generate the matched dataframe
index <- match(names(mb.vsd), rownames(test.pData)) 
matched.test.pData <- test.pData[index[!is.na(index)],] 
is.vector(matched.test.pData)

#as.data.frame(matched.test.pData) ### added 070917 in attempt to avoid downstream" $ atomic in a vector" error
matched.goi.vsd <- goi.vsd[!is.na(index)] ### linear expression
# matched.goi.vsd.cat <- ifelse(matched.goi.vsd>median(goi.vsd, na.rm = T), "high","low") 

### restrict dataframe to Group3
G3 <- matched.test.pData$meth =="G3" 

G3.group <- matched.test.pData [G3, ]

### restrict dataframe to Group4
G4 <- matched.test.pData$meth =="G4"
G4.group <- matched.test.pData [G4, ]
### restrict dataframe to combined Group3&4

G3G4.df <- rbind(G3.group, G4.group)

# nrow(G3G4.df)

# index.incl <- match(names(goi.vsd), rownames(G3G4.match.df)) 
# matched.G3G4.incl.pData <- G3G4.match.df[index.incl[!is.na(index.incl)],]

### running the GuiltbyAssociation
### this is the function

guiltByAssociation <-function(data, associated.gene, cores = 10){
  library(foreach)
  library(tictoc)
  library(parallel)
  library(doParallel)
  registerDoParallel(cores)
  
  index <- match(names(data), names(associated.gene))
  matched.data <- data[,!is.na(index)]
  matched.associated.gene <- associated.gene[index[!is.na(index)]] 
  
  guiltAssociation <- function(x,y){
    return(c(cor = cor.test(x,y)$estimate,
             p.val = cor.test(x,y)$p.value))
  }
  
  res <- foreach(i = 1:nrow(data), .combine = rbind)%dopar%{
    x <- as.numeric(matched.data[i,])
    return(guiltAssociation(x,matched.associated.gene))
  }
  
  rownames(res) <- rownames(data)
  return(res)
}

### change mb.vsd to G3.group or G3G4.df

MYC <- as.numeric(mb.vsd["ENSG00000136997.15_1",])
names(MYC) <- colnames(mb.vsd)
guilt.res.MYC <- guiltByAssociation(mb.vsd, MYC)

### run for Group 3 - will need to create a matched goi.vsd for each of these G3.group and for the combined G3G4.df


### run for Group 3&4


### then find out which gene the top transcript is associated with 

annot <- annotate.HTseq.IDs(rownames(guilt.res.MYC)) ### note this worked when I loaded the annotate.HTseq.IDs function again (clinical_data_functions_master.R)
cbind(guilt.res.MYC, adj.p.val=p.adjust(guilt.res.MYC[,2], method = "BH"), annot) -> guilt.res.MYC
guilt.res.MYC[!is.na(guilt.res.MYC[,1]),] -> guilt.res.MYC
head(guilt.res.MYC[order(guilt.res.MYC[,1]),],20) ### get the first 20 associated transcripts 
tail(guilt.res.MYC[order(guilt.res.MYC[,1]),],20) ### get the last 20 associated transcripts



