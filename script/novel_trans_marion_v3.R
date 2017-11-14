# setwd="/home/nmm199/R/MB_RNAseq/RNA_classes/"

#### first time you run you need to install this package
#install.packages(c("foreach", "doParallel"))


library(foreach)
library(doParallel)

#### first time you run you need to install this package
source("https://bioconductor.org/biocLite.R")
#?biocLite
##for the first time installing bioCLite packages
#biocLite(c("rtracklayer","GenomicFeatures","zoo","BSgenome.Hsapiens.UCSC.hg19","scatterplot3d")) OR can use
#biocLite(pkgs=c("rtracklayer", "GenomicFeatures", "zoo", "BSgenome.Hsapiens.UCSC.hg19", "scatterplot3d"), suppressUpdates=FALSE, suppressAutoUpdate=FALSE)


library(rtracklayer)
library(GenomicFeatures)
library(zoo)
library(BSgenome.Hsapiens.UCSC.hg19)
library(scatterplot3d)

### means it sets up all parallel functions to use 10 cores
registerDoParallel(cores = 10)
### load in description of all novel transcripts
#novel.gtf.old <- import("/media/microarray/output_cufflinks/cuffmerge/merged_asm/unknown.merged.gtf", format ="gtf")
novel.gtf <- import ("/home/dan/novel.mb.merged.gtf", form = "gtf") ### this is the updated file August 2017


### you can subset like a matrix e.g.
### novel.gtf[1:3,]
### novel.gtf[novel.gtf$gene_id=="XLOC_000022",]

### 1/5/17: query for Dan regarding txdb script line below as there is an error 'transcripts$tx_strand' must be "+" or "-"
# txdb <- makeTxDbFromGFF("/media/microarray/output_cufflinks/cuffmerge/merged_asm/unknown.merged.gtf")
### load refgene and lincRNAs
### 1/5/17: query for Dan regarding txdb script as there are warning messages regarding several parameters not used in query, ignore
txdb <- makeTxDbFromGFF(file = "/home/dan/mygit/rna_seq_mb/lincRNAs.ucsc.hg19.gtf")

txdb.refgene <- makeTxDbFromUCSC(genome="hg19", tablename="refGene")

### filter out any refgene and lincRNAs which overlap
novel.gtf[overlapsAny(novel.gtf,exons(txdb),ignore.strand=FALSE, maxgap=1000, minoverlap=1L),] -> tcons.linc.gtf
novel.gtf[overlapsAny(novel.gtf,exons(txdb.refgene),ignore.strand=FALSE,maxgap=1000, minoverlap=1L),] -> tcons.refgene.gtf
unique(tcons.linc.gtf$transcript_id) -> tcons.linc.remove
unique(tcons.refgene.gtf$transcript_id) -> tcons.refgene.remove
c(tcons.linc.remove,tcons.refgene.remove) -> tcons.all.remove
### 1/5/17 confirm with Dan: think this means that those that are being removed are removed into filt.novel.gtf and that novel.gtf contains only novel transcripts now
novel.gtf[-which(novel.gtf$transcript_id%in%tcons.all.remove),] -> filt.novel.gtf

################################################################################################
### this section below deals with assigning levels to each transcript and then determining sequence strings, which are then saved as separate files that are run against the PFAM moieties ()
### need to run the section below that contains "levs" and "nonrandom.novel.gtf" even if already have PFAM file that is generated using the RNA files contained in /home/R/MB_RNAseq/RNA_classes/working_files/split_fasta

###############################################################################################

### assign levels to each transcript (isoform)
levels(as.factor(novel.gtf$transcript_id)) -> levs

### return the number of exons
ex.no <- foreach(i = levs)%dopar%{
  novel.gtf[novel.gtf$transcript_id==i,] -> temp.grange
  return(max(temp.grange$exon_number))
}
ex.no <- unlist(ex.no)

###  extract each sequence as a DNAString
### need to exclude those transcripts that do not have a defined position within GenomicRanges. This is done by grepping for *random within the required gtf file
random.gtf <- grep(pattern = "*random", novel.gtf)   

### can treat gtf like a dataframe, by removing lines that are indicated by the random.gtf file
nonrandom.novel.gtf <- novel.gtf[-random.gtf,]

########
levs.nonrandom <- levels(as.factor(nonrandom.novel.gtf$transcript_id))  ### have now replaced everything that was levs with levs.nonrandom 1/11/17

seqs <- foreach(i = levs.nonrandom)%dopar%{
  nonrandom.novel.gtf[nonrandom.novel.gtf$transcript_id==i,] -> temp.grange ### used to be novel.gtf, now replaced with nonrandom.novel.gtf 29/8/17
  temp.seqs <- getSeq(Hsapiens, temp.grange)
  unlist(temp.seqs) -> temp.seqs
  return(temp.seqs)
}

### make a DNAStringSet
seqs.set <- DNAStringSet(seqs) 
names(seqs.set) <- levs.nonrandom


########################################
### if you already have the PFAM data against the sequence strings, then do not run this section below (line 93 -130)

### write out various versions of the sequence strings
# old script: writeXStringSet(seqs.set, file="./temp.seqs.fasta", format = "fasta", width=80)
writeXStringSet(seqs.set, file = "/home/nmm199/R/MB_RNAseq/RNA_classes/working_files/temp.seqs.fasta", format = "fasta", width=80)
writeXStringSet(seqs.set[ex.no>1], file="/home/nmm199/R/MB_RNAseq/RNA_classes/working_files/temp.seqs.ex.no.fasta", format = "fasta", width=80) 
### note that ex.no here is created using the levs file (i.e without randoms removed)
### ex.no is regenerated later on , as "no.ex" using levs.nonrandom, line 383

#dir.create("/home/nmm199/R/MB_RNAseq/RNA_classes/working_files/split_fasta")
1:length(seqs.set) -> d
split(d, ceiling(seq_along(d)/100)) -> splits

for(i in 1:length(splits)){
writeXStringSet(seqs.set[splits[[i]]], file=paste0("/home/nmm199/R/MB_RNAseq/RNA_classes/working_files/split_fasta/temp.seqs.",names(splits)[i],".fasta"), format = "fasta", width=80)
}


### run this to break down the larger seqs.set file to run PFAM more easily 31/8/17, replaced by for loop above 31/8/17
#writeXStringSet(seqs.set[1:100], file="/home/nmm199/R/MB_RNAseq/RNA_classes/working_files/temp.seqs.100.fasta", format = "fasta", width=80)
#writeXStringSet(seqs.set[1:2000], file="/home/nmm199/R/MB_RNAseq/RNA_classes/working_files/temp.seqs.2000.fasta", format = "fasta", width=80)
#writeXStringSet(seqs.set[2001:2000], file="/home/nmm199/R/MB_RNAseq/RNA_classes/working_files/temp.seqs.2000.fasta", format = "fasta", width=80)
#writeXStringSet(seqs.set[1:1000], file="/home/nmm199/R/MB_RNAseq/RNA_classes/working_files/temp.seqs.1000.fasta", format = "fasta", width=80)
#writeXStringSet(seqs.set[1:1000], file="/home/nmm199/R/MB_RNAseq/RNA_classes/working_files/temp.seqs.1000.fasta", format = "fasta", width=80)
#writeXStringSet(seqs.set[1:1000], file="/home/nmm199/R/MB_RNAseq/RNA_classes/working_files/temp.seqs.1000.fasta", format = "fasta", width=80)
#writeXStringSet(seqs.set[1:1000], file="/home/nmm199/R/MB_RNAseq/RNA_classes/working_files/temp.seqs.1000.fasta", format = "fasta", width=80)

### could run this first to check it is working before running it overnight
#writeXStringSet(seqs.set[1:10], file="./temp.seqs.10.fasta", format = "fasta", width=80) 


##################################################################################
##################################################################################

#### NOW GO AND PERFORM PFAM SEARCH ON SERVER ###
### tips: make sure use script created within atom to ensure runs currently on server using screen session. see macbook   Desktop/Bioinformatics/CPAT/PFAM/pfamscan.sh
### have created PFAM file now (5/10/17) via using FMS cluster. See files on macbook for further information (under Desktop/Bioinformatics/CPAT/PFAM). 
### this file is then used in the trinotate annotation script section below
### running here from 5/10/17 onwards, on NICR compsvr using screen

### pfam_scan.pl -translate orf -cpu 10 –fasta home/nmm199/R/MB_RNAseq/RNA classes/working files/temp.seqs.100.fasta -dir ./dat -outfile /home/nmm199/R/MB_RNAseq/RNA classes/working files/outfile.pfam


##################################################################################
##################################################################################

### read back in pfam results
#pfam <- read.table(file="/home/dan/pfam/pfam.res.10.txt", header = TRUE)   ### this is not needed 30/8/17
pfam <- read.table(file="/home/nmm199/R/MB_RNAseq/RNA_classes/working_files/all.seqs.fasta.no.head.res", header = TRUE) ### filename updated 5/10/17
colnames(pfam) <- c("seq_id", "alignment_start", "alignment_end", "envelope_start", "envelope_end", "hmm__acc", "hmm_name", "type", "hmm_start", "hmmend", "hmmlength", "bitscore", "E-value", "significance", "clan", "strand", "nt_start", "nt_end")


##################################################################################
#### RUN CPAT ON MAC ###

### note on file input: export current FASTA file (temp.seqs.fasta) from R studio on Macbook into Macbook downloads, and then move to folder /Users/Marion/Desktop/Bioinformatics/CPAT on macbook. 
### then run terminal on Macbook
### see shell script (cpat.sh) files on Atom on macbook
### then reexport files from Macbook onto R studio and read back into script below (currently exported temp.seqs.results.txt from Nov 14 2017 run)


##################################################################################
### read back in CPAT results

cpat.results <- read.delim(file = "/home/nmm199/R/MB_RNAseq/RNA_classes/working_files/temp.seqs.results.txt")   

#### set probability limit as defined in Landscape of drawing some graphs to visualise results
prob.limit = 0.5242
scatterplot3d(cpat.results$ORF_size,cpat.results$Fickett_score,cpat.results$Hexamer_score,
              pch = 19, ylim = c(0,1),zlim=c(0,1), color = ifelse(cpat.results$coding_prob>prob.limit,"red","blue"))


cpat.results[cpat.results$coding_prob>prob.limit,] -> cod.cpat.results
cod.cpat.results[order(cod.cpat.results$coding_prob, decreasing = TRUE),] -> order.cod.cpat.results
write.table(order.cod.cpat.results, file= "/home/nmm199/R/MB_RNAseq/RNA_classes/order.cod.cpat.results.txt", sep = "\t")


#toi <- "TCONS_00318898"
#as.character(seqs.set[toi])
#novel.gtf[novel.gtf$transcript_id==toi,]


### load in phastCons46way bigwig files
list.files("/home/dan/cons/phastCons46way", pattern = "*.bigWig", full.names = TRUE) -> phast.cons.files
### remove spurious files "gl"
### grep, grepl, regexpr, gregexpr and regexec search for matches to argument pattern within each element of a character vector: they differ in the format of and amount of detail in the results.
phast.cons.files[-grep("gl",phast.cons.files)] -> phast.cons.files


### script below imports big wig files, takes overnight, run over several cores
phast.cons<- mclapply(phast.cons.files, function(i){return(import(i, format="bw"))}, mc.cores = 32)

### sub and gsub perform replacement of the first and all matches respectively.
gsub(".phastCons46way.bigWig","",basename(phast.cons.files)) -> names(phast.cons)

### extract phast scores into vector

phast.cons.scores <- foreach(i = 1:length(phast.cons))%do%{
  start(phast.cons[[i]]) -> temp.starts
  rep(NA,max(temp.starts)) -> temp.vector
  temp.vector[temp.starts] <- phast.cons[[i]]$score
  return(temp.vector)
}
names(phast.cons) -> names(phast.cons.scores)

### phast.cons takes a few hours with 32 cores

# saveRDS(phast.cons, "/home/nmm199/R/MB_RNAseq/RNA_classes/working_files/May_2017/phast.cons.rds") ### previous file
# saveRDS(phast.cons, "/home/nmm199/R/MB_RNAseq/RNA_classes/working_files/phast.cons.rds") ### takes a while to save as RDS object
# load("/home/nmm199/R/MB_RNAseq/RNA_classes/working_files/May_2017/phast.cons.rds")
# load("/home/nmm199/R/MB_RNAseq/RNA_classes/working_files/phast.cons.rds")

###  5/5/17 phast.window script obsolete, thus deleted

#### function to make a progress bar compatible with parallel processing

f <- function(n){
  pb <- txtProgressBar(min=1, max=n-1,style=3)
  count <- 0
  function(...) {
    count <<- count + length(list(...)) - 1
    setTxtProgressBar(pb,count)
    Sys.sleep(0.01)
    flush.console()
    c(...)
  }
}


ptm <- proc.time()
library(utils)
pb <- txtProgressBar(min = 1, max = length(levs.nonrandom), style=3)
registerDoParallel(cores = 6)

### calculate maximum phast score within 200bp rolling window, script takes a few hours

phast.window <- foreach(i = 1:length(levs.nonrandom), .combine=c)%dopar%{
  setTxtProgressBar(pb, i)
  nonrandom.novel.gtf[nonrandom.novel.gtf$transcript_id==levs.nonrandom[i],] -> temp.grange ### changed from novel.gtf to nonrandom.novel.gtf 5/10/17
  as.character(seqnames(temp.grange))[1] -> chromosome
  temp.index <- foreach(j = 1:length(temp.grange), .combine=c)%do%{
    return(start(temp.grange)[j]:end(temp.grange)[j])
  }
  
  phast.cons.scores[[chromosome]][temp.index] -> temp.score
  TS <- zoo(c(temp.score))
  if(length(TS)>200){
    return(max(rollapply(TS, width = 200, FUN = mean, align = "left", na.rm = TRUE)))
  }else{
    return(NA)
  }
}

close(pb)
proc.time() - ptm
names(phast.window) <- levs.nonrandom


############ error solving, task failed at 5527, NA/NaN argument

# phast.window <- foreach(i = 5526:5528, .combine=c)%dopar%{
 # setTxtProgressBar(pb, i)
 # nonrandom.novel.gtf[nonrandom.novel.gtf$transcript_id==levs.nonrandom[i],] -> temp.grange ### changed from novel.gtf to nonrandom.novel.gtf 5/10/17
 # as.character(seqnames(temp.grange))[1] -> chromosome
 # temp.index <- foreach(j = 1:length(temp.grange), .combine=c)%do%{
 #   return(start(temp.grange)[j]:end(temp.grange)[j])
 # }
  
 # phast.cons.scores[[chromosome]][temp.index] -> temp.score
 # TS <- zoo(c(temp.score))
 # if(length(TS)>200){
 #   return(max(rollapply(TS, width = 200, FUN = mean, align = "left", na.rm = TRUE)))
 # }else{
 #   return(NA)
 # }
#}

# close(pb)
# proc.time() - ptm
# names(phast.window) <- levs.nonrandom
################################
#### save phast window scores

save(phast.window,file = "/home/nmm199/R/MB_RNAseq/RNA_classes/working_files/phast.window.Rdata")

### remove large memory object
# rm(phast.cons.scores)
# gc()

### call gc() after a large object has been removed, as this may prompt R to return memory to the operating system
rm(phast.cons)
gc()


##################################################################################################
### list all bigwig phyloP scores
list.files("/home/dan/cons/phyloP46way", pattern = "*.bigWig", full.names = TRUE) -> phyloP.files

#### remove spurious chromosomes
phyloP.files[-grep("gl",phyloP.files)] -> phyloP.files
phyloP.files[-grep("hap",phyloP.files)] -> phyloP.files

library(foreach)
### read in bigwig files
phyloP.cons <- foreach(i = phyloP.files)%dopar%{  ### updated to %dopar%  November 9 2017
  return(import(i, format="bw"))
}
gsub(".phyloP46way.bigWig","",basename(phyloP.files)) -> names(phyloP.cons)

### turn bigwig into vectors
phyloP.scores <- foreach(i = 1:length(phyloP.cons))%do%{  
  start(phyloP.cons[[i]]) -> temp.starts
  rep(NA, max(temp.starts)) -> temp.vector
  temp.vector[temp.starts] <- phyloP.cons[[i]]$score
  return(temp.vector)
}

names(phyloP.cons) -> names(phyloP.scores)
### remove large memory object
rm(phyloP.cons)
gc()

### if running this script from lines 301, will need to regenerate lines 28-87, as require levs.nonrandom, nonrandom.novel.gtf
#### return fraction of conserved bases score greater than 2
ptm <- proc.time()
registerDoParallel(cores = 6)
phyloP.fraction <- foreach(i = levs.nonrandom, .combine=c)%dopar%{   ### have changed to levs.nonrandom from levs (2/11/17)
  nonrandom.novel.gtf[nonrandom.novel.gtf$transcript_id==i,] -> temp.grange  ### changed to nonrandom.novel.gft from novel.gtf 5/10/17
  as.character(seqnames(temp.grange))[1] -> chromosome
  temp.index <- foreach(j = 1:length(temp.grange), .combine=c)%do%{
    return(start(temp.grange)[j]:end(temp.grange)[j])
  }
  phyloP.scores[[chromosome]][temp.index] -> temp.score
  return(length(which(temp.score > 2))/length(temp.score))
}
proc.time() - ptm

names(phyloP.fraction) <- levs.nonrandom
saveRDS(phyloP.fraction, file="/home/nmm199/R/MB_RNAseq/RNA_classes/working_files/phyloP.fraction.Rdata")
rm(phyloP.scores)
gc()



####################################################
### loading in all of the data, note PhyloP data from home/nmm199/R

load(file = "/home/dan/mygit/rna_seq_mb/ALL.specific.novel.Rdata")
load(file = "/home/dan/mygit/rna_seq_mb/ALL.specific.novel.isoforms.Rdata")
phyloP.fraction <- readRDS(file = "/home/nmm199/R/MB_RNAseq/RNA_classes/working_files/phyloP.fraction.Rdata")
load(file = "/home/nmm199/R/MB_RNAseq/RNA_classes/working_files/phast.window.Rdata")
load(file = "/home/dan/mygit/rna_seq_mb/filtered.gene_exp.Rdata")
load(file = "/home/dan/mygit/rna_seq_mb/filtered.gene_isoforms.Rdata")



### if running script from line 345, loading in generated data, need to run lines 28-87, as require levs.nonrandom, nonrandom.novel.gtf below

#### asign the specific group to genes and isoforms
rep(c("GRP3","GRP4","SHH","WNT"), unlist(lapply(filtered.gene_exp[[1]],length))) -> gene.groups
rep(c("GRP3","GRP4","SHH","WNT"), unlist(lapply(filtered.isoforms_exp[[1]],length))) -> isoform.groups

names(gene.groups) <- ALL.specific.novel
names(isoform.groups) <- ALL.specific.novel.isoforms

# read.delim(file="/home/nmm199/R/MB_Data/temp.seqs.results.txt") -> cpat.results ### older file version pre 9/11/17

cpat.results <- read.delim(file = "/home/nmm199/R/MB_RNAseq/RNA_classes/working_files/temp.seqs.results.txt")  

registerDoParallel(cores = 16) ### was 6

### check length of transcripts
### use nonrandom.novel.gtf file instead of novel.gtf for the below script, therefore changed below 5/10/17
length.transcript <- foreach(i = 1:length(levs.nonrandom), .combine=c)%dopar%{
  nonrandom.novel.gtf[nonrandom.novel.gtf$transcript_id==levs.nonrandom[i],] -> temp.grange ### changed here 5/10/17
  temp.index <- foreach(j = 1:length(temp.grange), .combine=c)%do%{
    return(start(temp.grange)[j]:end(temp.grange)[j])
  }
  return(length(temp.index))
}

### calculate number of exons
no.ex <- foreach(i = 1:length(levs.nonrandom), .combine=c)%dopar%{
  nonrandom.novel.gtf[nonrandom.novel.gtf$transcript_id==levs.nonrandom[i],] -> temp.grange  ### changed here 5/10/17
  return(length(temp.grange$exon_number))
}


### return name of gene for each transcript
gene.id <- foreach(i = 1:length(levs.nonrandom), .combine=c)%dopar%{
  nonrandom.novel.gtf[nonrandom.novel.gtf$transcript_id==levs.nonrandom[i],] -> temp.grange
  return(temp.grange$gene_id[1])
}

### load Matt's Trinotate annotation (note default home folder for dan is home/dan/mygit/rna_seq_mb/)
### check: do I need this section  5/10/17 
trinotate.annot <- read.delim(file="/home/dan/mygit/rna_seq_mb/trinotate_annotation_report.tab") ### this is where you update the name of the PFAM results file to use 30/8/17

### check if it has a pfam entry
any.pfam <- foreach(i = 1:length(levs.nonrandom), .combine=c)%dopar%{
  trinotate.annot[trinotate.annot$transcript_id==levs.nonrandom[i],] -> temp.trin
  return(any(temp.trin$Pfam!="."))
}

### return which pfam it is
specific.pfam <- foreach(i = 1:length(levs.nonrandom), .combine=c)%dopar%{
  trinotate.annot[trinotate.annot$transcript_id==levs.nonrandom[i],] -> temp.trin
  return(paste(temp.trin$Pfam, collapse = ";"))
}

### makes some logical vectors describing whether pass cutoff
(no.ex > 1) -> no.ex.index
(length.transcript > 250) -> length.transcript.index
(phyloP.fraction>0.0947) -> conservation.index
(phast.window>0.9986) -> uce.index
cpat.results$coding_prob>0.5242 -> cpat.index
# levs.nonrandom%in%ALL.specific.novel.isoforms ### relevance?


#### put all into one dataframe
conservation.scores <- data.frame(length.transcript = length.transcript,
                                  phyloP.fraction =  phyloP.fraction,
                                  phast.window=phast.window,
                                  cpat.prob = cpat.results$coding_prob,
                                  no.ex = no.ex,
                                  length.transcript.index  = length.transcript.index,
                                  conservation.index = conservation.index,
                                  uce.index = uce.index,
                                  cpat.index = cpat.index,
                                  no.ex.index = no.ex.index,
                                  any.pfam = any.pfam,
                                  gene.id = gene.id,
                                  subgroup.spec = levs.nonrandom%in%ALL.specific.novel.isoforms,
                                  specific.pfam = specific.pfam,
                                  subgroup.isoforms = isoform.groups[levs.nonrandom],
                                  subgroup.genes = gene.groups[gene.id]
)


save (conservation.scores, file = "/home/nmm199/R/MB_RNAseq/RNA_classes/working_files/conservation.scores.Rdata") ### if using saveRDS, use with .rds

### CPAT results updated so that dataframe now works, with 36071 transcripts in dataframe. 

### remove overlaps with refgene/LINCRNAs
conservation.scores[-which(rownames(conservation.scores)%in%tcons.all.remove),] -> conservation.scores.remove

dim(conservation.scores.remove)
length(unique(conservation.scores.remove$gene.id))

conservation.scores.remove[conservation.scores.remove$subgroup.spec&(conservation.scores.remove$cpat.index|conservation.scores.remove$any.pfam),] -> subgroup.specific.cons

subgroup.specific.cons[order(subgroup.specific.cons$subgroup.isoforms),c(-14,-16)]


nonrandom.novel.gtf[nonrandom.novel.gtf$transcript_id%in%rownames(subgroup.specific.cons),] ### changed 5/10/17

# nonrandom.novel.gtf[nonrandom.novel.gtf$transcript_id=="TCONS_01254472",]

### only transcripts greater than 250bp
conservation.scores.remove[conservation.scores.remove$length.transcript.index,] -> conservation.scores.remove.length

#### None overlapping with GENOCDE/LINC/REFGENE
nrow(conservation.scores.remove.length) -> total.iso
length(unique(conservation.scores.remove.length$gene.id))  -> total.gen
c(total.iso,total.gen)


### Number of exons
conservation.scores.remove.length[which(conservation.scores.remove.length$no.ex.index),] -> conservation.scores.remove.length.no.ex
nrow(conservation.scores.remove.length.no.ex) -> total.multi.iso
length(unique(conservation.scores.remove.length.no.ex$gene.id)) -> total.multi.gen
total.iso - total.multi.iso -> total.single.iso
total.gen-total.multi.gen -> total.single.gen
### Multi-exonic
c(total.multi.iso,total.multi.gen)
### single-exonic
c(total.single.iso,total.single.gen )


### Coding Potential CPAT ### CPAT file updated 14/11/17, note that TCONS names are generated from cuffdiff-level, and thus same TCONS will be used in PFAM and CPAT files.

conservation.scores.remove.length.no.ex[which(conservation.scores.remove.length.no.ex$cpat.index),] -> cpat.table
nrow(cpat.table) -> total.cpat.iso
length(unique(cpat.table$gene.id)) -> total.cpat.gen
c(total.cpat.iso,total.cpat.gen)
### Coding Potential PFAM
conservation.scores.remove.length.no.ex[which(conservation.scores.remove.length.no.ex$any.pfam),] -> pfam.table
nrow(pfam.table) -> total.pfam.iso
length(unique(pfam.table$gene.id)) -> total.pfam.gen
c(total.pfam.iso,total.pfam.gen)
### Coding Potential TOTAL
conservation.scores.remove.length.no.ex[which(conservation.scores.remove.length.no.ex$cpat.index|conservation.scores.remove.length.no.ex$any.pfam),] -> tucp.table
conservation.scores.remove.length.no.ex[-which((conservation.scores.remove.length.no.ex$cpat.index|conservation.scores.remove.length.no.ex$any.pfam)),] -> linc.table
nrow(tucp.table) -> total.tucp.iso
length(unique(tucp.table$gene.id)) -> total.tucp.gen
c(total.tucp.iso,total.tucp.gen)

### LINC RNAs
total.multi.iso - total.tucp.iso -> total.linc.iso
total.multi.gen - total.tucp.gen -> total.linc.gen
c(total.linc.iso,total.linc.gen)

### subgroup specific TUCP
# dir.create("/home/nmm199/R/MB_RNAseq/RNA_classes/subgroup.tucp")
tucp.table[which(tucp.table$subgroup.spec),] -> subgroup.tucp.table
write.table(subgroup.tucp.table, file = "/home/nmm199/R/MB_RNAseq/RNA_classes/subgroup.tucp/subgroup.tucp.table.txt", sep = "\t")
nrow(subgroup.tucp.table) -> total.subgroup.tucp.iso
length(unique(subgroup.tucp.table$gene.id)) -> total.subgroup.tucp.gen
c(total.subgroup.tucp.iso,total.subgroup.tucp.gen)
### subgroup specific LINC ***** OUTPUT THIS ONE FOR JANET ****
linc.table[which(linc.table$subgroup.spec),] -> subgroup.linc.table
# dir.create("/home/nmm199/R/MB_RNAseq/RNA_classes/subgroup.linc")
write.table(subgroup.linc.table, file = "/home/nmm199/R/MB_RNAseq/RNA_classes/subgroup.linc/subgroup.linc.table.txt", sep = "\t")
nrow(subgroup.linc.table) -> total.subgroup.linc.iso
length(unique(subgroup.linc.table$gene.id)) -> total.subgroup.linc.gen
c(total.subgroup.linc.iso,total.subgroup.linc.gen)

#### conserved elements LINC UCE
subgroup.linc.table[which(subgroup.linc.table$uce.index),] -> uce.linc.table
nrow(uce.linc.table) -> total.uce.linc.iso
length(unique(uce.linc.table$gene.id)) -> total.uce.linc.gen
c(total.uce.linc.iso,total.uce.linc.gen)
#### conserved elements TUCP UCE
subgroup.tucp.table[which(subgroup.tucp.table$uce.index),] -> uce.tucp.table
nrow(uce.tucp.table) -> total.uce.tucp.iso
length(unique(uce.tucp.table$gene.id)) -> total.uce.tucp.gen
c(total.uce.tucp.iso,total.uce.tucp.gen)
#### conserved elements LINC CON
subgroup.linc.table[which(subgroup.linc.table$conservation.index),] -> conservation.linc.table
nrow(conservation.linc.table) -> total.conservation.linc.iso
length(unique(conservation.linc.table$gene.id)) -> total.conservation.linc.gen
c(total.conservation.linc.iso,total.conservation.linc.gen)
#### conserved elements TUCP CON
subgroup.tucp.table[which(subgroup.tucp.table$conservation.index),] -> conservation.tucp.table
nrow(conservation.tucp.table) -> total.conservation.tucp.iso
length(unique(conservation.tucp.table$gene.id)) -> total.conservation.tucp.gen
c(total.conservation.tucp.iso,total.conservation.tucp.gen)
#### conserved elements TUCP TOTAL CON
subgroup.tucp.table[which(subgroup.tucp.table$uce.index|subgroup.tucp.table$conservation.index),] -> cons.tucp.table
nrow(cons.tucp.table) -> total.cons.tucp.iso
length(unique(cons.tucp.table$gene.id)) -> total.cons.tucp.gen
c(total.cons.tucp.iso,total.cons.tucp.gen)
#### conserved elements LINC TOTAL CON
subgroup.linc.table[which(subgroup.linc.table$uce.index|subgroup.linc.table$conservation.index),] -> cons.linc.table
nrow(cons.linc.table) -> total.cons.linc.iso
length(unique(cons.linc.table$gene.id)) -> total.cons.linc.gen
c(total.cons.linc.iso,total.cons.linc.gen)


### put it all together in a table and output

dir.create("/home/nmm199/R/MB_RNAseq/RNA_classes/results")

nrow(conservation.scores.remove.length.no.ex)

rbind(
  paste(length(which(conservation.scores.remove.length.no.ex$cpat.index)),length(which(!conservation.scores.remove.length.no.ex$cpat.index)), sep = ","),
  paste(length(which(conservation.scores.remove.length.no.ex$conservation.index)),length(which(!conservation.scores.remove.length.no.ex$conservation.index)), sep = ","),
  paste(length(which(conservation.scores.remove.length.no.ex$uce.index)),length(which(!conservation.scores.remove.length.no.ex$uce.index)), sep = ","),
  paste(length(which(conservation.scores.remove.length.no.ex$any.pfam)),length(which(!conservation.scores.remove.length.no.ex$any.pfam)), sep = ","),
  paste(length(which(conservation.scores.remove.length.no.ex$subgroup.spec)),length(which(!conservation.scores.remove.length.no.ex$subgroup.spec)), sep = ",")
) -> total.counts
rownames(total.counts)<- c("cpat", "cons","UCE","pfam","subgroup")
colnames(total.counts)<- "total"

rbind(
  tapply(conservation.scores.remove.length.no.ex$cpat.index,conservation.scores.remove.length.no.ex$subgroup.isoforms,function(x){paste(length(which(x)),length(which(!x)), sep = ",")}),
  tapply(conservation.scores.remove.length.no.ex$conservation.index,conservation.scores.remove.length.no.ex$subgroup.isoforms,function(x){paste(length(which(x)),length(which(!x)), sep = ",")}),
  tapply(conservation.scores.remove.length.no.ex$uce.index,conservation.scores.remove.length.no.ex$subgroup.isoforms,function(x){paste(length(which(x)),length(which(!x)), sep = ",")}),
  tapply(conservation.scores.remove.length.no.ex$any.pfam,conservation.scores.remove.length.no.ex$subgroup.isoforms,function(x){paste(length(which(x)),length(which(!x)), sep = ",")}),
  tapply(conservation.scores.remove.length.no.ex$subgroup.spec,conservation.scores.remove.length.no.ex$subgroup.isoforms,function(x){paste(length(which(x)),length(which(!x)), sep = ",")})
) -> subgroup.counts
rownames(subgroup.counts)<- c("cpat", "cons","UCE","pfam","subgroup")

write.table(cbind(total.counts, subgroup.counts), file = "/home/nmm199/R/MB_RNAseq/RNA_classes/results/counts.table.novel.trans.txt", sep = "\t")

###########################################################################
### output as venn diagram
###########################################################################

table(conservation.scores.remove.length.no.ex[,c(7,8,9,11,13)])
# install.packages ('gplots')
# install.packages('VennDiagram')

library(gplots)
library('VennDiagram')

pdf(file = "/home/nmm199/R/MB_RNAseq/RNA_classes/results/venn.diagram.novel.trans.pdf")
v.table<-venn( list("conserved" = rownames(conservation.scores.remove.length.no.ex[conservation.scores.remove.length.no.ex[,7]==TRUE,]),
                    "UCE" = rownames(conservation.scores.remove.length.no.ex[which(conservation.scores.remove.length.no.ex[,8]==TRUE),]),
                    "CPAT" = rownames(conservation.scores.remove.length.no.ex[conservation.scores.remove.length.no.ex[,9]==TRUE,]),
                    "PFAM" = rownames(conservation.scores.remove.length.no.ex[conservation.scores.remove.length.no.ex[,11]==TRUE,]),
                    "subgroup specific" = rownames(conservation.scores.remove.length.no.ex[conservation.scores.remove.length.no.ex[,13]==TRUE,])
) )
dev.off()
str(v.table)
print(v.table)


########################################################################
### there may be errors with the CPAT file not aligning with the previously generated PFAM file etc, therefore to view fasta file:
### Nov 14 2017

library("Biostrings")
temp.seqs = readDNAStringSet("/home/nmm199/R/MB_RNAseq/RNA_classes/working_files/temp.seqs.fasta")

### other options that did not work:
# install.packages("seqinr", repos="http://R-Forge.R-project.org")
# read.fasta(file = system.file("sequences/ct.fasta.gz", package = "seqinr"), 
   #         seqtype = c("DNA", "AA"), as.string = FALSE, forceDNAtolower = TRUE,
   #         set.attributes = TRUE, legacy.mode = TRUE, seqonly = FALSE, strip.desc = FALSE,
    #        bfa = FALSE, sizeof.longlong = .Machine$sizeof.longlong,
     #       endian = .Platform$endian, apply.mask = TRUE)

# temp.seqs <- read.fasta(file = system.file ("/home/nmm199/R/MB_RNAseq/RNA_classes/working_files/temp.seqs.fasta", package = "seqinr"), 
 #         seqtype = c("DNA", "AA"), as.string = FALSE, forceDNAtolower = TRUE, 
 #         set.attributes = TRUE)

###

