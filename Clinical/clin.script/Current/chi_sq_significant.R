### Script for running of significant chi squared results after run biomarker assessment folder
### This has been removed for the moment due to the slowness of the biomarker_assess_09_08_2017.R
### this can be added back in later if required
### have prioritised chi_sq_relapse to remain within the main biomarker_assess_09_08_2017.R 

### Author: Dr Marion Mateos
### Date: August 15 2017

### Companion files to run beforehand: 
# /home/nmm199/R/MB_RNAseq/Clinical/clin.script/current/prepare_pData_table.R
# /home/nmm199/R/MB_RNAseq/Clinical/clin.script/current/biomarker_discovery_functions.R
# /home/nmm199/R/MB_RNAseq/Clinical/clin.script/current/Biomarker_assess_09_08_2017.R


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

