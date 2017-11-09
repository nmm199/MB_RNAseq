### some additional script that is helpful for troubleshooting 29/8/17##

grep(pattern="chr17_gl000204_random", novel.gtf)

head(novel.gtf, 15655) ### this shows what the file looks like from beginning of file, skipping, then including up until row 15655 

### getwd() if not sure of where files are being saved

### if getting errors within a function then enter the actual function and add in each variable as named variable, run each line of script separately to see where error lies.



#################################
# try (script, silent = T )
### if ongoing errors try this:

# tryCatch{ #### This catches the error and outputs it to screen but allows the program to continue running
# results.master <- foreach(i = 1:nrow(mb.vsd))%dopar%{
#  as.numeric(mb.vsd [i,]) -> x
# names(x) <- colnames(mb.vsd)
#  names(x) <- gsub("T","",names(mb.vsd)) ### check that this is correct
#  error=function(e){cat("ERROR :",conditionMessage(e), "\n")} 
#  return(clinPathAssess(test.pData,x)) 
# }
# }  #### prints the error message to screen


### 9/11/17 when determining where the error lies for a script

library(foreach)
extract.multivar.cox <- function(results.master, subset.index){
  cox.dataframe (pval <- unlist(foreach(i = 1:length(results.master))%do%{return(extract.multivar.cox.pval(results.master[[i]], subset.index = subset.index))}),
                 Zscore <- unlist(foreach(i = 1:length(results.master))%do%{return(extract.multivar.cox.Zscore(results.master[[i]], subset.index = subset.index))})      , 
                 Zscore = lapply(results.master, extract.multivar.cox.Zscore, subset.index = subset.index),
                 HR = lapply(results.master, extract.multivar.cox.HR,subset.index = subset.index ),
                 L95CI = lapply(results.master, extract.multivar.cox.L95CI, subset.index = subset.index),
                 U95CI = lapply(results.master, extract.multivar.cox.U95CI, subset.index = subset.index )
  )               
}

### can use unlist and foreach, %do% (or %dopar%)