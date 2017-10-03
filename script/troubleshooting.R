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
