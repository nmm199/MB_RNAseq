### subsetting guide

### In the clinical_data_extract_DW file we have converted individual script from clinical_data_extract into a generic function. 
### In order to test that we have not lost capability or are subsetting unusual results, this script below is used to generate the same dataframes for comparison

### avoids duplicating objects within the same workspace
### currently there is one generic function to extract cox pvalue, Z score, HR, U95CI and L95CI, however depending on the statistical test, there may be subsetting "out of bounds errors"
### such as for SHH.old group



###########################################################################################################

###############################################

### extract cox regression p value and Z scores into individual dataframes

#################################################

### Cox PFS overall
### cox.U95CI.PFS.HR.all : this is the prototype of the function that works to extract elements when there is no value for some elements in that position

#########################################################

### cox PFS p value, categorical biomarker, entire cohort

extract.coxpval.test<- function (x){
  return(ifelse(length(x[[4]])< 3, NA,
                ifelse(length(x[[4]][[3]])<6, NA, 
                       ifelse(length(x[[4]][[3]][[1]])<1, NA,
                              x[[4]][[3]][[1]]))))
}


cox.PFS.cat.pval.all.test <- lapply(results.master, extract.coxpval.test)

### cox Z score for categorical biomarker

extract.PFS.cat.cox.Zscore.test <- function(x){
  return(ifelse(length(x[[4]])<3, NA,
                ifelse(length(x[[4]][[3]])<6, NA, 
                       ifelse(length(x[[4]][[3]][[5]])<1, NA,
                              x[[4]][[3]][[5]]))))
}

cox.PFS.cat.Zscore.all.test <- lapply(results.master, extract.PFS.cat.cox.Zscore.test) ### where biomarker is categorical variable

### cox PFS HR, categorical

extract.cox.PFS.cat.HR.test <- function(x){
  ifelse(length(x[[4]])<3, NA,
         ifelse(length(x[[4]][[3]])<6, NA, 
                ifelse(length(x[[4]][[3]][[2]])<1, NA,
                       x[[4]][[3]][[2]])))
}

cox.PFS.cat.HR.all.test <- lapply (results.master, extract.cox.PFS.cat.HR.test)

### cox PFS L95CI

extract.cox.cat.L95CI.HR.test <- function(x){
  ifelse(length(x[[4]])<3, NA,  ### get error message if make this < 2
         ifelse(length(x[[4]][[3]])<6, NA,  ### original script has <4 for this parameter, updated. ### 12109 is where the error is occurring
                ifelse(length(x[[4]][[3]][[3]])<1, NA,
                       x[[4]][[3]][[3]])))
}

cox.L95CI.PFS.cat.HR.all.test <- lapply (results.master, extract.cox.cat.L95CI.HR.test)


### cox PFS U95CI

extract.function.test <- function(x){
  return(ifelse(length(x[[4]]) < 3, NA,
                ifelse(length(x[[4]][[3]])<6, NA, ### changed from <4, to < 6. May need to alter for other RDS input file
                       ifelse(length(x[[4]][[3]][[4]])<1, NA, ### updated
                              x[[4]][[3]][[4]]))))
}

cox.U95CI.PFS.cat.HR.all.test <- lapply(results.master, extract.function.test)

cox.PFS.cat.all.df.test <- cox.dataframe(pval = cox.PFS.cat.pval.all.test, Zscore = cox.PFS.cat.Zscore.all.test, HR = cox.PFS.cat.HR.all.test, L95CI = cox.L95CI.PFS.cat.HR.all.test, U95CI = cox.U95CI.PFS.cat.HR.all.test)

colnames(cox.PFS.cat.all.df.test) <- c("cox.PFS.pval.cat.all", "cox.PFS.adj.pval.cat.all", "cox.PFS.Zscore.cat.all", "cox.PFS.HR.cat.all", "cox.PFS.L95CI.cat.all", "cox.PFS.U95CI.cat.all")

significant.cox.PFS.cat.all.test <- cox.PFS.cat.all.df.test[which(cox.PFS.cat.all.df.test[, 2]<0.05),]


### what is the best way to annotate these genes within R. I can do data manipulation in excel if needed

# annotate.cox.PFS.all <- annotate.HTseq.IDs(significant.cox.PFS.all) ### did not work "row names were found from a short variable and have been discarded"
# names(annotate.cox.PFS.all)<- gsub ("T", "", names(significant.cox.PFS.all))

try(annotate.cox.PFS.cat.all <- annotate.HTseq.IDs(rownames(significant.cox.PFS.cat.all)), silent = T)

write.csv(significant.cox.PFS.cat.all.test, file = "/home/nmm199/R/MB_RNAseq/Clinical/clin.results/significant.cox.PFS.cat.complete.allgroups.csv")





### cox PFS for G3G4 for categorical variable
extract.coxpval.PFS.cat.G3G4<- function (x){
  return(ifelse(length(x[[4]])< 5, NA, ###  works with length < 3 above for PFS for all groups, I had to change it to 5 or 6 for it to work for G3G4 PFS
                ifelse(length(x[[4]][[5]])<6, NA, 
                       ifelse(length(x[[4]][[5]][[1]])<1, NA,
                              x[[4]][[5]][[1]]))))
}

cox.PFS.pval.cat.G3G4 <- lapply(results.master, extract.coxpval.PFS.cat.G3G4)

### HR

extract.HR.PFS.cat.G3G4<- function (x){
  return(ifelse(length(x[[4]])< 5, NA, 
                ifelse(length(x[[4]][[5]])<6, NA, 
                       ifelse(length(x[[4]][[5]][[2]])<1, NA,
                              x[[4]][[5]][[2]]))))
}

cox.PFS.HR.cat.G3G4 <- lapply (results.master, extract.HR.PFS.cat.G3G4)

### Z score

extract.Zscore.PFS.cat.G3G4 <- function(x){
 return (ifelse(length(x[[4]])<5, NA,
               ifelse(length(x[[4]][[5]])<6, NA,
                     ifelse(length(x[[4]][[5]][[5]])<1, NA, 
                           x[[4]][[5]][[5]]))))
}

cox.PFS.Zscore.cat.G3G4 <- lapply (results.master, extract.Zscore.PFS.cat.G3G4)

### 95CI

 extract.L95CI.PFS.cat.G3G4 <- function(x){
 return (ifelse(length(x[[4]])<5, NA,
     ifelse(length(x[[4]][[5]])<6, NA,
            ifelse(length(x[[4]][[5]][[3]])<1, NA, 
                  x[[4]][[5]][[3]]))))
 }

 extract.U95CI.PFS.cat.G3G4 <- function(x){
 return (ifelse(length(x[[4]])<5, NA,
   ifelse(length(x[[4]][[5]])<6, NA,
        ifelse(length(x[[4]][[5]][[4]])<1, NA, 
             x[[4]][[5]][[4]]))))
 }

 cox.L95CI.PFS.HR.cat.G3G4 <- lapply(results.master, extract.L95CI.PFS.cat.G3G4)
 cox.U95CI.PFS.HR.cat.G3G4 <- lapply(results.master, extract.U95CI.PFS.cat.G3G4)

### cox dataframe
 cox.PFS.cat.G3G4.df <- cox.dataframe (pval = cox.PFS.pval.cat.G3G4, Zscore = cox.PFS.Zscore.cat.G3G4, HR = cox.PFS.HR.cat.G3G4, L95CI = cox.L95CI.PFS.HR.cat.G3G4 , U95CI = cox.U95CI.PFS.HR.cat.G3G4)
 colnames(cox.PFS.cat.G3G4.df) <- c("cox.PFS.pval.cat.G3G4", "cox.PFS.adj.pval.cat.G3G4", "cox.PFS.Zscore.cat.G3G4", "cox.PFS.HR.cat.G3G4", "cox.PFS.L95CI.cat.G3G4", "cox.PFS.U95CI.cat.G3G4")

 significant.cox.PFS.cat.G3G4 <- cox.PFS.cat.G3G4.df[which(cox.PFS.cat.G3G4.df[, 2]<0.05),]

 write.csv(significant.cox.PFS.cat.G3G4, file = "/home/nmm199/R/MB_RNAseq/Clinical/clin.results/significant.cox.PFS.cat.complete.G3G4.csv")


 
 
 
 ### Cox PFS for G3G4
 ### p value for continuous biomarker
 
 extract.coxpval.PFS.cont.G3G4<- function (x){
   return(ifelse(length(x[[4]])< 5, NA, ###  works with length < 3 above for PFS for all groups, I had to change it to 5 or 6 for it to work for G3G4 PFS
                 ifelse(length(x[[4]][[6]])<6, NA, 
                        ifelse(length(x[[4]][[6]][[1]])<1, NA,
                               x[[4]][[6]][[1]]))))
 }
 
 cox.PFS.pval.cont.G3G4 <- lapply(results.master, extract.coxpval.PFS.cont.G3G4)
 
 ### HR
 
 extract.HR.PFS.cont.G3G4<- function (x){
   return(ifelse(length(x[[4]])< 5, NA, 
                 ifelse(length(x[[4]][[6]])<6, NA, 
                        ifelse(length(x[[4]][[6]][[2]])<1, NA,
                               x[[4]][[6]][[2]]))))
 }
 
 cox.PFS.HR.cont.G3G4 <- lapply (results.master, extract.HR.PFS.cont.G3G4)
 
 ### Z score
 
 extract.Zscore.PFS.cont.G3G4 <- function(x){
   return (ifelse(length(x[[4]])<5, NA,
                  ifelse(length(x[[4]][[6]])<6, NA,
                         ifelse(length(x[[4]][[6]][[5]])<1, NA, 
                                x[[4]][[6]][[5]]))))
 }
 
 cox.PFS.Zscore.cont.G3G4 <- lapply (results.master, extract.Zscore.PFS.cont.G3G4)
 
 ### 95CI
 
 extract.L95CI.PFS.cont.G3G4 <- function(x){
   return (ifelse(length(x[[4]])<5, NA,
                  ifelse(length(x[[4]][[6]])<6, NA,
                         ifelse(length(x[[4]][[6]][[3]])<1, NA, 
                                x[[4]][[6]][[3]]))))
 }
 
 extract.U95CI.PFS.cont.G3G4 <- function(x){
   return (ifelse(length(x[[4]])<5, NA,
                  ifelse(length(x[[4]][[6]])<6, NA,
                         ifelse(length(x[[4]][[6]][[4]])<1, NA, 
                                x[[4]][[6]][[4]]))))
 }
 
 cox.L95CI.PFS.HR.cont.G3G4 <- lapply(results.master, extract.L95CI.PFS.cont.G3G4)
 cox.U95CI.PFS.HR.cont.G3G4 <- lapply(results.master, extract.U95CI.PFS.cont.G3G4)
 
 ### cox dataframe
 cox.PFS.cont.G3G4.df <- cox.dataframe (pval = cox.PFS.pval.cont.G3G4, Zscore = cox.PFS.Zscore.cont.G3G4, HR = cox.PFS.HR.cont.G3G4, L95CI = cox.L95CI.PFS.HR.cont.G3G4 , U95CI = cox.U95CI.PFS.HR.cont.G3G4)
 colnames(cox.PFS.cont.G3G4.df) <- c("cox.PFS.pval.cont.G3G4", "cox.PFS.adj.pval.cont.G3G4", "cox.PFS.Zscore.cont.G3G4", "cox.PFS.HR.cont.G3G4", "cox.PFS.L95CI.cont.G3G4", "cox.PFS.U95CI.cont.G3G4")
 
 significant.cox.PFS.cont.G3G4 <- cox.PFS.cont.G3G4.df[which(cox.PFS.cont.G3G4.df[, 2]<0.05),]
 
 write.csv(significant.cox.PFS.cont.G3G4, file = "/home/nmm199/R/MB_RNAseq/Clinical/clin.results/significant.cox.PFS.cont.complete.G3G4.csv")
 
 ### p value for categorical biomarker
 
 extract.coxpval.PFS.cat.SHH<- function (x){
   return(ifelse(length(x[[4]])< 3, NA, ###  works with length < 3 above for PFS for all groups, I had to change it to 5 or 6 for it to work for G3G4 PFS
                 ifelse(length(x[[4]][[7]])<6, NA, 
                        ifelse(length(x[[4]][[7]][[1]])<1, NA,
                               x[[4]][[7]][[1]]))))
 }
 
 cox.PFS.pval.cat.SHH <- lapply(results.master, extract.coxpval.PFS.cat.SHH)
 
 ### HR
 
 extract.HR.PFS.cat.SHH<- function (x){
   return(ifelse(length(x[[4]])< 3, NA, 
                 ifelse(length(x[[4]][[7]])<6, NA, 
                        ifelse(length(x[[4]][[7]][[2]])<1, NA,
                               x[[4]][[7]][[2]]))))
 }
 
 cox.PFS.HR.cat.SHH <- lapply (results.master, extract.HR.PFS.cat.SHH)
 
 ### Z score
 
 extract.Zscore.PFS.cat.SHH <- function(x){
   return (ifelse(length(x[[4]])<5, NA,
                  ifelse(length(x[[4]][[7]])<6, NA,
                         ifelse(length(x[[4]][[7]][[5]])<1, NA, 
                                x[[4]][[7]][[5]]))))
 }
 
 cox.PFS.Zscore.cat.SHH <- lapply (results.master, extract.Zscore.PFS.cat.SHH)
 
 ### 95CI
 
 extract.L95CI.PFS.cat.SHH <- function(x){
   return (ifelse(length(x[[4]])<5, NA,
                  ifelse(length(x[[4]][[7]])<6, NA,
                         ifelse(length(x[[4]][[7]][[3]])<1, NA, 
                                x[[4]][[7]][[3]]))))
 }
 
 extract.U95CI.PFS.cat.SHH <- function(x){
   return (ifelse(length(x[[4]])<5, NA,
                  ifelse(length(x[[4]][[7]])<6, NA,
                         ifelse(length(x[[4]][[7]][[4]])<1, NA, 
                                x[[4]][[7]][[4]]))))
 }
 
 cox.L95CI.PFS.HR.cat.SHH <- lapply(results.master, extract.L95CI.PFS.cat.SHH)
 cox.U95CI.PFS.HR.cat.SHH <- lapply(results.master, extract.U95CI.PFS.cat.SHH)
 
 ### cox dataframe
 cox.PFS.cat.SHH.df <- cox.dataframe (pval = cox.PFS.pval.cat.SHH, Zscore = cox.PFS.Zscore.cat.SHH, HR = cox.PFS.HR.cat.SHH, L95CI = cox.L95CI.PFS.HR.cat.SHH , U95CI = cox.U95CI.PFS.HR.cat.SHH)
 colnames(cox.PFS.cat.SHH.df) <- c("cox.PFS.pval.cat.SHH", "cox.PFS.adj.pval.cat.SHH", "cox.PFS.Zscore.cat.SHH", "cox.PFS.HR.cat.SHH", "cox.PFS.L95CI.cat.SHH", "cox.PFS.U95CI.cat.SHH")
 
 significant.cox.PFS.cat.SHH <- cox.PFS.cat.SHH.df[which(cox.PFS.cat.SHH.df[, 2]<0.05),]
 
 write.csv(significant.cox.PFS.cat.SHH, file = "/home/nmm199/R/MB_RNAseq/Clinical/clin.results/significant.cox.PFS.cat.complete.SHH.csv")
 
 #########################################################################################
 
 ### cox PFS (relapse) for SHH, continuous variable
 
 
 extract.coxpval.PFS.cont.SHH<- function (x){
   return(ifelse(length(x[[4]])< 3, NA, ###  works with length < 3 above for PFS for all groups, I had to change it to 5 or 6 for it to work for G3G4 PFS
                 ifelse(length(x[[4]][[8]])<6, NA, 
                        ifelse(length(x[[4]][[8]][[1]])<1, NA,
                               x[[4]][[8]][[1]]))))
 }
 
 cox.PFS.pval.cont.SHH <- lapply(results.master, extract.coxpval.PFS.cont.SHH)
 
 ### HR
 
 extract.HR.PFS.cont.SHH<- function (x){
   return(ifelse(length(x[[4]])< 3, NA, 
                 ifelse(length(x[[4]][[8]])<6, NA, 
                        ifelse(length(x[[4]][[8]][[2]])<1, NA,
                               x[[4]][[8]][[2]]))))
 }
 
 cox.PFS.HR.cont.SHH <- lapply (results.master, extract.HR.PFS.cont.SHH)
 
 ### Z score
 
 extract.Zscore.PFS.cont.SHH <- function(x){
   return (ifelse(length(x[[4]])<5, NA,
                  ifelse(length(x[[4]][[8]])<6, NA,
                         ifelse(length(x[[4]][[8]][[5]])<1, NA, 
                                x[[4]][[8]][[5]]))))
 }
 
 cox.PFS.Zscore.cont.SHH <- lapply (results.master, extract.Zscore.PFS.cont.SHH)
 
 ### 95CI
 
 extract.L95CI.PFS.cont.SHH <- function(x){
   return (ifelse(length(x[[4]])<5, NA,
                  ifelse(length(x[[4]][[8]])<6, NA,
                         ifelse(length(x[[4]][[8]][[3]])<1, NA, 
                                x[[4]][[8]][[3]]))))
 }
 
 extract.U95CI.PFS.cont.SHH <- function(x){
   return (ifelse(length(x[[4]])<5, NA,
                  ifelse(length(x[[4]][[8]])<6, NA,
                         ifelse(length(x[[4]][[8]][[4]])<1, NA, 
                                x[[4]][[8]][[4]]))))
 }
 
 cox.L95CI.PFS.HR.cont.SHH <- lapply(results.master, extract.L95CI.PFS.cont.SHH)
 cox.U95CI.PFS.HR.cont.SHH <- lapply(results.master, extract.U95CI.PFS.cont.SHH)
 
 ### cox dataframe
 cox.PFS.cont.SHH.df <- cox.dataframe (pval = cox.PFS.pval.cont.SHH, Zscore = cox.PFS.Zscore.cont.SHH, HR = cox.PFS.HR.cont.SHH, L95CI = cox.L95CI.PFS.HR.cont.SHH , U95CI = cox.U95CI.PFS.HR.cont.SHH)
 colnames(cox.PFS.cont.SHH.df) <- c("cox.PFS.pval.cont.SHH", "cox.PFS.adj.pval.cont.SHH", "cox.PFS.Zscore.cont.SHH", "cox.PFS.HR.cont.SHH", "cox.PFS.L95CI.cont.SHH", "cox.PFS.U95CI.cont.SHH")
 
 significant.cox.PFS.cont.SHH <- cox.PFS.cont.SHH.df[which(cox.PFS.cont.SHH.df[, 2]<0.05),]
 
 write.csv(significant.cox.PFS.cont.SHH, file = "/home/nmm199/R/MB_RNAseq/Clinical/clin.results/significant.cox.PFS.cont.complete.SHH.csv")
 
 
 
 ############################################################################
 
 
 
 ### cox PFS (relapse) for SHH.old, categorical 
 ### 11/10/17, tried to create a separate "extract.cox.SHH.old" function, pulling on individual aspects of the cox regression, worked once [[subset.index]] was included in function
 ### original script is below
 
 ### p value for categorical biomarker
 
 extract.coxpval.PFS.cat.SHH.old<- function (x){
   return(ifelse(length(x[[4]])< 10, NA, ###  works with length < 3 above for PFS for all groups, I had to change it to 10 for it to work for G3G4 PFS
                 ifelse(length(x[[4]][[9]])<6, NA, 
                        ifelse(length(x[[4]][[9]][[1]])<1, NA,
                               x[[4]][[9]][[1]]))))
 }
 
 cox.PFS.pval.cat.SHH.old <- lapply(results.master, extract.coxpval.PFS.cat.SHH.old)
 
 ### HR
 
 extract.HR.PFS.cat.SHH.old <- function (x){
   return(ifelse(length(x[[4]])< 10, NA, 
                 ifelse(length(x[[4]][[9]])<6, NA, 
                        ifelse(length(x[[4]][[9]][[2]])<1, NA,
                               x[[4]][[9]][[2]]))))
 }
 
 cox.PFS.HR.cat.SHH.old <- lapply (results.master, extract.HR.PFS.cat.SHH.old)
 
 ### Z score
 
 extract.Zscore.PFS.cat.SHH.old <- function(x){
   return (ifelse(length(x[[4]])<10, NA,
                  ifelse(length(x[[4]][[9]])<6, NA,
                         ifelse(length(x[[4]][[9]][[5]])<1, NA, 
                                x[[4]][[9]][[5]]))))
 }
 
 cox.PFS.Zscore.cat.SHH.old <- lapply (results.master, extract.Zscore.PFS.cat.SHH.old)
 
 ### 95CI
 
 extract.L95CI.PFS.cat.SHH.old <- function(x){
   return (ifelse(length(x[[4]])<10, NA,
                  ifelse(length(x[[4]][[9]])<6, NA,
                         ifelse(length(x[[4]][[9]][[3]])<1, NA, 
                                x[[4]][[9]][[3]]))))
 }
 
 extract.U95CI.PFS.cat.SHH.old <- function(x){
   return (ifelse(length(x[[4]])<10, NA,
                  ifelse(length(x[[4]][[9]])<6, NA,
                         ifelse(length(x[[4]][[9]][[4]])<1, NA, 
                                x[[4]][[9]][[4]]))))
 }
 
 cox.L95CI.PFS.HR.cat.SHH.old <- lapply(results.master, extract.L95CI.PFS.cat.SHH.old)
 cox.U95CI.PFS.HR.cat.SHH.old <- lapply(results.master, extract.U95CI.PFS.cat.SHH.old)
 
 
 ### cox PFS (relapse) for SHH.old, continuous variable
 
 
 extract.coxpval.PFS.cont.SHH.old <- function (x){
   return(ifelse(length(x[[4]])< 3, NA, ###  works with length < 3 above for PFS for all groups, I had to change it to 5 or 6 for it to work for G3G4 PFS
                 ifelse(length(x[[4]][[10]])<6, NA, 
                        ifelse(length(x[[4]][[10]][[1]])<1, NA,
                               x[[4]][[10]][[1]]))))
 }
 
 cox.PFS.pval.cont.SHH.old <- lapply(results.master, extract.coxpval.PFS.cont.SHH.old)
 
 ### HR
 
 extract.HR.PFS.cont.SHH.old <- function (x){
   return(ifelse(length(x[[4]])< 3, NA, 
                 ifelse(length(x[[4]][[10]])<6, NA, 
                        ifelse(length(x[[4]][[10]][[2]])<1, NA,
                               x[[4]][[10]][[2]]))))
 }
 
 cox.PFS.HR.cont.SHH.old <- lapply (results.master, extract.HR.PFS.cont.SHH.old)
 
 ### Z score
 
 extract.Zscore.PFS.cont.SHH.old <- function(x){
   return (ifelse(length(x[[4]])<5, NA,
                  ifelse(length(x[[4]][[10]])<6, NA,
                         ifelse(length(x[[4]][[10]][[5]])<1, NA, 
                                x[[4]][[10]][[5]]))))
 }
 
 cox.PFS.Zscore.cont.SHH.old <- lapply (results.master, extract.Zscore.PFS.cont.SHH.old)
 
 ### 95CI
 
 extract.L95CI.PFS.cont.SHH.old <- function(x){
   return (ifelse(length(x[[4]])<5, NA,
                  ifelse(length(x[[4]][[10]])<6, NA,
                         ifelse(length(x[[4]][[10]][[3]])<1, NA, 
                                x[[4]][[10]][[3]]))))
 }
 
 extract.U95CI.PFS.cont.SHH.old <- function(x){
   return (ifelse(length(x[[4]])<5, NA,
                  ifelse(length(x[[4]][[10]])<6, NA,
                         ifelse(length(x[[4]][[10]][[4]])<1, NA, 
                                x[[4]][[10]][[4]]))))
 }
 
 cox.L95CI.PFS.HR.cont.SHH.old <- lapply(results.master, extract.L95CI.PFS.cont.SHH.old)
 cox.U95CI.PFS.HR.cont.SHH.old <- lapply(results.master, extract.U95CI.PFS.cont.SHH.old)
 
 ### cox dataframe
 cox.PFS.cont.SHH.old.df <- cox.dataframe (pval = cox.PFS.pval.cont.SHH.old, Zscore = cox.PFS.Zscore.cont.SHH.old, HR = cox.PFS.HR.cont.SHH.old, L95CI = cox.L95CI.PFS.HR.cont.SHH.old , U95CI = cox.U95CI.PFS.HR.cont.SHH.old)
 colnames(cox.PFS.cont.SHH.old.df) <- c("cox.PFS.pval.cont.SHH.old", "cox.PFS.adj.pval.cont.SHH.old", "cox.PFS.Zscore.cont.SHH.old", "cox.PFS.HR.cont.SHH.old", "cox.PFS.L95CI.cont.SHH.old", "cox.PFS.U95CI.cont.SHH.old")
 
 
 
 ##################################
 
 ### p value
 
 extract.coxpval.OS.cat.all<- function (x){
   return(ifelse(length(x[[4]])< 11, NA, ###  works with length < 3 above for PFS for all groups, I had to change it to 12 on11/10/17
                 ifelse(length(x[[4]][[11]])<6, NA, 
                        ifelse(length(x[[4]][[11]][[1]])<1, NA,
                               x[[4]][[11]][[1]]))))
 }
 
 cox.OS.pval.cat.all <- lapply(results.master, extract.coxpval.OS.cat.all)
 
 ### HR
 
 extract.HR.OS.cat.all<- function (x){
   return(ifelse(length(x[[4]])< 11, NA, ###  works with length < 3 above for PFS for all groups, I had to change it to 5 or 6 for it to work for G3G4 PFS
                 ifelse(length(x[[4]][[11]])<6, NA, 
                        ifelse(length(x[[4]][[11]][[2]])<1, NA,
                               x[[4]][[11]][[2]]))))
 }
 
 cox.OS.HR.cat.all <- lapply(results.master, extract.HR.OS.cat.all)
 
 ### Z score
 
 extract.Zscore.OS.cat.all<- function (x){
   return(ifelse(length(x[[4]])< 11, NA, ###  works with length < 3 above for PFS for all groups, I had to change it to 5 or 6 for it to work for G3G4 PFS
                 ifelse(length(x[[4]][[11]])<6, NA, 
                        ifelse(length(x[[4]][[11]][[5]])<1, NA,
                               x[[4]][[11]][[5]]))))
 }
 
 cox.OS.Zscore.cat.all <- lapply(results.master, extract.Zscore.OS.cat.all)
 
 ### HR 95 CI
 
 extract.L95CI.OS.cat.all<- function (x){
   return(ifelse(length(x[[4]])< 11, NA, ###  works with length < 3 above for PFS for all groups, I had to change it to 5 or 6 for it to work for G3G4 PFS
                 ifelse(length(x[[4]][[11]])<6, NA, 
                        ifelse(length(x[[4]][[11]][[3]])<1, NA,
                               x[[4]][[11]][[3]]))))
 }
 
 
 extract.U95CI.OS.cat.all<- function (x){
   return(ifelse(length(x[[4]])< 11, NA, ###  works with length < 3 above for PFS for all groups, I had to change it to 5 or 6 for it to work for G3G4 PFS
                 ifelse(length(x[[4]][[11]])<6, NA, 
                        ifelse(length(x[[4]][[11]][[4]])<1, NA,
                               x[[4]][[11]][[4]]))))
 }
 
 cox.L95CI.OS.HR.cat.all<- lapply(results.master, extract.L95CI.OS.cat.all)
 cox.U95CI.OS.HR.cat.all<- lapply(results.master, extract.U95CI.OS.cat.all)
 
 ### cox OS dataframe
 
 cox.OS.cat.all.df <- cox.dataframe(pval = cox.OS.pval.cat.all , Zscore = cox.OS.Zscore.cat.all, HR = cox.OS.HR.cat.all, L95CI = cox.L95CI.OS.HR.cat.all, U95CI = cox.U95CI.OS.HR.cat.all)
 colnames(cox.OS.cat.all.df) <- c("cox.OS.pval.cat.all", "cox.OS.adj.pval.cat.all", "cox.OS.Zscore.cat.all", "cox.OS.HR.cat.all","cox.L95CI.OS.HR.cat.all", "cox.U95CI.OS.HR.cat.all")
 
 
 ########################################
 
 ### cox OS for all groups, continuous variable
 
 ### p value
 
 extract.coxpval.OS.cont.all<- function (x){
   return(ifelse(length(x[[4]])< 12, NA, ###  works with length < 3 above for PFS for all groups, I had to change it to 5 or 6 for it to work for G3G4 PFS
                 ifelse(length(x[[4]][[12]])<6, NA, 
                        ifelse(length(x[[4]][[12]][[1]])<1, NA,
                               x[[4]][[12]][[1]]))))
 }
 
 cox.OS.pval.cont.all <- lapply(results.master, extract.coxpval.OS.cont.all)
 
 
 ### HR
 
 extract.HR.OS.cont.all<- function (x){
   return(ifelse(length(x[[4]])< 12, NA, ###  works with length < 3 above for PFS for all groups, I had to change it to 5 or 6 for it to work for G3G4 PFS
                 ifelse(length(x[[4]][[12]])<6, NA, 
                        ifelse(length(x[[4]][[12]][[2]])<1, NA,
                               x[[4]][[12]][[2]]))))
 }
 
 cox.OS.HR.cont.all <- lapply(results.master, extract.HR.OS.cont.all)
 
 ### Z score
 
 extract.Zscore.OS.cont.all<- function (x){
   return(ifelse(length(x[[4]])< 12, NA, ###  works with length < 3 above for PFS for all groups, I had to change it to 5 or 6 for it to work for G3G4 PFS
                 ifelse(length(x[[4]][[12]])<6, NA, 
                        ifelse(length(x[[4]][[12]][[5]])<1, NA,
                               x[[4]][[12]][[5]]))))
 }
 
 cox.OS.Zscore.cont.all <- lapply(results.master, extract.Zscore.OS.cont.all)
 
 ### HR 95 CI
 
 extract.L95CI.OS.cont.all<- function (x){
   return(ifelse(length(x[[4]])< 12, NA, ###  works with length < 3 above for PFS for all groups, I had to change it to 5 or 6 for it to work for G3G4 PFS
                 ifelse(length(x[[4]][[12]])<6, NA, 
                        ifelse(length(x[[4]][[12]][[3]])<1, NA,
                               x[[4]][[12]][[3]]))))
 }
 
 cox.L95CI.OS.HR.cont.all<- lapply(results.master, extract.L95CI.OS.cont.all)
 
 
 extract.U95CI.OS.cont.all<- function (x){
   return(ifelse(length(x[[4]])< 12, NA, ###  works with length < 3 above for PFS for all groups, I had to change it to 5 or 6 for it to work for G3G4 PFS
                 ifelse(length(x[[4]][[12]])<6, NA, 
                        ifelse(length(x[[4]][[12]][[4]])<1, NA,
                               x[[4]][[12]][[4]]))))
 }
 
 cox.U95CI.OS.HR.cont.all<- lapply(results.master, extract.U95CI.OS.cont.all)
 
 ###########################################
 
 ### cox OS dataframe
 
 cox.OS.all.cont.df <- cox.dataframe(pval = cox.OS.pval.cont.all, Zscore = cox.OS.Zscore.cont.all, HR = cox.OS.HR.cont.all, L95CI = cox.L95CI.OS.HR.cont.all, U95CI = cox.U95CI.OS.HR.cont.all)
 colnames(cox.OS.all.cont.df) <- c("cox.OS.pval.cont.all", "cox.OS.adj.pval.cont.all", "cox.OS.Zscore.cont.all", "cox.OS.HR.cont.all","cox.L95CI.OS.HR.cont.all", "cox.U95CI.OS.HR.cont.all")
 
 
 ###
 
 ### Cox OS for G3G4, categorical 
 
 extract.coxpval.OS.cat.G3G4<- function (x){
   return(ifelse(length(x[[4]])< 13, NA, ###  works with length < 3 above for PFS for all groups, I had to change it to 5 or 6 for it to work for G3G4 PFS
                 ifelse(length(x[[4]][[13]])<6, NA, 
                        ifelse(length(x[[4]][[13]][[1]])<1, NA,
                               x[[4]][[13]][[1]]))))
 }
 
 cox.OS.pval.cat.G3G4 <- lapply(results.master, extract.coxpval.OS.cat.G3G4)
 
 ### HR 
 extract.HR.OS.cat.G3G4 <- function (x){
   return(ifelse(length(x[[4]])< 13, NA, ###  works with length < 3 above for PFS for all groups, I had to change it to 5 or 6 for it to work for G3G4 PFS
                 ifelse(length(x[[4]][[13]])<6, NA, 
                        ifelse(length(x[[4]][[13]][[2]])<1, NA,
                               x[[4]][[13]][[2]]))))
 }
 
 cox.OS.HR.cat.G3G4 <- lapply(results.master, extract.HR.OS.cat.G3G4)
 
 ### Z score
 
 extract.Zscore.OS.cat.G3G4<- function (x){
   return(ifelse(length(x[[4]])< 13, NA, ###  works with length < 3 above for PFS for all groups, I had to change it to 5 or 6 for it to work for G3G4 PFS
                 ifelse(length(x[[4]][[13]])<6, NA, 
                        ifelse(length(x[[4]][[13]][[5]])<1, NA,
                               x[[4]][[13]][[5]]))))
 }
 
 cox.OS.Zscore.cat.G3G4 <- lapply(results.master, extract.Zscore.OS.cat.G3G4)
 
 ### 95 CI
 
 extract.L95CI.OS.cat.G3G4<- function (x){
   return(ifelse(length(x[[4]])< 13, NA, ###  works with length < 3 above for PFS for all groups, I had to change it to 5 or 6 for it to work for G3G4 PFS
                 ifelse(length(x[[4]][[13]])<6, NA, 
                        ifelse(length(x[[4]][[13]][[3]])<1, NA,
                               x[[4]][[13]][[3]]))))
 }
 
 cox.L95CI.OS.HR.cat.G3G4 <- lapply(results.master, extract.L95CI.OS.cat.G3G4)
 
 extract.U95CI.OS.cat.G3G4<- function (x){
   return(ifelse(length(x[[4]])< 13, NA, ###  works with length < 3 above for PFS for all groups, I had to change it to 5 or 6 for it to work for G3G4 PFS
                 ifelse(length(x[[4]][[13]])<6, NA, 
                        ifelse(length(x[[4]][[13]][[4]])<1, NA,
                               x[[4]][[13]][[4]]))))
 }
 
 cox.U95CI.OS.HR.cat.G3G4 <- lapply(results.master, extract.U95CI.OS.cat.G3G4)
 
 
 ### cox dataframe for OS G3G4
 cox.OS.cat.G3G4.df <- cox.dataframe(pval = cox.OS.pval.cat.G3G4, Zscore = cox.OS.Zscore.cat.G3G4, HR = cox.OS.HR.cat.G3G4, L95CI = cox.L95CI.OS.HR.cat.G3G4, U95CI = cox.U95CI.OS.HR.cat.G3G4)
 colnames(cox.OS.cat.G3G4.df) <- c("cox.OS.pval.cat.G3G4", "cox.OS.adj.pval.cat.G3G4", "cox.OS.Zscore.cat.G3G4", "cox.OS.HR.cat.G3G4","cox.L95CI.OS.HR.cat.G3G4", "cox.U95CI.OS.HR.cat.G3G4")
 
 ################################
 
 ### Cox OS for G3G4, continuous
 
 extract.coxpval.OS.cont.G3G4<- function (x){
   return(ifelse(length(x[[4]])< 14, NA, ###  works with length < 3 above for PFS for all groups, I had to change it to 5 or 6 for it to work for G3G4 PFS
                 ifelse(length(x[[4]][[14]])<6, NA, 
                        ifelse(length(x[[4]][[14]][[1]])<1, NA,
                               x[[4]][[14]][[1]]))))
 }
 
 cox.OS.pval.cont.G3G4 <- lapply(results.master, extract.coxpval.OS.cont.G3G4)
 
 ### HR 
 extract.HR.OS.cont.G3G4 <- function (x){
   return(ifelse(length(x[[4]])< 14, NA, ###  works with length < 3 above for PFS for all groups, I had to change it to 5 or 6 for it to work for G3G4 PFS
                 ifelse(length(x[[4]][[14]])<6, NA, 
                        ifelse(length(x[[4]][[14]][[2]])<1, NA,
                               x[[4]][[14]][[2]]))))
 }
 
 cox.OS.HR.cont.G3G4 <- lapply(results.master, extract.HR.OS.cont.G3G4)
 
 ### Z score
 
 extract.Zscore.OS.cont.G3G4<- function (x){
   return(ifelse(length(x[[4]])< 14, NA, ###  works with length < 3 above for PFS for all groups, I had to change it to 5 or 6 for it to work for G3G4 PFS
                 ifelse(length(x[[4]][[14]])<6, NA, 
                        ifelse(length(x[[4]][[14]][[5]])<1, NA,
                               x[[4]][[14]][[5]]))))
 }
 
 cox.OS.Zscore.cont.G3G4 <- lapply(results.master, extract.Zscore.OS.cont.G3G4)
 
 ### 95 CI
 
 extract.L95CI.OS.cont.G3G4<- function (x){
   return(ifelse(length(x[[4]])< 14, NA, ###  works with length < 3 above for PFS for all groups, I had to change it to 5 or 6 for it to work for G3G4 PFS
                 ifelse(length(x[[4]][[14]])<6, NA, 
                        ifelse(length(x[[4]][[14]][[3]])<1, NA,
                               x[[4]][[14]][[3]]))))
 }
 
 cox.L95CI.OS.HR.cont.G3G4 <- lapply(results.master, extract.L95CI.OS.cont.G3G4)
 
 extract.U95CI.OS.cont.G3G4<- function (x){
   return(ifelse(length(x[[4]])< 14, NA, ###  works with length < 3 above for PFS for all groups, I had to change it to 5 or 6 for it to work for G3G4 PFS
                 ifelse(length(x[[4]][[14]])<6, NA, 
                        ifelse(length(x[[4]][[14]][[4]])<1, NA,
                               x[[4]][[14]][[4]]))))
 }
 
 cox.U95CI.OS.HR.cont.G3G4 <- lapply(results.master, extract.U95CI.OS.cont.G3G4)
 
 
 ### cox dataframe for OS G3G4
 cox.OS.cont.G3G4.df <- cox.dataframe(pval = cox.OS.pval.cont.G3G4, Zscore = cox.OS.Zscore.cont.G3G4, HR = cox.OS.HR.cont.G3G4, L95CI = cox.L95CI.OS.HR.cont.G3G4, U95CI = cox.U95CI.OS.HR.cont.G3G4)
 colnames(cox.OS.cont.G3G4.df) <- c("cox.OS.pval.cont.G3G4", "cox.OS.adj.pval.cont.G3G4", "cox.OS.Zscore.cont.G3G4", "cox.OS.HR.cont.G3G4","cox.L95CI.OS.HR.cont.G3G4", "cox.U95CI.OS.HR.cont.G3G4")
 
 
 
 ####################################################
 ### cox EFS for categorical, all patients (in survival cohort)
 
 try(cox.EFS.pval.all <- lapply(results.master, function(x){return(x[[4]][[1]][[1]])}), silent = T)
 try (cox.EFS.Zscore.all <- lapply(results.master, function(x){return(x[[4]][[1]][[5]])}), silent = T)
 try(cox.EFS.HR.all <- lapply(results.master, function(x){return(x[[4]][[1]][[2]])}), silent = T)
 try(cox.U95CI.EFS.HR.all <- lapply(results.master, function(x){return(x[[4]][[1]][[4]])}), silent = T)  
 try(cox.L95CI.EFS.HR.all <- lapply(results.master, function(x){return(x[[4]][[1]][[3]])}), silent = T)
 
 try(cox.EFS.all.df <- cox.dataframe(pval = cox.EFS.pval.all, Zscore = cox.EFS.Zscore.all, HR = cox.EFS.HR.all, L95CI = cox.L95CI.EFS.HR.all, U95CI = cox.U95CI.EFS.HR.all), silent = T)
 try(colnames(cox.EFS.all.df) <- c("cox.EFS.pval.all", "cox.EFS.adj.pval.all", "cox.EFS.Zscore.all", "cox.EFS.HR.all", "cox.U95CI.EFS.HR.all", "cox.L95CI.EFS.HR.all"), silent = T)
 
 
 
 ### cox EFS for G3G4
 
 try(cox.EFS.pval.G3G4  <- lapply(results.master, function(x){return(x[[4]][[2]][[1]])}), silent = T)
 try(cox.EFS.Zscore.G3G4 <- lapply(results.master, function(x){return(x[[4]][[2]][[5]])}), silent = T)
 try(cox.EFS.HR.G3G4 <- lapply(results.master, function(x){return(x[[4]][[2]][[2]])}), silent = T)
 try(cox.U95CI.EFS.HR.G3G4 <- lapply(results.master, function(x){return(x[[4]][[2]][[4]])}), silent = T)  
 try(cox.L95CI.EFS.HR.G3G4 <- lapply(results.master, function(x){return(x[[4]][[2]][[3]])}), silent = T)
 
 try(cox.EFS.G3G4.df <- cox.dataframe(pval = cox.EFS.pval.G3G4, Zscore = cox.EFS.Zscore.G3G4, HR = cox.EFS.HR.G3G4, L95CI = cox.L95CI.EFS.HR.G3G4, U95CI = cox.U95CI.EFS.HR.G3G4), silent = T)
 try(colnames(cox.EFS.G3G4.df) <- c("cox.EFS.pval.G3G4", "cox.EFS.adj.pval.G3G4", "cox.EFS.Zscore.G3G4", "cox.EFS.HR.G3G4", "cox.U95CI.EFS.HR.G3G4", "cox.L95CI.EFS.HR.G3G4"), silent = T)
 