#######################

### File information: this is unused script taken from clinical data script (clinical_data_3.R)
### File aim: 
### 1. Provide back up information and code in case this is required for further analysis
### 2. Files that were used to cross-reference RNA seq cohort (NMB numbers)
### 3. Files that were used to run the cytogenetic arm change analysis

### Author: Dr Marion Mateos
### Date: 3/7/17




### input files:
file = "/home/nmm199/R/MB_RNAseq/Input data/arm_calls_clean280617.txt"




### Functions:

### frequencies of cytogenetic arm changes

cytogen.freq <- function (x){
  table.temp <- table (x)
  return (table.temp)
}

# x <- chrosomomal arm
### returns frequency of gain, loss, neutral for each chromosomal arm included in Ed Schwalbe's arm analysis



########################################
### matching NMB numbers to determine which included in RNA seq for reporting purposes

# test.pData$RNAdata <- mb.vsd[match(rownames(test.pData), colnames(mb.vsd)), ]

### can save it as a vector instead
RNAseq <- colnames (mb.vsd)
RNAseq

### save as readable text file, which can then be uploaded to computer harddrive and matched to existing clinical database in excel
write.table(RNAseq, "/home/nmm199/R/MB_RNAseq/Clinical/RNAseq.txt", sep="\t")

#########################################


### Distribution of the data


### note that s.d. is only used with mean, with normally distributed symmetric data
# var(data)
# sd(age)
# range(age)
# head(age)
# min(age)
# max(age)

##################################################

### other ways of constructing chi squared table and chisq test

# table(as.factor(MYC.cat),as.factor(relapse))
# chisq.test(table(as.factor(MYC.cat),as.factor(relapse)))
# tab.MYC.relapse <- table(MYC.cat, relapse)



######################################################

### pearson's correlation coefficient
# cor.test(x, y)
# plot(x,y)
# abline (lm(x ~y)



### if relationship does not appear linear then use spearman correlation coefficient

# cor.test(x, y, method = 'spearman')


#########################################


### linear regression
# this may not be required, could use correlation above

age.reg <- lm(matched.test.pData$age.cont ~ matched.goi.vsd)
age.reg
summary(age.reg)



#########################################



### Mann-whitney U (aka wilcoxon rank sum test) for non-parametric data
### note that if the sample size is large enough (e.g > 30) then could theoretically use t-test, even if the data is non-parametric
### however we will use the Wilcoxon rank sum test



#######################################



### assessing the frequency for each arm changes (e.g p1 is short arm on chromosome 1)

cytogen.p1 <- cytogen.freq(x = cytogen$p1)
cytogen.p2 <- cytogen.freq(x = cytogen$p2)
cytogen.p3 <- cytogen.freq(x = cytogen$p3)
cytogen.p4 <- cytogen.freq(x= cytogen$p4)
cytogen.p5 <- cytogen.freq(x= cytogen$p5)
cytogen.p6 <- cytogen.freq(x= cytogen$p6)
cytogen.p7 <- cytogen.freq(x= cytogen$p7)
cytogen.p8 <- cytogen.freq(x= cytogen$p8)
cytogen.p9 <- cytogen.freq(x= cytogen$p9)
cytogen.p10 <- cytogen.freq(x= cytogen$p10)
cytogen.p11 <- cytogen.freq(x= cytogen$p11)
cytogen.p12 <- cytogen.freq(x= cytogen$p12)
cytogen.p16 <- cytogen.freq(x= cytogen$p16)
cytogen.p17 <- cytogen.freq(x= cytogen$p17)
cytogen.p18 <- cytogen.freq(x= cytogen$p18)
cytogen.p19 <- cytogen.freq(x= cytogen$p19)
cytogen.p20 <- cytogen.freq(x= cytogen$p20)
cytogen.p21 <- cytogen.freq(x= cytogen$p21)


cytogen.q1 <- cytogen.freq(x = cytogen$q1)
cytogen.q2 <- cytogen.freq(x = cytogen$q2)
cytogen.q3 <- cytogen.freq(x = cytogen$q3)
cytogen.q4 <- cytogen.freq(x= cytogen$q4)
cytogen.q5 <- cytogen.freq(x= cytogen$q5)
cytogen.q6 <- cytogen.freq(x= cytogen$q6)
cytogen.q7 <- cytogen.freq(x= cytogen$q7)
cytogen.q8 <- cytogen.freq(x= cytogen$q8)
cytogen.q9 <- cytogen.freq(x= cytogen$q9)
cytogen.q10 <- cytogen.freq(x= cytogen$q10)
cytogen.q11 <- cytogen.freq(x= cytogen$q11)
cytogen.q12 <- cytogen.freq(x= cytogen$q12)
cytogen.q13 <- cytogen.freq(x= cytogen$q13)
cytogen.q14 <- cytogen.freq(x= cytogen$q14)
cytogen.q15 <- cytogen.freq(x= cytogen$q15)
cytogen.q16 <- cytogen.freq(x= cytogen$q16)
cytogen.q17 <- cytogen.freq(x= cytogen$q17)
cytogen.q18 <- cytogen.freq(x= cytogen$q18)
cytogen.q19 <- cytogen.freq(x= cytogen$q19)
cytogen.q20 <- cytogen.freq(x= cytogen$q20)
cytogen.q21 <- cytogen.freq(x= cytogen$q21)
cytogen.q22 <- cytogen.freq(x= cytogen$q22)


# cytogen.df <- as.data.frame(cytogen.p1, 
#                         cytogen.p2, 
#                         cytogen.p3,
#                         cytogen.p4, 
#                         cytogen.p5,
#                         cytogen.p6,
#                         cytogen.p7,
#                         cytogen.p8, 
#                         cytogen.p9, 
#                         cytogen.p10,
#                         cytogen.p11,
#                         cytogen.p12,
#                         cytogen.p16,
#                         cytogen.p17, 
#                         cytogen.p18,
#                         cytogen.p19, 
#                         cytogen.p20,
#                         cytogen.p21,
#                         cytogen.q1, 
#                         cytogen.q2, 
#                         cytogen.q3, 
#                         cytogen.q4, 
#                         cytogen.q5, 
#                         cytogen.q6, 
#                         cytogen.q7)

# rownames(cytogen.df) <- colnames(cytogen)                        

# nrow(cytogen)




##############################################




