### Multivariate analysis working document

### Code extracted from larger script file clinical_data_4.R
### created July 4th 2017
### Author: Dr Marion Mateos


### multivariate analysis
### to use R+, M+, LCA , sex, MYCMYCN.cat, TP53 across entire group, restricted to age 3-16 yo
### cross reference with dataset of curative intent (CSI and max resection)


### for relapse

log.reg.multi <- glm(matched.test.incl.pData$relapse ~ matched.test.incl.pData$resection 
                     + matched.test.incl.pData$mstatus 
                     + matched.test.incl.pData$sex 
                     + matched.test.incl.pData$resection
                     + matched.test.incl.pData$histopath
                     + matched.test.incl.pData$MYC.cat
                     + matched.test.incl.pData$MYCN.cat
                     + matched.test.incl.pData$TP53.cat, family=binomial (link = 'logit'),
                     data = matched.test.incl.pData)

summary(log.reg.multi)

### compare to model with inclusion of biomarker

log.reg.multi.goi <- glm(matched.test.incl.pData$relapse ~ matched.test.incl.pData$resection 
                         + matched.test.incl.pData$mstatus 
                         + matched.test.incl.pData$sex 
                         + matched.test.incl.pData$resection
                         + matched.test.incl.pData$histopath
                         + matched.test.incl.pData$MYC.cat
                         + matched.test.incl.pData$MYCN.cat
                         + matched.test.incl.pData$TP53.cat
                         + matched.goi.vsd.cat.incl, 
                         family=binomial (link = 'logit'),
                         data = matched.test.incl.pData)

summary(log.reg.multi.goi)

### consider breaking down into subgroups, as the above may overwhelm the dataset with missing values

message ("processing subgroup specific biomarkers")


##########################################

### backward elimination to remove least significant variable at each stage, model with G3G4
### assign levels to relevant factors, below function does not work

MYC.G3G4.cat <- relevel (matched.G3G4.incl.pData$MYC.cat, ref="MYC non ampl")
mstatus.G3G4 <- relevel (matched.G3G4.incl.pData$mstatus, ref = "M0/M1")
resection.G3G4 <- relevel (matched.G3G4.incl.pData$resection, ref = "Gross total resection")
histopath.G3G4 <- relevel (matched.G3G4.incl.pData$histopath, ref = "FAV")
sex.G3G4 <- relevel (matched.G3G4.incl.pData$sex, ref = "female")
MYCN.G3G4.cat <- relevel (matched.G3G4.incl.pData$MYCN.cat, ref = "MYCN non ampl")

# relevel.G3G4 <- function (x,y){
# relevel.temp <- relevel (x, ref = "y")
# return (relevel.temp)
# }

# MYC.G3G4.cat <- relevel.G3G4(x = matched.G3G4.incl.pData$MYC.cat, y = "MYC non ampl")


log.reg.multi.G3G4.1 <- glm(matched.G3G4.incl.pData$relapse ~ resection.G3G4 
                            + mstatus.G3G4
                            + sex.G3G4
                            + resection.G3G4
                            + histopath.G3G4
                            + MYC.G3G4.cat
                            + MYCN.G3G4.cat,
                            family=binomial (link = 'logit'),
                            data = matched.G3G4.incl.pData)

summary(log.reg.multi.G3G4.1)
confint (log.reg.multi.G3G4.1)
exp(coef(log.reg.multi.G3G4.1))
exp(cbind(OR = coef(log.reg.multi.G3G4.1), confint(log.reg.multi.G3G4.1)))


### question: proceed in this manner, and do manually at the time? 


### then compare to the addition of the biomarker, perform backward elimination


log.reg.multi.G3G4.1A <- glm(matched.G3G4.incl.pData$relapse ~ resection.G3G4 
                             + mstatus.G3G4
                             + sex.G3G4
                             + resection.G3G4
                             + histopath.G3G4
                             + MYC.G3G4.cat
                             family=binomial (link = 'logit'),
                             data = matched.G3G4.incl.pData)

summary(log.reg.multi.G3G4.1A)
confint (log.reg.multi.G3G4.1A)
exp(coef(log.reg.multi.G3G4.1A))
exp(cbind(OR = coef(log.reg.multi.G3G4.1A), confint(log.reg.multi.G3G4.1A)))



log.reg.multi.G3G4.1B <- glm(matched.G3G4.incl.pData$relapse ~ matched.G3G4.incl.pData$resection 
                             + matched.G3G4.incl.pData$mstatus 
                             + matched.G3G4.incl.pData$sex 
                             + matched.G3G4.incl.pData$resection
                             + matched.G3G4.incl.pData$histopath
                             + matched.G3G4.incl.pData$MYC.cat,
                             family=binomial (link = 'logit'),
                             data = matched.G3G4.incl.pData)


summary(log.reg.multi.G3G4.1B)


### compare with addition of biomarker

log.reg.multi.goi.G3G4 <- glm(matched.G3G4.incl.pData$relapse ~ matched.G3G4.incl.pData$resection 
                              + matched.G3G4.incl.pData$mstatus 
                              + matched.G3G4.incl.pData$sex 
                              + matched.G3G4.incl.pData$resection
                              + matched.G3G4.incl.pData$histopath
                              + matched.G3G4.incl.pData$MYC.cat
                              + matched.goi.vsd.cat.G3G4.incl, 
                              family=binomial (link = 'logit'),
                              data = matched.G3G4.incl.pData)


### other subgroups e.g SHH TP53




# summary(summary.cox)

### function not currently running

### multivariate cox regression analysis


time <- matched.test.incl.pData$PFS
# event <- ifelse(matched.test.incl.pData$relapse=="relapse", 1, 0)
event <- relapse.bin.incl
y <-matched.test.incl.pData$histopath
a <-matched.test.incl.pData$mstatus
b <-matched.test.incl.pData$resection
c <-matched.test.incl.pData$MYC.cat 
d <-matched.test.incl.pData$MYCN.cat

View(test.pData)

cox.result <- function (time, event, y, a, b, c, d)  
{
  cox.temp <- coxph (Surv(time, event)~y + a + b + c + d, data= matched.test.incl.pData)
  summary.cox <- summary(cox.temp)$logtest
  return (summary(summary.cox))
}


cox.result.1 <- cox.result (time = matched.test.incl.pData$PFS, event = relapse.bin.incl, 
                            y = matched.test.incl.pData$histopath,
                            a = matched.test.incl.pData$mstatus,
                            b = matched.test.incl.pData$resection,
                            c = matched.test.incl.pData$MYC.cat,
                            d = matched.test.incl.pData$MYCN.cat
                            #  )
                            
                            #summary(cox.result.1)
                            
                            ### only input values below pre-defined p threshold for multivariate model
                            ### let's use p<0.2 from the univariate logistic regression
                            
                            cox.multi <- coxph (Surv(matched.test.incl.pData$PFS, relapse.bin.incl) ~ )
                            summary(cox.multi)$logtest
                            summary(cox.multi)
                            
                            
                            cox.multi.marker <- coxph (Surv(matched.test.incl.pData$PFS, relapse.bin.incl) ~ matched.goi.vsd.cat.incl)
                            summary(cox.multi)$logtest
                            summary(cox.multi)
                            
                            
                            ### also compare the model to that with the addition of a single biomarker
                            
                            
                            
                            
                            
                            
                            ###################################
                            ### AUC and ROC analysis
                            
                            
                            
                            
                            
                            
                            