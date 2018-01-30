### List of functions for clinical data extract
### note that this is saved on notepad also 

### kaplan meier
extract.km.OS.pval
extract.km.EFS.pval
extract.km.PFS.pval

### cox-related functions
cox.dataframe
extract.cox ### note that functions within this function are located below, # extract.cox.pval (as an example)
extract.cox.SHH.old
extract.cox.OS  ### needed to develop separate functions as there were subsetting ("subscripts out of bounds" errors)

# extract.cox.pval
# extract.cox.Zscore
# extract.cox.HR
# extract.cox.L95CI.HR
# extract.cox.U95CI.HR

# extract.cox.pval.SHH.old
# extract.cox.Zscore.SHH.old
# extract.cox.HR.SHH.old
# extract.cox.L95CI.SHH.old
# extract.cox.U95CI.SHH.old

### there are additional functions that the extract.cox.OS function relies on: 
# extract.cox.OS.pval
# extract.cox.OS.HR
# extract.cox.OS.Zscore
# extract.cox.OS.L95CI
# extract.cox.OS.U95CI

### multivariate cox regression
extract.multivar.cox
# extract.multivar.cox.pval 
# extract.multivar.cox.Zscore
# extract.multivar.cox.L95CI
# extract.multivar.cox.U95CI 

extract.multivar.cox.PFS
# extract.multivar.cox.PFS.pval
# extract.multivar.cox.PFS.HR
# extract.multivar.cox.PFS.L95CI
# extract.multivar.cox.PFS.U95CI
# extract.multivar.cox.PFS.Zscore

extract.multivar.cox.PFS.SHH
# extract.multivar.cox.PFS.SHH.pval
# extract.multivar.cox.PFS.SHH.HR
# extract.multivar.cox.PFS.SHH.L95CI
# extract.multivar.cox.PFS.SHH.U95CI
# extract.multivar.cox.PFS.SHH.Zscore

### logistic regression
log.reg.dataframe

### chi-square
extract.chi.all

########################################################################################
########################################################################################

### List of outputs
km.OS.all.results 
km.OS.G3G4.results 
km.OS.SHH.results 
km.OS.SHH.old.results 

### EFS p values
km.EFS.all.results 
km.EFS.G3G4.results 

### PFS p values
km.PFS.all.results 
km.PFS.G3G4.results
km.PFS.SHH.results
km.PFS.SHH.old.results

### significant kaplan meier survival results (univariate)

significant.km.EFS.all 
significant.km.EFS.G3G4 

significant.km.OS.all 
significant.km.OS.G3G4 
significant.km.OS.SHH 
significant.km.OS.SHH.old 

significant.km.PFS.all              
significant.km.PFS.G3G4         
significant.km.PFS.SHH 
significant.km.PFS.SHH.old 

###
cox.PFS.cat.all.df 
cox.PFS.cont.all.df 
cox.PFS.cat.G3G4.df
cox.PFS.cont.G3G4.df 
cox.PFS.cat.SHH.df
cox.PFS.cont.SHH.df 
cox.PFS.cat.SHH.old.df                                     
cox.PFS.cont.SHH.old.df 

### defined threshold : adjusted p<0.05

sig.cox.PFS.cat.all 
sig.cox.PFS.cont.all 
sig.cox.PFS.cat.G3G4 
sig.cox.PFS.cont.G3G4
sig.cox.PFS.cat.SHH 
sig.cox.PFS.cont.SHH 
sig.cox.PFS.cat.SHH.old 
sig.cox.PFS.cont.SHH.old 

###########################################################################################
### annotated which means that ensembl ID and gene names are included in the object
annot.cox.PFS.cont.all   ### useful for characterising hits later on within multivar cox
annot.sig.cox.PFS.cont.all 
annot.cox.PFS.cat.all
annot.sig.cox.PFS.cat.all
annot.cox.PFS.cat.G3G4
annot.sig.cox.PFS.cat.G3G4

# annot.cox.PFS.cont.SHH ### these objects may not exist
# annot.sig.cox.PFS.cont.SHH  ### these objects may not exist
# annot.cox.PFS.cat.SHH   ### these objects may not exist

annot.cox.PFS.cont.SHH.old 
annot.sig.cox.PFS.cont.SHH.old 
annot.cox.PFS.cat.SHH.old 
annot.sig.cox.PFS.cat.SHH.old 

clean.annot.sig.cox.PFS.cont.all ### with NAs removed   

########################################################################################
########################################################################################
### cox survival modelling objects

cox.OS.cat.all.df 
cox.OS.cont.all.df 
cox.OS.cat.G3G4.df 
cox.OS.cont.G3G4.df 
cox.OS.cat.SHH.df 
cox.OS.cont.SHH.df
cox.OS.cat.SHH.old.df 
cox.OS.cont.SHH.old.df 

##################

sig.cox.OS.cat.all 
annot.sig.cox.OS.cat.all
sig.cox.OS.cont.all
annot.sig.cox.OS.cont.all
sig.cox.OS.cat.G3G4 
annot.sig.cox.OS.cat.G3G4 
annot.sig.cox.OS.cat.G3G4
sig.cox.OS.cont.G3G4 
annot.sig.cox.OS.cat.G3G4
sig.cox.OS.cont.G3G4
sig.cox.OS.cat.SHH 
sig.cox.OS.cont.SHH.df 
sig.cox.OS.cat.SHH.old 
sig.cox.OS.cont.SHH.old.df 


cox.EFS.cat.all.df 
cox.EFS.cat.G3G4.df 
cox.EFS.cat.G3G4.v2.df
sig.cox.EFS.cat.all 
sig.cox.EFS.cat.G3G4
sig.cox.EFS.cat.G3G4.v2 

#################################################################################################
#################################################################################################

### extract logistic regression dataframe objects, with pvalue, adjusted p value, odds ratio(OR), 95CI (upper and lower)

extract.logreg.LCA.df 
extract.logreg.relapse.df 
extract.logreg.mstatus.df 
extract.logreg.age.cat.df 
extract.logreg.meth.df 
extract.logreg.meth7.df 
extract.logreg.MYC.df 
extract.logreg.MYCN.df 
extract.logreg.MYCMYCN.df 
extract.logreg.resection.df 
extract.logreg.sex.df 
extract.logreg.TERT.df 
extract.logreg.TP53.df 

logistic.reg.results ### list of all logistic regression results

#################################################################################################
#################################################################################################
### chi square p value, adjusted p value

chi.age.cat.infant.result 
chi.CSI.result 
chi.LCA.result 
chi.meth4.result 
chi.meth7.result 
chi.mstatus.result 
chi.MYC.result 
chi.MYCMYCN.result 
chi.MYCN.result 
chi.q13loss.result 
chi.relapse.result 
chi.resection.result 
chi.RTX.result 
chi.sex.result 
chi.TERT.result 
chi.TP53.result 

### extract adj p <0.05 for relapse, mstatus, MYC, MYCN, MYCMYCN

significant.chi.relapse ### n=4388 4/12/17 for mb.vsd
significant.chi.mstatus ### n=3875, 4/12/17 for mb.vsd
significant.chi.MYC   ### n=4640, 4/12/17 for mb.vsd
significant.chi.MYCN 
significant.chi.MYCMYCN ### n=214 4/12/17 for mb.vsd


########################################################################
########################################################################

###  multivariate cox, looking for transcripts that are significant beyond either the current PNET5, the Lancet oncology paper (Schwalbe et al 2017) or a combined model taking both models together

multivar.cox.OS.combined.cat.df  ### updated so that p value is for biomarker not overall model p val 21/11/17
multivar.cox.OS.combined.cont.df 
multivar.cox.OS.lancetG3G4.cat.df 
multivar.cox.OS.lancetG3G4.cont.df 
multivar.cox.OS.PNET5.cat.df 
multivar.cox.OS.PNET5.cont.df 
multivar.cox.OS.SHHold.cat.df 
multivar.cox.OS.SHHold.cont.df 
multivar.cox.PFS.combined.cat.df 
multivar.cox.PFS.combined.cont.df 
multivar.cox.PFS.lancetG3G4.cat.df 
multivar.cox.PFS.lancetG3G4.cont.df 
multivar.cox.PFS.PNET5.cat.df 
multivar.cox.PFS.PNET5.cont.df 
multivar.cox.PFS.SHHold.cat.df 
multivar.cox.PFS.SHHold.cont.df 

### generating significant dataframes for the multivariate cox modelling, ie transcripts that perform about and beyond current clinical risk models

significant.multivar.cox.OS.combined.cat 
significant.multivar.cox.OS.combined.cont  ###n=13 21/12/17 for all transcripts mb.vsd
significant.multivar.cox.OS.lancetG3G4.cat 
significant.multivar.cox.OS.lancetG3G4.cont
significant.multivar.cox.OS.PNET5.cat       ### n=43 for mb.vsd
significant.multivar.cox.OS.PNET5.cont      ### n=30 for mb.vsd
significant.multivar.cox.OS.SHHold.cat 
significant.multivar.cox.OS.SHHold.cont 