

### Documented created July 13 2017 to determine best way to label graphics from clinical_data_goi.R 
### for example how to label axes of heatmap from chi.sq function
### For exploratory purposes only
### to be used in conjunction with clinical_data_goi.R file


### Author: Dr Marion Mateos


### Libraries to be used
# install.packages('gplots')
# install.packages('survival')
# install.packages("ComplexHeatmap")

library(NMF)
library(gplots)
library(car)
library(stats)
library(survival)
library(ComplexHeatmap)

### Functions used

source(file = "/home/nmm199/R/MB_RNAseq/Clinical/clin.script/clinical_data_functions_run.R")

### names of functions for info on function see source file
### "chi.sq"
### "cor.result"
### "lin.reg"
### "km.log.test"
### "km.log.test.OS"
### "cox.result.OS"
### "km.log.test.EFS"
### "updatepData"



list.age.cat.infant <- chi.sq(x = matched.test.pData$age.cat.infant, y = matched.goi.vsd.cat)

list.sex <- chi.sq (x = matched.test.pData$sex, y= matched.goi.vsd.cat)

list.mstatus <- chi.sq(x = matched.test.pData$mstatus, y=matched.goi.vsd.cat) 



### experimenting with chi.sq function

chi.sq <- function(x,y, LabRow, LabCol, Title){
  table.temp <- table(x, y) 
  table.temp.perc <- prop.table(table.temp)*100
  summary.table(table.temp)
  chi.test.temp <- chisq.test(table.temp) 
  chi.test.temp.stat <- c(stat=chi.test.temp$statistic, p.value=chi.test.temp$p.value) 
  chi.test.temp.res <- chi.test.temp$residuals
  aheatmap(chi.test.temp$residuals, Rowv=NA, Colv = NA, labRow = LabRow, labCol = LabCol, main = Title)
  list.temp <- list  (table.temp, 
                      table.temp.perc,
                      chi.test.temp,
                      chi.test.temp.res
  )
  
  return(list.temp)
}



list.age.cat.infant <- chi.sq(x = matched.test.pData$age.cat.infant, y = matched.goi.vsd.cat, LabRow = "Gender", LabCol = "Biomarker expression", Title = "Age and Biomarker expression (low vs high)")


### 

# mat = readRDS(paste0(system.file("extdata", package = "ComplexHeatmap"), "/measles.rds"))
# ha1 = HeatmapAnnotation(dist1 = anno_barplot(colSums(mat), bar_width = 1, gp = gpar(col = NA, fill = "#FFE200"), 
     #                                        border = FALSE, axis = TRUE))
# ha2 = rowAnnotation(dist2 = anno_barplot(rowSums(mat), bar_width = 1, gp = gpar(col = NA, fill = "#FFE200"), 
      #                                   border = FALSE, which = "row", axis = TRUE), width = unit(1, "cm"))
# ha_column = HeatmapAnnotation(cn = function(index) {
  # year = as.numeric(colnames(mat))
  # which_decade = which(year %% 10 == 0)
  # grid.text(year[which_decade], which_decade/length(year), 1, just = c("center", "top"))
#})
#Heatmap(mat, name = "cases", col = colorRamp2(c(0, 800, 1000, 127000), c("white", "cornflowerblue", "yellow", "red")),
  #      cluster_columns = FALSE, show_row_dend = FALSE, rect_gp = gpar(col= "white"), show_column_names = FALSE,
   #     row_names_side = "left", row_names_gp = gpar(fontsize = 10),
    #    column_title = 'Measles cases in US states 1930-2001\nVaccine introduced 1961',
     #   top_annotation = ha1, top_annotation_height = unit(1, "cm"),
      #  bottom_annotation = ha_column, bottom_annotation_height = grobHeight(textGrob("1900"))) + ha2
#
decorate_heatmap_body("cases", {
  i = which(colnames(mat) == "1961")
  x = i/ncol(mat)
  grid.lines(c(x, x), c(0, 1), gp = gpar(lwd = 2))
  grid.text("Vaccine introduced", x, unit(1, "npc") + unit(5, "mm"))
})