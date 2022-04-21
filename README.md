# finalProjectMUMME
Final Project for BMI585-Spring2022

Hope Mumme (hmumme@emory.edu)

## Installation:
```
devtools::install_github("hmumme/finalProjectMUMME")
```
Required Libraries:
- Seurat
- ArhcR
- dplyr
- fgsea

## Functions:
```
buildCITE(rna, adt)
getTFs(object, rna.markers, tf.db, num = all)
tfEnrichment(atac.markers, db.filt)
plotEn(atac.markers, db.filt, TF.name)
getProteins(adt.markers, db.filt)
getScale(object, ident)
estProtein(x, object, db.filt, TF.name, ident)
tfBar(est, TF.name)
```
## Demo:
See demo.Rmd, demo.RData, and demo.pdf for example workthrough of analysis.
