# finalProjectMUMME
Final Project for BMI585-Spring2022

Hope Mumme (hmumme@emory.edu)

Provides functions to predict transcription factor and target activities using Seurat (CITE) and ArchR (ATAC)

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
See demo folder with demo.Rmd, demo.html, and demo.pdf for an example workthrough of analysis.

The following data is needed to run the demo:
- Gene Set Library ENCODE Transcription Factor Targets Dataset from the [Harmonizome portal](https://maayanlab.cloud/Harmonizome/dataset/ENCODE+Transcription+Factor+Targets)
- fragments files and rds files for BMMC_D1T1, BMMC_D1T2, MPAL1_T1, MPAL1_T2, MPAL2_T1, and MPAL2_T2 samples from GSE139369 on the [GEO accession viewer](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE139369)
