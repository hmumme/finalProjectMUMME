### data.R ###
## Hope Mumme (hmumme@emory.edu)
## Dataset:
    ## GEO repository for "Single-cell multiomic analysis identifies regulatory programs in mixed-phenotype acute leukemia" publication by Greenleaf group
    ## https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE139369 

## Download Data:
# (1.1) MPAL1 CITE
# wget https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM4138nnn/GSM4138878/suppl/GSM4138878_scRNA_MPAL1_T1.rds.gz
# wget https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM4138nnn/GSM4138878/suppl/GSM4138878_scADT_MPAL1_T1.rds.gz
# wget https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM4138nnn/GSM4138879/suppl/GSM4138879_scRNA_MPAL1_T2.rds.gz 
# wget https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM4138nnn/GSM4138879/suppl/GSM4138879_scADT_MPAL1_T2.rds.gz

# (1.2) MPAL1 ATAC
# wget https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM4138nnn/GSM4138898/suppl/GSM4138898_scATAC_MPAL1_T1.fragments.tsv.gz
# wget https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM4138nnn/GSM4138899/suppl/GSM4138899_scATAC_MPAL1_T2.fragments.tsv.gz

# (2.1) MPAL2 CITE
# wget https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM4138nnn/GSM4138880/suppl/GSM4138880_scRNA_MPAL2_T1.rds.gz
# wget https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM4138nnn/GSM4138880/suppl/GSM4138880_scADT_MPAL2_T1.rds.gz

# (2.2) MPAL2 ATAC
# wget https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM4138nnn/GSM4138900/suppl/GSM4138900_scATAC_MPAL2_T1.fragments.tsv.gz
# wget https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM4138nnn/GSM4138901/suppl/GSM4138901_scATAC_MPAL2_T2.fragments.tsv.gz

## Load Data into R:
library(Seurat)
library(ArchR)