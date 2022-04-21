## Data Processing ##

#' function to build Seurat Object with CITE seq files
#' @param rna gz file name for RNA generated using CITE-seq
#' @param adt gz file name for ADT generated using CITE-seq
#' @return seurat object with RNA and ADT assays
#' @export
buildCITE = function(rna, adt) {
  # unzip .gz files
  GEOquery::gunzip(rna)
  GEOquery::gunzip(adt)

  # get string without .gz
  ind_r = gregexpr(pattern = '.gz', rna)[[1]][1]
  rna_gz = substr(rna, 1, ind_r-1)
  ind_a = gregexpr(pattern = '.gz', adt)[[1]][1]
  adt_gz = substr(adt, 1, ind_a-1)

  r = Seurat::CreateSeuratObject(readRDS(rna_gz)) # build Seurat object for RNA
  a = Seurat::CreateAssayObject(readRDS(adt_gz)) # build assay object for RNA
  r[["ADT"]] = a # add ADT into Seurat object

  return(r) # return seurat object
}
