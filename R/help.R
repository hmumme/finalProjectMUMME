## Helper Functions for findTFs.R and plotting.R

#' Function to assist with voltf() and plotProteins()
#' @param df dataframe output from getTargets()
#' @param log2FC threshold for significance
#' @param pT FDR threshold for significance
#' @return a dataframe with columns added SIG (UP, DOWN, NS) and show (label names)
#'
getSIG = function(df, logT, pT) {
  df['SIG'] = "NS" # assign all as not significant to start
  df$SIG[df$log2FC > logT & df$pval < pT] <- "UP"
  df$SIG[df$log2FC < -logT & df$pval < pT] <- "DOWN"
  df['show'] = NA
  df$show[df$SIG != "NS"] = df$target[df$SIG != "NS"]
  return(df)
}

#' Function to return list of features found in assay for Seurat object
#' @param obj Seurat object
#' @param features features to check for
#' @param assay assay in Seurat object to check for features
#' @return filtered list of features with only features found in object's assay
#'
checkFeatures = function(obj, features, assay) {
  out = c()
  for (f in features) {
    if (any(rownames(obj[[assay]]) == f)) out = c(out,f)
  }
  return(out)
}

# go to gene information on ma'ayan_lab website
geneInfo = function(gene) {
  browseURL(paste0('https://maayanlab.cloud/archs4/gene/',gene))
}
