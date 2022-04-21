## Functions to Generate Plots for TF Analysis

#' Function to plot heatmap for top 10 TFs found differentially expressed (RNA) per cluster/grouping
#' @param obj Seurat object with RNA data, active ident should be what markers are comparing (ex: seurat_clusters)
#' @param markers differentially expressed genes by cluster or group
#' @param type type of markers to plot; options are upregulated (up), downregulated (down), or default both
#' @param db TF database with names(db) as TFs
#' @return ggplot object with Heatmap
#' @export
#' @example
#' tfMap(obj, markers, db, type = "up")
#'
tfMap = function(obj, markers, tfs, type="both") {
  markers = markers[markers$gene %in% names(tfs),] # filter to only include genes that are in TF database
  n = length(levels(Seurat::Idents(obj)))
  features = NULL
  for (i in 1:n) {
    ident = levels(Seurat::Idents(obj))[i]
    identMarks = markers[markers$cluster == ident,]
    if (type == "up") {
      identMarks = identMarks[identMarks$avg_log2FC > 0.25,]
      finalMarks = head(identMarks[order(identMarks$avg_log2FC, decreasing = T),"gene"],10)
      title = "Top Upregulated TFs"
    } else if (type == "down") {
      identMarks = identMarks[identMarks$avg_log2FC < -0.25,]
      finalMarks = head(identMarks[order(identMarks$avg_log2FC, decreasing = F),"gene"],10)
      title = "Top Downregulated TFs"
    } else {
      downMarks = identMarks[identMarks$avg_log2FC < -0.25,]
      upMarks = identMarks[identMarks$avg_log2FC > 0.25,]
      finalMarks = c(head(upMarks[order(upMarks$avg_log2FC, decreasing = T),"gene"],5),head(downMarks[order(downMarks$avg_log2FC,decreasing=F),"gene"],5))
      title = "Top Up/Down-regulated TFs"
    }
    features = c(features,finalMarks)
  }
  plot = Seurat::DoHeatmap(obj, features = features) + ggplot2::ggtitle(title)
  return(plot)
}

#' Function to plot a matrix with Jaccard Similarity values based on how similar up-regulated TFs are between clusters/groups
#' @param obj Seurat object with RNA data, active ident should be what markers are comparing (ex: seurat_clusters)
#' @param markers differentially expressed genes by cluster or group
#' @param db TF database with names(db) as TFs
#' @param plot logical to indicate whether to return plot (TRUE) or return similarity matrix (FALSE)
#' @return either ggplot object or Jaccard Similarity matrix
#' @export
#' @example
#' tfSimilarity(obj, markers, db, plot = T) # returns similarity plot
#'
tfSimilarity = function(obj, markers, db, plot = T) {
  markers = markers[markers$gene %in% tfs,] # filter to only include genes that are in TF database
  n = length(levels(Seurat::Idents(obj)))

  # set up similarity matrix
  sim = matrix(data = NA, nrow = n, ncol = n, dimnames = NULL)
  rownames(sim) = factor(levels(Seurat::Idents(obj)))
  colnames(sim) = factor(levels(Seurat::Idents(obj)))

  # find similarity between clusters/groups based on upregulated marker lists
  for (i in 1:nrow(sim)) {
    ident_i = levels(Seurat::Idents(obj))[i]
    iMarks = markers[markers$cluster == ident_i,]
    iMarks = iMarks[iMarks$avg_log2FC > 0.25,"gene"]
    for (j in 1:nrow(sim)) {
      ident_j = levels(Seurat::Idents(obj))[j]
      jMarks = markers[markers$cluster == ident_j,]
      jMarks = jMarks[jMarks$avg_log2FC > 0.25,"gene"]
      int_ij = length(intersect(iMarks,jMarks)) # intersection between i and j
      u_ij = length(unique(c(iMarks,jMarks))) # union of i and j
      sim[i,j] = as.integer(int_ij)/as.integer(u_ij) # jaccard similarity = intersection/union
    }
  }

  if (plot) {
    plot = ggcorrplot::ggcorrplot(sim) + ggplot2::ggtitle("Jaccard Similarity Among Upregulated TFs in Clusters/Groups")
    print("Returning Jaccard Similarity Plot ...")
    return(plot)
  } else {
    print("Returning Jaccard Similarity Matrix ...")
    return(sim)
  }
}

#' Function to make volcano plot for each top upregulated TF in a group/cluster
#' @param targets dataframe output from getTargets()
#' @param ident ident to plot TF-targets for
#' @param filename filename for saving pdf with all plots
#' @param logT log2FC threshold of significance for volcano plot (default = 0.5)
#' @param pT FDR threshold of significance for volcano plot (default = 0.01)
#' @return saves pdf with filename
#' @export
#' @example
#' voltf(targets = df, ident = "4", filename = "test/testVolcano.pdf", logT = 1, pT = 0.01)
#'
voltf = function(targets, db, tfs, ident, filename, logT = 0.5, pT = 0.01) {
  # check if there are atac.targets and tfs present for ident
  if (!any(names(targets) == ident) & !any(names(tfs) == ident)) {
    print("ident not found ...")
    return(NULL)
  }

  identTargets = data.frame(targets[[ident]])
  n = length(tfs[[ident]])
  plots = list()
  for (i in 1:n) {
    tf_i = tfs[[ident]][i]
    targets_i = db[[tf_i]]
    subTargets = identTargets[identTargets$name %in% targets_i,]
    subTargets = getSIG(subTargets, logT, pT)
    plots[[i]] = ggplot(subTargets, aes(x = Log2FC, y = -log(FDR), color = SIG, label = show)) +
      geom_point() + geom_text() +
      scale_color_manual(values=c("blue", "black", "red")) +
      geom_vline(xintercept=c(-logT, logT), col="black") +
      geom_hline(yintercept=-log10(pT), col="black") +
      ggtitle(paste("Targets for TF", tf_i, "- Ident", ident))
  }
  out = gridExtra::marrangeGrob(plots,nrow=1, ncol=2)
  ggsave(filename, out, width = 12, height = 6)
}

#' Function to plot protein levels of significant TF targets in different clusters/groups
#'
#' @param targets
#' @param type
#' @param ident
#' @param logT
#' @param pT
#' @export
#'
plotProteins = function(targets, type, obj, filename, ident, logT = 0.5, pT = 0.01) {
  if (!any(ident == targets$ident)) {
    return("ident not found ...")
  }

  plot_1 = FALSE
  plot_2 = FALSE

  identTargets = targets[targets$ident == ident,]
  n = length(unique(identTargets$TF))
  plots = list()

  if (type == "TF") {
    features = checkFeatures(obj, unique(identTargets$TF), "ADT")

    if (length(features) > 0) {
      plots = VlnPlot(obj, assay = "ADT", features = features)
    }

  } else if (type == "targets") {

    for (i in 1:n) {
      tf_i = unique(identTargets$TF)[i]
      subTargets = identTargets[identTargets$TF == tf_i,]
      subTargets = getSIG(subTargets, logT, pT)
      upTargets = checkFeatures(obj, subTargets$target[subTargets$SIG == "UP"], "ADT")
      downTargets = checkFeatures(obj, subTargets$target[subTargets$SIG == "DOWN"], "ADT")

      if (length(upTargets > 0)) {
        pi_1 = Seurat::VlnPlot(obj, assay = "ADT", features = upTargets, combine = F)
        plot_1 = TRUE
        num_1 = length(pi_1)
        pi_1 = gridExtra::marrangeGrob(pi_1,nrow=1, ncol=num_1, top = paste("Upregulated targets of TF", tf_i, "- Ident", ident))
      }

      if (length(downTargets > 0)) {
        pi_2 = Seurat::VlnPlot(obj, assay = "ADT", features = downTargets, combine = F)
        plot_2 = TRUE
        num_2 = length(pi_2)
        pi_2 = gridExtra::marrangeGrob(pi_2,nrow=1, ncol=num_2, top = paste("Downregulated targets of TF", tf_i, "- Ident", ident))
      }

      if (plot_1) {plots = c(plots,pi_1)}
      if (plot_2) {plots = c(plots,pi_2)}
    }
  }


  if (length(plots) == 0) {
    return("No targets found in protein data ...")
  }

  if (type == "TF") {
    out = plots
  } else {
    out = gridExtra::marrangeGrob(plots,nrow=1, ncol=1)
  }
  ggsave(filename, out, width = 12, height = 6)

}

#' Plots enrichment of TF gene set
#' @param markers atac DE markers
#' @param db.filt TF-target database
#' @param path TF of interest
#' @return enrichment plot generated using fgsea
#' @export

plotEn = function(markers, db.filt, path) {
  # rank all genes based on their fold change in ATAC markers
  genes = markers$Log2FC
  names(genes) = markers$name

  # plot enrichment
  p = fgsea::plotEnrichment(db.filt[[path]], genes) + ggplot2::ggtitle(paste(path, "Target Enrichment"))

  return(p)
}

vol = function(targets, logT, pT) {
  targets = getSIG(targets, logT, pT)
  plots = list()
  i = 1
  for (tf in unique(targets$TF)) {
    plots[[i]] = ggplot(targets, aes(x = log2FC, y = -log(pval), color = SIG, label = show)) +
      geom_point() + geom_text() +
      scale_color_manual(values=c("blue", "black", "red")) +
      geom_vline(xintercept=c(-logT, logT), col="black") +
      geom_hline(yintercept=-log10(pT), col="black") +
      ggtitle(paste("Targets for",tf))
  }
  return(plots)
}

#' function to generate boxplots of estimated protein expression of TF targets
#' @param est estimated protein expressin output from estProteins() function
#' @param tf transcription factor of interest
#' @return boxplot with expression values and signficance from ANOVA or t.test
#' @export
tfBar = function(est, tf) {
  if (length(unique(est$ident)) > 2) {
    compare = "ANOVA"
  } else {
    compare = "t.test"
  }
  p = ggplot2::ggplot(est, ggplot2::aes(x = ident, y = e.protein, fill = ident)) + ggplot2::geom_boxplot() + ggplot2::ggtitle(paste(tf, "target protein estimated expression")) +
    ggpubr::stat_compare_means(method = compare)
  return(p)
}
