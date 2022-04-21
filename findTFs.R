## Functions to Generate top TFs or TF Targets

#' Function to retrieve top upregulated TFs for RNA clusters/groups
#' @param obj Seurat object with RNA data, active ident should be what markers are comparing (ex: seurat_clusters)
#' @param markers differentially expressed genes by cluster or group
#' @param db TF-target database with names(db) as TFs
#' @param num number of TFs to export per cluster/group, default is 20
#' @export
getTFs = function(obj, markers, tfs, num = 20) {
  n = length(levels(Seurat::Idents(obj)))
  out = vector("list", n) # create empty list of length n to hold top TFs per cluster/group
  if ("cluster" %in% colnames(markers)) {
    markers = markers[markers$gene %in% names(tfs),] # filter to only include genes that are in TF database
    for (i in 1:n) {
      ident = levels(Seurat::Idents(obj))[i]
      identMarks = markers[markers$cluster == ident,]
      identMarks = identMarks[identMarks$avg_log2FC > 0.25,]
      if (num == "all") {
        out[[i]] = identMarks[order(identMarks$avg_log2FC, decreasing = T),"gene"]
      } else {
        out[[i]] = head(identMarks[order(identMarks$avg_log2FC, decreasing = T),"gene"], num)
      }

    }
  } else {
    markers = markers[rownames(markers) %in% names(tfs),] # filter to only include genes that are in TF database
    ident1_markers = markers[markers$avg_log2FC > 0.25,]
    ident2_markers = markers[markers$avg_log2FC < 0.25,]

    if (num == "all") {
      out[[1]] = rownames(ident1_markers[order(ident1_markers$avg_log2FC, decreasing = T),])
      out[[2]] = rownames(ident2_markers[order(ident2_markers$avg_log2FC, decreasing = T),])
    } else {
      out[[1]] = rownames(head(ident1_markers[order(ident1_markers$avg_log2FC, decreasing = T),], num))
      out[[2]] = rownames(head(ident2_markers[order(ident2_markers$avg_log2FC, decreasing = T),], num))
    }
  }

  names(out) = levels(Seurat::Idents(obj))
  return(out)
}

#' Function to retrieve Log2FC values for TF targets in ATAC based on input topTFs for cluster/groups in RNA
#'
#' @param tfs output from getTFs, list of top TFs per cluster/group
#' @param markers ATAC top markers per cluster/group, output from ArchR::getMarkers() function
#' @param db TF-target database with names(db) as TFs
#' @export
getTargets = function(tfs, markers, db) {
  out = data.frame(ident = NA, TF = NA, target = NA, atac_log2FC = NA, atac_FDR = NA)
  out_row = 2
  n = length(tfs) # number of idents

  for (i in 1:n) {
    ident_i = names(tfs)[i]
    m = length(tfs[[i]]) # number of TFs for ith ident
    if (any(names(markers) == ident_i)) {
      atac_i = markers[[ident_i]]
      if (nrow(atac_i) == 0) {print(paste("ident", ident_i, "not found in ATAC top markers..."))}
      for (j in 1:m) {
        tf_j = tfs[[i]][j] # extract jth TF for ith ident TFs
        tf_targets = db[names(db) == tf_j][[1]] # find targets for tf_j
        tf_targets = tf_targets[tf_targets %in% atac_i$name] # filter to only include targets present in atac de

        for (k in tf_targets) {
          out[out_row,"ident"] = ident_i
          out[out_row,"TF"] = tf_j
          out[out_row,"target"] = k
          out[out_row,"atac_log2FC"] = atac_i$Log2FC[atac_i$name == k]
          out[out_row,"atac_FDR"] = atac_i$FDR[atac_i$name == k]
          out_row = out_row + 1
        }
      }
    } else {
      print(paste("ident", ident_i, "not found in ATAC top markers..."))
    }

  }
  return(out[-1,])
}

#' @export
filterTargets = function(tfs, markers, db) {
    for (i in names(markers)) {
    tfs_i = unique(unlist(tfs[[i]]))
    ts = unlist(unique(db[names(db) %in% tfs_i]))
    markers_i = markers[[i]]
    markers[[i]] = markers_i[markers_i$name %in% ts,]
  }
  return(markers)
}

#' @export
tfEnrichment = function(markers, db.filt) {
  # rank all genes based on their fold change in ATAC markers
  genes = markers$Log2FC
  names(genes) = markers$name

  # run fgsea
  en = fgsea::fgsea(pathways = db.filt, stats = genes)
  return(en)
}

#' @export
getProteins = function(markers, db.filt) {
  out = data.frame(target = NA, log2FC = NA, pval = NA, TF = NA)
  i = 1
  for (tf in names(db.filt)) {
    targets_tf = db.filt[[tf]]
    targets_tf = targets_tf[targets_tf %in% rownames(markers)]
    for (t in targets_tf) {
      out[i,'target'] = t
      out[i,'log2FC'] = markers$avg_log2FC[rownames(markers) == t]
      out[i,'pval'] = markers$p_val_adj[rownames(markers) == t]
      out[i,'TF'] = tf
      i = i + 1
    }
  }
  return(out)
}

#' Estimate protein expression by calculating scaling factor per cluster
#' @param obj
#' @param ident
#' @param sampleSize
#' @export
#'
getScale = function(obj, ident, sampleSize = 100) {
  # extract adt and rna
  Seurat::DefaultAssay(obj) = "ADT"
  adt = Seurat::GetAssayData(obj, slot = "counts")
  Seurat::DefaultAssay(obj) = "RNA"
  rna = Seurat::GetAssayData(obj, slot = "counts")

  proteins = rownames(adt)

  idents = unlist(unique(obj[[ident]]))
  ss = sampleSize

  Seurat::Idents(obj) = ident

  x = rep(0, length(idents))
  i = 1
  cluster_ids = list()

  for (c in idents) {
    print(paste("Calculating scale factor for ident", c, "..."))
    c_ids = Seurat::WhichCells(obj, idents = c)
    cluster_ids[[c]] = c_ids
    c_sample = sample(c_ids, size = ss)

    c_rna = base::rowMeans(rna[rownames(rna) %in% proteins, colnames(rna) %in% c_sample], na.rm = T)
    c_adt = base::rowMeans(adt[rownames(adt) %in% proteins, colnames(adt) %in% c_sample], na.rm = T)
    c_adt = c_adt[!is.na(c_adt)]

    c_sum = 0
    c_count = 0

    for (p in proteins) {
      if ((p %in% names(c_rna)) & (p %in% names(c_adt))) {
        c_sum = c_sum + unlist(c_rna[p] / c_adt[p])
        c_count = c_count + 1
      }
    }

    x[i] = c_sum / c_count

    i = i + 1
  }
  names(x) = idents
  return(x)
}

#' Estimate target protein expression for TF
#' @param x vector with protein-gene scaling factor per cluster
#' @param obj
#' @param db.filt
#' @param tf transcription factor to estimate
#' @export
estProtein = function(x, obj, db.filt, tf, ident) {
  # extract adt
  Seurat::DefaultAssay(obj) = "ADT"
  adt = Seurat::GetAssayData(obj, slot = "counts")

  # get 100 random target genes for tf
  if (length(db.filt[[tf]]) > 100) {
    genes = sample(db.filt[[tf]],100)
  } else {
    genes = db.filt[[tf]]
  }

  rna.avg = Seurat::AverageExpression(obj, assays = "RNA", group.by = ident, features = genes)$RNA # average expression per cluster

  est = data.frame("target" = rep(rownames(rna.avg), length(x)))
  nGenes = nrow(rna.avg)
  est$rna.exp = NA
  est$e.protein = NA
  est$ident = NA

  # estimate protein values for 100 target genes of tf
  for (g in rownames(rna.avg)) {
    est$ident[est$target == g] = names(x)
    for (c in names(x)) {
      g_exp = rna.avg[rownames(rna.avg) == g, colnames(rna.avg) == c]
      est$rna.exp[est$target == g & est$ident == c] = g_exp
      est$e.protein[est$target == g & est$ident == c] = g_exp / x[names(x) == c]
    }
  }
  return(est)
}
