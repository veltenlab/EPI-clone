plot_tfbs_methylation <- function(seurat_rna,
                                  counts,
                                  tfbs,
                                  panel,
                                  point_size=.5){
  require(viridis)
  panel <- read.table(panel, sep='\t')
  sel_amplis <- intersect(row.names(panel[!is.na(panel[, tfbs]), ]), colnames(counts))
  print(paste('Using', length(sel_amplis), 'amplicons'))
  meth_tfbs <- rowMeans(counts[, sel_amplis, drop=FALSE])/rowMeans(counts)
  to_plot <- as.data.frame(seurat_rna[['umap']]@cell.embeddings)
  to_plot <- data.frame(to_plot, TFBS=meth_tfbs)
  colnames(to_plot)[ncol(to_plot)] <- tfbs
  plot <- ggplot(to_plot, aes_string(x='UMAP_1', y='UMAP_2', color=tfbs))+geom_point(size=.35, stroke=.35)+
    plot_theme+scale_color_viridis(option='mako', direction=-1, name="Relative\nMethylation")
  
  return(plot)
}

plot_type_methylation <- function(seurat_rna,
                                  counts,
                                  type,
                                  panel,
                                  point_size=.5){
  require(viridis)
  panel <- read.table(panel, sep='\t')
  sel_amplis <- intersect(row.names(panel[which(panel$Type==type), ]), colnames(counts))
  print(paste('Using', length(sel_amplis), 'amplicons'))
  meth_tfbs <- rowMeans(counts[, sel_amplis, drop=FALSE])/rowMeans(counts)
  to_plot <- as.data.frame(seurat_rna[['umap']]@cell.embeddings)
  to_plot <- data.frame(to_plot, TFBS=meth_tfbs)
  colnames(to_plot)[ncol(to_plot)] <- type
  plot <- ggplot(to_plot, aes_string(x='UMAP_1', y='UMAP_2', color=type))+geom_point(size=.35, stroke=.35)+
    plot_theme+scale_color_viridis(option='mako', direction=-1, name="Relative\nMethylation")
  
  return(plot)
}
run_marker_enrichment <- function(markers, panel, fc.cut=1.0){
  pos.markers <- markers[markers$avg_log2FC>fc.cut&markers$p_val_adj<0.05, ]
  if(nrow(pos.markers)<1){
    pos_homer <- NA
    pos_tfbs <- NA
    meth_enhancers <- NA
  }else{
    pos.markers <- panel[row.names(pos.markers), ]
    pos_homer <- NA
    pos_tfbs <- plyr::count(unlist(pos.markers[, tfbs_columns]))
    pos_tfbs$enrichment <- c(enrich_tfbs(pos_tfbs, panel, nrow(pos.markers)), NA)
    pos_tfbs[which(pos_tfbs$freq<=3), 'enrichment'] <- 1
    meth_enhancers <- pos.markers$enhancer_annotation_gene_name[!is.na(pos.markers$enhancer_annotation_gene_name)]
  }
  neg.markers <- markers[markers$avg_log2FC<(-fc.cut)&markers$p_val_adj<0.05, ]
  if(nrow(neg.markers)<1){
    neg_homer <- NA
    neg_tfbs <- NA
    demeth_enhancers <- NA
  }else{
    neg.markers <- panel[row.names(neg.markers), ]
    neg_homer <- NA
    neg_tfbs <- plyr::count(unlist(neg.markers[, tfbs_columns]))
    neg_tfbs[which(neg_tfbs$freq<=3), 'enrichment'] <- 1
    neg_tfbs$enrichment <- c(enrich_tfbs(neg_tfbs, panel, nrow(neg.markers)), NA)
    demeth_enhancers = neg.markers$enhancer_annotation_gene_name[!is.na(neg.markers$enhancer_annotation_gene_name)]
  }
  return(list(Positive=list(HOMER=pos_homer, TFBS=na.omit(pos_tfbs), Meth_enhancers=meth_enhancers),
              Negative=list(HOMER=neg_homer, TFBS=na.omit(neg_tfbs), Demeth_enhancers=demeth_enhancers)))
}
#tfbs_columns <- c(26:73, 79:96)
tfbs_columns <- c(19:66, 69:86)

enrich_tfbs <- function(tfbs_count, panel, num){
  if(nrow(tfbs_count)<1){
    return(NULL)
  }
  p_vals <- c()
  for(i in 1:nrow(tfbs_count)){
    tfbs <- tfbs_count[i, 'x']
    if(is.na(tfbs)) next
    tfbs <- gsub('-', '.', tfbs)
    tp <-  tfbs_count[i, 'freq']
    bkg <- sum(!is.na(panel[, tfbs]))
    if(bkg<=1){
      p_vals <- c(p_vals, NA)
      next
    }
    uni <- nrow(panel)-bkg
    p_val <- 1-phyper(tp, bkg, uni, num)
    p_vals <- c(p_vals, p_val)
  }
  return(p_vals)
}
differential_test <- function(seurat_rna,
                              ident.1,
                              ident.2=NULL,
                              logfc.threshold=0.1){
  counts <- t(as.matrix(GetAssayData(seurat_rna, assay = 'DNAm', slot='counts')))
  selected_cells <- Idents(seurat_rna)%in%ident.1
  if(is.null(ident.2)){
    other_cells <- !selected_cells
  }else{
    other_cells <- Idents(seurat_rna)%in%ident.2
  }
  res <- sapply(colnames(counts), function(amplicon){
    p_val <- wilcox.test(counts[selected_cells, amplicon], counts[other_cells, amplicon])$p.value
    log2_fc <- log2(mean(counts[selected_cells, amplicon])/mean(counts[other_cells, amplicon]))
    return(c(p_val, log2_fc))
  })
  res <- as.data.frame(t(res))
  colnames(res) <- c('p_val', 'avg_log2FC')
  res$p_val_adj <- p.adjust(res$p_val, method='fdr')
  return(res[abs(res$avg_log2FC)>logfc.threshold, ])
}
