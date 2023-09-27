require(Seurat)
require(ggplot2)
require(plyr)
require(reshape2)
require(Matrix)
require(infotheo)
require(parallel)
require(pheatmap)
require(RColorBrewer)
require(gridExtra)
require(fossil)
require(ROCR)

epiclone <- function(seurat_obj, trueClone = "cleanClone", batch = "Sample", celltype = "seurat_clusters",
                  plotFolder = "plots", tuneParams = F, tuneParams.cores = 6,
                  methylation.assay.name = "DNAm", protein.assay.name = "AB",
                  lower.thr.methrate = 0.25, upper.thr.methrate = 0.9, thr.protein.ass = NULL,selected.CpGs = NULL,
                  ncells.bigClone = 5, npcs.bigCloneSelection = 100, k.bigCloneSelection = 5, thr.bigCloneSelection =1, smoothen.bigCloneSelection = NULL,
                  npcs.Clustering = 100, k.Clustering = 25, res.Clustering=5, returnIntermediateSeurat = F) {
  suppressMessages({
    
    if (!is.null(trueClone)) seurat_obj@meta.data[,trueClone] <- as.character(seurat_obj@meta.data[,trueClone])
    
    if (is.null(trueClone) & tuneParams) {
      cat("Cannot tune parameters when no clone information is available\n")
      tuneParams <- F
    }
    
    if(!is.null(plotFolder)) {
      if (!dir.exists(plotFolder)) dir.create(plotFolder, recursive = T)
    }
    
    plots <- list()
    results <- list()
    DNAm <- seurat_obj@assays[[methylation.assay.name]]@data
    if (!is.null(protein.assay.name)) protein <- seurat_obj@assays[[protein.assay.name]]@data
    
    
    if (is.null(selected.CpGs)) {
    cat("Selecting CpGs...\n")
    if (!is.null(protein.assay.name)) {
      pvals <- apply(DNAm,1, function(met) {
        
        apply(protein, 1, function(prot) {
          use <- !is.na(prot)
          a <- prot[use][met[use]==1]
          b <- prot[use][met[use]==0]
          if (length(a) < 3 | length(b) < 3) return(1) else return(ks.test(a,b)$p.value)
        })
      })
      min_pval <- apply(pvals,2,min)
      if (is.null(thr.protein.ass)) thr.protein.ass <- 1/(nrow(protein) * nrow(DNAm))
    } else {
      cat("No protein data available. Assuming Idents() slot of Seurat object carries cell state annotation...\n")
      globalid <- Idents(seurat_obj)
      min_pval <- p.adjust(apply(DNAm,1,function(met) {
        chisq.test(table(met,globalid))$p.value
      }),method = "bonferroni")
      if (is.null(thr.protein.ass)) thr.protein.ass <- 0.001
    }
    
    
    avg_meth_rate <- apply(DNAm, 1, mean)

    #if true clonal info is present, compute association with clone
    if(!is.null(trueClone)) {
      #for this use only clones > 5
      true_clone <- seurat_obj@meta.data[,trueClone]
      for_prediction <- seurat_obj
      
      for_prediction$use <- !is.na(true_clone)
      for_prediction <- subset(for_prediction, use)
      a <-table(for_prediction@meta.data[,trueClone])
      use_for_prediction <- names(a)[a > ncells.bigClone]
      for_prediction$use <- for_prediction@meta.data[,trueClone] %in% use_for_prediction
      for_prediction <- subset(for_prediction,use)
      
      
      cloneid <- factor(for_prediction@meta.data[,trueClone], levels = unique(for_prediction@meta.data[,trueClone]))
      pvals_cloneass <- p.adjust(apply(for_prediction@assays[[methylation.assay.name]]@data,1,function(met) {
        chisq.test(table(met,cloneid))$p.value
      }),method = "bonferroni")
      plots[["CpGSelection"]] <- qplot(x = avg_meth_rate, y = ifelse(min_pval<1e-30, 1e-30, min_pval), color = -log10(pvals_cloneass+1e-50), log ="y")+ 
        geom_hline(yintercept = thr.protein.ass) + geom_vline(xintercept = c(lower.thr.methrate,upper.thr.methrate)) +
        theme_bw() + theme(panel.grid = element_blank()) + xlab("Avg methylation rate") + ylab(ifelse(is.null(thr.protein.ass), "p value cell state association", "Min. p value protein association")) +
        scale_color_gradientn(colours = c("black","blue","red"), name = "-log10 p-val\nClone association")
      results[["CpGSelection"]] <- data.frame(CpG = names(pvals_cloneass), avg_meth_rate, min_pval, pvals_cloneass)
    } else {
      plots[["CpGSelection"]] <- qplot(x = avg_meth_rate, y = ifelse(min_pval<1e-30, 1e-30, min_pval), log ="y")+ 
        geom_hline(yintercept = thr.protein.ass) + geom_vline(xintercept = c(lower.thr.methrate,upper.thr.methrate)) +
        theme_bw() + theme(panel.grid = element_blank()) + xlab("Avg methylation rate") + ylab(ifelse(is.null(thr.protein.ass), "p value cell state association", "Min. p value protein association"))
    }
    
    
    selected_not_protein <- names(min_pval)[min_pval > thr.protein.ass & avg_meth_rate < upper.thr.methrate & avg_meth_rate > lower.thr.methrate]
    cat(length(selected_not_protein), "CpGs selected\n")
    results[["selectedCpGs"]] <- selected_not_protein
    } else selected_not_protein <- selected.CpGs
    
    
    npcs.compute <- max(c(100, npcs.bigCloneSelection, npcs.Clustering))
    if (npcs.compute > length(selected_not_protein)) {
      cat("Can only compute ", length(selected_not_protein), "Principal components, adjusting parameters\n")
      npcs.compute <- length(selected_not_protein)-3
      if (npcs.bigCloneSelection > npcs.compute) npcs.bigCloneSelection <- npcs.compute
      if (npcs.Clustering > npcs.compute) npcs.Clustering <- npcs.compute
    }
    
    
    #### 3. PC and clustering based on these CpGs #####
    cat("Identifying big clones...\n")
    
    seurat_obj <- ScaleData(seurat_obj, assay = methylation.assay.name, features = selected_not_protein)
    seurat_obj <- RunPCA(seurat_obj, assay = methylation.assay.name, features = selected_not_protein, reduction.name = "clonePCA", reduction.key = "CLONEPC_", npcs = npcs.compute)
    seurat_obj <- RunUMAP(seurat_obj,reduction = "clonePCA", dims = 1:npcs.bigCloneSelection, reduction.name = "cloneUMAP", reduction.key = "CLONEUMAP_" )
    
    
    #### 4. compute point density ####
    seurat_obj <- FindNeighbors(seurat_obj, reduction = "clonePCA", dims = 1:npcs.bigCloneSelection, k.param = k.bigCloneSelection, return.neighbor=T, graph.name = "clone.neighbors" )
    seurat_obj$avgNNdist <- apply(seurat_obj@neighbors$clone.neighbors@nn.dist,1,function(x) mean(x[x>0]))
    
    to.summarise <- seurat_obj@meta.data[, c(batch ,celltype,"avgNNdist",sprintf("nFeature_%s", methylation.assay.name))]
    colnames(to.summarise)[1] <- "Sample"
    colnames(to.summarise)[2] <- "celltype"
    colnames(to.summarise)[4] <- "nFeature"
    if (length(unique(to.summarise$Sample)) == 1) m <- lm(avgNNdist ~ nFeature + celltype, data = to.summarise) else m <- lm(avgNNdist ~ nFeature + Sample + celltype, data = to.summarise)
    seurat_obj$avgNNres <- residuals(m)
    
    if (!is.null(smoothen.bigCloneSelection)) {
      seurat_obj <- FindNeighbors(seurat_obj, reduction = "clonePCA", dims = 1:npcs.bigCloneSelection, k.param = smoothen.bigCloneSelection, return.neighbor=T, graph.name = "clone.neighbors" )
      seurat_obj$avgNNres <- apply(seurat_obj@neighbors$clone.neighbors@nn.idx,1, function(x) mean(seurat_obj$avgNNres[x]))
    }
    
    to.summarise$avgNNres <- seurat_obj$avgNNres
    
    if(!is.null(trueClone)) {
      to.summarise$cleanClone <- seurat_obj@meta.data[,trueClone]
      csize <- table(to.summarise$cleanClone)
      seurat_obj$csize <- csize[as.character(seurat_obj@meta.data[,trueClone])] -> to.summarise$csize
      to.summarise$csize[is.na(to.summarise$csize)] <- 0
      
      
      plots[["BigCloneSelection.ParameterSelection"]] <- ggplot(aes(x = csize+1, y = avgNNres), data = to.summarise) + scale_x_log10() + geom_point(position = position_jitter(width=0.1),size=.5) + ylab("Distance to nearest neighbors") + xlab("LARRY clone size") + geom_hline(yintercept = thr.bigCloneSelection, color = "red")
      
    }
    plots[["BigCloneSelection.Density"]] <- FeaturePlot(seurat_obj, "avgNNres",reduction = "cloneUMAP") + scale_color_gradientn(colours = rev(c("red","blue","black")), name = "Distance to NN")#the worse cells all scatter more widely :-(
    
    seurat_obj$selected <- seurat_obj$avgNNres < thr.bigCloneSelection
    if (returnIntermediateSeurat) return(seurat_obj)
    
    plots[["BigCloneSelection.Selected"]] <- DimPlot(seurat_obj, group.by = "selected", reduction="cloneUMAP") 
    plots[["BigCloneSelection.Selected.CellState"]] <- DimPlot(seurat_obj, group.by = "selected") 
    
    if(!is.null(trueClone)) {
      seurat_obj$true_big_clone <-  seurat_obj$csize > ncells.bigClone
      plots[["BigCloneSelection.Cotrol"]] <-  DimPlot(seurat_obj, group.by = "true_big_clone", reduction="cloneUMAP") 
      
      forAUC <-data.frame( true_small_clones = seurat_obj$csize <= ncells.bigClone, 
                           predictor = seurat_obj$avgNNres)
      forAUC <- na.omit(forAUC)
      predictions <- prediction( forAUC$predictor, forAUC$true_small_clones)
      perf <- performance(predictions, measure = "tpr", x.measure = "fpr" )
      results[["bigClone.ROC"]] <- perf
      if(!is.null(plotFolder)) {
        png(sprintf("%s/BigCloneSelection.ROC.png",plotFolder),width=640,height=480,res=150)
        plot(perf)
        abline(0,1)
        dev.off()
      }
      auc <- performance(predictions, measure = "auc")
      
      prerec <- performance(predictions, measure = "prec", x.measure = "rec" )
      results[["bigClone.PRC"]] <- prerec
      if(!is.null(plotFolder)) {
        png(sprintf("%s/BigCloneSelection.PRC.png",plotFolder),width=640,height=480,res=150)
        plot(prerec)
        dev.off()
      }
      integral <- sum(0.5 * (prerec@y.values[[1]][-1] + prerec@y.values[[1]][-length(prerec@y.values[[1]])]) * (prerec@x.values[[1]][-1] - prerec@x.values[[1]][-length(prerec@y.values[[1]])]), na.rm = T)
      results[["bigClone.stats"]] <- list(auc = auc@y.values[[1]], auprc = integral)
      
      cat(sprintf("Completed selection of big clones, AUC: %.3f , AUPRC: %.3f\n", auc@y.values[[1]], integral))
      
    } else {
      cat("Completed selection of big clones\n")
    }
    if(tuneParams) {
      cat("Tuning parameters for big clone Selection...")
      fortuning <- seurat_obj
      aucs <- mclapply(seq(30,npcs.compute,by=10), function(npcs) {
        out <- lapply(c(5,7,9,15,20), function(k) {
          out <- lapply(c(1,5,9,15,20,25,30), function(ksmooth) {
            cat(k,ksmooth,"\n")
            fortuning <- FindNeighbors(fortuning, reduction = "clonePCA", dims = 1:npcs, k.param = k, return.neighbor=T, graph.name = "clone.neighbors" )
            fortuning$avgNNdist <- apply(fortuning@neighbors$clone.neighbors@nn.dist,1,function(x) mean(x[x>0], na.rm=T))
            to.summarise <- fortuning@meta.data[, c(trueClone,batch,celltype,"avgNNdist",sprintf("nFeature_%s", methylation.assay.name))]
            
            colnames(to.summarise)[2] <- "Sample"
            colnames(to.summarise)[3] <- "celltype"
            colnames(to.summarise)[5] <- "nFeature"
            if (length(unique(to.summarise$Sample)) == 1) m <- lm(avgNNdist ~ nFeature + celltype, data = to.summarise) else m <- lm(avgNNdist ~ nFeature + Sample + celltype, data = to.summarise)
            fortuning$avgNNres <- residuals(m)
            
            if (ksmooth > 1) {
              fortuning <- FindNeighbors(fortuning, reduction = "clonePCA", dims = 1:npcs, k.param = ksmooth, return.neighbor=T, graph.name = "clone.neighbors" )
              fortuning$avgNNres <- apply(fortuning@neighbors$clone.neighbors@nn.idx,1, function(x) mean(fortuning$avgNNres[x]))
            }
            
            to.summarise$avgNNres <- fortuning$avgNNres

            csize <- table(to.summarise[,trueClone])
            fortuning$csize <- csize[as.character(fortuning@meta.data[,trueClone])] -> to.summarise$csize
            forAUC <-data.frame( true_small_clones = fortuning$csize < ncells.bigClone,
                                 predictor = fortuning$avgNNres)
            forAUC <- na.omit(forAUC)
            predictions <- prediction( forAUC$predictor,forAUC$true_small_clones)
            auc <- performance(predictions, measure = "auc")
            prerec <- performance(predictions, measure = "prec", x.measure = "rec" )
            integral <- sum(0.5 * (prerec@y.values[[1]][-1] + prerec@y.values[[1]][-length(prerec@y.values[[1]])]) * (prerec@x.values[[1]][-1] - prerec@x.values[[1]][-length(prerec@y.values[[1]])]), na.rm = T)
            data.frame(npcs = npcs, k = k, ksmooth = ksmooth, auc = auc@y.values[[1]], auprc = integral)
          })
          do.call(rbind,out)
        })
        do.call(rbind,out)
      }, mc.cores=tuneParams.cores)
      results[["bigClone.tuning"]] <- do.call(rbind, aucs)
      plots[["BigCloneSelection.tuning"]] <- qplot(x = npcs, y = auprc, data = results[["bigClone.tuning"]], color = k) + xlab("Number of PCs") + ylab("AUPRC")
      cat("Complete. Results are in slot bigClone.tuning. Rerun EPI-Clone with your preferred values.\n")
    }
    
    
    
    ### 6. cluster only the selected cells, compare to LARRY labels. ####
    cat("Now clustering clones...\n")
    for_clustering <- subset(seurat_obj, selected)
    for_clustering <- ScaleData(for_clustering, assay = methylation.assay.name, features = selected_not_protein)
    for_clustering <- RunPCA(for_clustering, assay = methylation.assay.name, features = selected_not_protein, reduction.name = "clonePCA", reduction.key = "CLONEPC_", npcs = npcs.compute)
    for_clustering <- RunUMAP(for_clustering,reduction = "clonePCA", dims = 1:npcs.Clustering, reduction.name = "cloneUMAP", reduction.key = "CLONEUMAP_" )
    
    for_clustering <- FindNeighbors(for_clustering,reduction = "clonePCA", dims = 1:npcs.Clustering, k.param =k.Clustering, graph.name = "clone_nn")
    for_clustering <- FindClusters(for_clustering, resolution=res.Clustering, graph.name = "clone_nn")
    
    plots[["finalClustering"]] <- DimPlot(for_clustering, reduction="cloneUMAP", label =T)  + NoLegend()
    
    if (!is.null(trueClone)) {
      plots[["finalClustering.groundTruth"]] <- DimPlot(for_clustering, reduction="cloneUMAP", group.by = trueClone)  + NoLegend()
      results[["clustering.adjRandIndex"]] <- adj.rand.index(as.integer(as.factor(Idents(for_clustering))), as.integer(as.factor(for_clustering@meta.data[,trueClone])))
      cat(sprintf("Completes clustering, adjusted rand index: %.3f\n", results[["clustering.adjRandIndex"]] ))
    }
    
    results[["finalSeurat"]] <- for_clustering
    if(tuneParams) {
      cat("Tuning parameters for Clustering")
      param_search <- mclapply(seq(20,npcs.compute,by=20), function(npcpar) {
        out <- lapply(c(5,9,15,25), function(kpar) {
          cat(npcpar,",",kpar,"\n")
          for_clustering <- FindNeighbors(for_clustering,reduction = "clonePCA", dims = 1:npcpar, k.param =kpar, graph.name = "clone_nn")
          out <- lapply(c(3,5,7,9), function(respar) {
            for_clustering <- FindClusters(for_clustering, resolution=respar, graph.name = "clone_nn")
            data.frame(PCs = npcpar, k = kpar, res = respar,
                       rand = adj.rand.index(as.integer(as.factor(Idents(for_clustering))), as.integer(as.factor(for_clustering@meta.data[,trueClone]))),
                       mutinf = mutinformation(Idents(for_clustering), for_clustering@meta.data[,trueClone]),
                       ncl = length(unique(Idents(for_clustering))))
          })
          do.call(rbind, out)
        })
        do.call(rbind,out)
      },mc.cores=6)
      results[["clustering.tuning"]] <- do.call(rbind, param_search)
      cat("Complete. Results are in slot clustering.tuning. Rerun EPI-Clone with your preferred values.\n")
    }
    
    
    ### 7. nicer visualizations ####
    if (!is.null(trueClone) & !is.null(plotFolder)) {
      cloneSizes <- table(for_clustering@meta.data[,trueClone])
      for_clustering$CloneForRiver <- as.factor(ifelse(cloneSizes[for_clustering@meta.data[,trueClone]] <= 10, "small clone", for_clustering@meta.data[,trueClone]))
      
      forriver <- data.frame(LARRY = factor(for_clustering$CloneForRiver), unsuper = Idents(for_clustering))
      forriver <- acast(LARRY ~unsuper, data= forriver,fun.aggregate = length)
      results[["forheatmap"]] <- forriver
      forriver_norm <- t(t(forriver) / colSums(forriver))
      
      dist.col <- as.dist(1-cor(forriver_norm))
      dist.row <- as.dist(1-cor(t(forriver_norm)))
      
      maxis <- apply(forriver_norm[!rownames(forriver_norm) %in% c("NA","small clone"),],2,max)
      rownames(forriver_norm)[!rownames(forriver_norm) %in% c("NA","small clone")] <- ""
      
      
      png(sprintf("%s/Cluster2clone.png",plotFolder),width=800,height=600,res=150)
      pheatmap(forriver_norm[,order(maxis)],cluster_cols = F,fontsize_col=6,clustering_distance_rows = dist(forriver_norm[,order(maxis)]),
               color = colorRampPalette(c("white","blue","purple","red"))(100))
      dev.off()
      
      png(sprintf("%s/Cluster2clone_cells.png",plotFolder),width=800,height=600,res=150)
      pheatmap(forriver[,order(maxis)],cluster_cols = F,fontsize_col=6,clustering_distance_rows = dist(forriver_norm[,order(maxis)]),
               color = colorRampPalette(c("white","blue","purple","red"))(100), show_rownames = F)
      dev.off()
    }
    
    results[["plots"]] <- plots
    
    if(!is.null(plotFolder)) {
      for (n in names(plots)) {
        png(sprintf("%s/%s.png", plotFolder, n), width=800,height=600,res=150)
        print(plots[[n]])
        dev.off()
      }
    }
  })
  
    return(results)
  
}
