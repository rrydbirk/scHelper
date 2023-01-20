#' @import magrittr dplyr conos
NULL

#' @export
mitoFraction <- function(con, species="human", missing = NA) {
  if (species == "human") symb <- "MT-" else if (species == "mouse") symb <- "mt-" else stop("Species must either be 'human' or 'mouse'.")
  
  out <- con$samples %>% 
    lapply(\(sample) {
      mt.no <- sample$counts %>% 
        colnames() %>% 
        grep(symb, .) %>% 
        length()
      
      if (mt.no == 1) {
        sample$counts[,grep(symb, colnames(sample$counts))] / sparseMatrixStats::rowSums2(sample$counts)
      } else if (mt.no == 0) {
        rep(missing, sample$counts %>% nrow()) %>% 
          setNames(sample$counts %>% rownames())
      } else {
        sparseMatrixStats::rowSums2(sample$counts[,grep(symb, colnames(sample$counts))]) / sparseMatrixStats::rowSums2(sample$counts)
      }
    }) %>% 
    Reduce(c, .) %>% 
    setNames(con$samples %>% lapply(`[[`, "counts") %>% sapply(rownames) %>% Reduce(c, .))
  
  if (any(is.na(out))) warning(paste0("Some samples did not contain mitochondrial expression. Values for cells in those samples set to ",missing,"."))
  
  return(out)
}

#' @export
addEmbeddingP2Web <- function(p2, con, embedding=NULL, name="UMAP") {
  if(is.null(embedding)) embedding <- con$embedding
  
  if(identical(dim(p2$originalP2object$embeddings$PCA[[1]]),dim(embedding))) {
    p2$originalP2object$embeddings$PCA[[name]] <- embedding
    return(p2)
  } else {
    stop("The embedding dimensions of the p2.web object and the input object are not identical.")
  }
}

#' @export
embedUMAP <- function(con,
                      min.dist=0.01,
                      spread=15,
                      min.prob.lower=1e-7,
                      method=leiden.community,
                      resolution=1,
                      min.group.size=25,
                      n.iterations = 1) {
  message("Creating UMAP embedding...")
  con$embedGraph(method="UMAP", 
                 min.dist=min.dist, 
                 spread=spread,
                 min.prob.lower=min.prob.lower)
  
  message("Estimating clusters...")
  con$findCommunities(method=leiden.community, resolution=resolution, min.group.size=min.group.size, n.iterations=n.iterations)
  
  return(con)
}

#' @export
buildConosGraph <- function(con,
                            k.conos=15, 
                            k.self=15, 
                            space='PCA', 
                            ncomps=40,
                            n.odgenes=2e3,
                            matching.method='mNN', 
                            metric='angular', 
                            score.component.variance=T,
                            alignment.strength=0,
                            min.dist=0.01, 
                            spread=15,
                            min.prob.lower=1e-3,
                            resolution=1,
                            min.group.size=25,
                            n.iterations = 1) {
  message("Building graph...")
  con$buildGraph(k=k.conos, 
                 k.self=k.self, 
                 space=space, 
                 ncomps=ncomps, 
                 n.odgenes=n.odgenes, 
                 matching.method=matching.method, 
                 metric=metric, 
                 verbose=T, 
                 score.component.variance=score.component.variance,
                 alignment.strength=alignment.strength)
  
  embedUMAP(con=con,
            min.dist=min.dist,
            spread=spread,
            min.prob.lower=min.prob.lower,
            method=leiden.community,
            resolution=resolution,
            min.group.size=min.group.size,
            n.iterations=n.iterations)
  
  return(con)
}

#' @export
quickConos <- function(cms, 
                       sample.names,
                       n.cores.p2,
                       n.cores.con,
                       n.odgenes=3e3, 
                       nPcs = 50, 
                       k.p2 = 30, 
                       perplexity = 50, 
                       log.scale = TRUE, 
                       trim = 10, 
                       keep.genes = NULL, 
                       min.cells.per.gene = 3, 
                       min.transcripts.per.cell = 200, 
                       get.largevis = F, 
                       get.tsne = F, 
                       make.geneknn = F,
                       k.conos=15, 
                       k.self=30, 
                       space='PCA', 
                       ncomps=40, 
                       matching.method='mNN', 
                       metric='angular', 
                       score.component.variance=T,
                       alignment.strength=0,
                       min.dist=0.01, 
                       spread=15,
                       n.iterations = 1) {
  if(length(cms)==length(sample.names)) {
    if(any(is.na(sample.names))) stop("Names contains NAs")
  }
  
  if(any(duplicated(unlist(lapply(cms,colnames))))) {
    cms <- sample.names %>% 
      lapply(\(sample) {
        cm <- cms[[sample]]
        colnames(cm) %<>% {paste0(sample,"!!",.)}
        return(cm)
      }) %>% 
      setNames(sample.names)
    
    message("Performing P2 processing...")
    panel.preprocessed <- lapply(cms, function(x) basicP2proc(x, n.cores = n.cores.p2,
                                                              n.odgenes = n.odgenes, 
                                                              nPcs = nPcs,
                                                              k = k.p2, 
                                                              perplexity = perplexity, 
                                                              log.scale = log.scale, 
                                                              trim = trim, 
                                                              keep.genes = keep.genes, 
                                                              min.cells.per.gene = min.cells.per.gene, 
                                                              min.transcripts.per.cell = min.transcripts.per.cell, 
                                                              get.largevis = get.largevis, 
                                                              get.tsne = get.tsne, 
                                                              make.geneknn = make.geneknn))
    
    names(panel.preprocessed) = sample.names
    con <- Conos$new(panel.preprocessed, n.cores=n.cores.con)
    
    con <- buildConosGraph(con=con,
                           k.conos=k.conos, 
                           k.self=k.self, 
                           space=space, 
                           ncomps=ncomps, 
                           n.odgenes=n.odgenes, 
                           matching.method=matching.method, 
                           metric=metric, 
                           score.component.variance=score.component.variance,
                           alignment.strength=alignment.strength,
                           min.dist=min.dist, 
                           spread=spread,
                           n.iterations=n.iterations)
    
    return(list(con=con, panel.preprocessed=panel.preprocessed))
  } else {
    stop("Sample names must match number of count matrices.")
  }
  
}

#' @export
collapseAnnotation <- function(anno, label) {
  anno %<>% factor
  idx <- grepl(label,levels(anno))
  cat(paste0("Collapsing ",sum(idx)," labels containing '",label,"' in their name into one label.\n"))
  levels(anno)[idx] <- c(label)
  anno %<>% factor
  return(anno)
}

#' @export
getConosDepth <- function(con) {
  lapply(con$samples, function(d) d$depth) %>% unlist %>% setNames(.,(strsplit(names(.), ".", T) %>% 
                                                                        sapply(function(d) d[2])))
}

#' @export
getConosCluster <- function(con, name="leiden") {
  con$clusters[[name]]$groups
}

#' @export
sortCPDB <- function(path) {
  pval_path <- paste0(path,"pvalues.txt")
  sigmean_path <- paste0(path,"significant_means.txt")
  message(paste0("Looking for the following files:\n",pval_path,"\n",sigmean_path))
  pval_full <- read.table(pval_path, sep="\t", header=T)
  sigmean_full <- read.table(sigmean_path, sep="\t", header=T)
  message("Input is ",nrow(pval_full)," interactions.")
  
  #Remove doublet interactions
  before <- nrow(pval_full)
  pval_full %<>% .[match(unique(.$id_cp_interaction), .$id_cp_interaction),]
  sigmean_full %<>% .[match(unique(.$id_cp_interaction), .$id_cp_interaction),]
  
  if(before > nrow(pval_full)) message("Removed ",before - nrow(pval_full)," doublet interactions.")
  
  sigmean <- sigmean_full[,13:ncol(sigmean_full)]
  rownames(sigmean) <- sigmean_full$id_cp_interaction
  
  #Check consistency between data
  message(paste0(length(setdiff(pval_full$id_cp_interaction,rownames(sigmean)))," pairs mismatching between matrices."))
  
  # Interactions W/O significant pairs
  rs <- rowSums(sigmean, na.rm = T)
  row_idx <- rs[rs == 0] %>% 
    names()
  
  cs <- colSums(sigmean, na.rm = T)
  col_idx <- cs[cs == 0] %>% 
    names()
  
  pval_full %<>% 
    .[,!colnames(.) %in% col_idx] %>% 
    .[!.$id_cp_interaction %in% row_idx,]
  sigmean_full %<>% 
    .[,!colnames(.) %in% col_idx] %>% 
    .[!.$id_cp_interaction %in% row_idx,]
  
  # Save cleaned tables
  message(paste0("Saving tables with ",nrow(sigmean)," interaction pairs."))
  write.table(pval_full, paste0(path,"pvalues_clean.txt"), sep="\t", col.names=T, row.names=F)
  sigmean_full[is.na(sigmean_full)] <- ""
  write.table(sigmean_full, paste0(path,"significant_means_clean.txt"), sep="\t", col.names=T, row.names=F)
  message("All done!")
}

#' @export
renameAnnotation <- function(annotation, old, new) {
  annotation %<>% factor()
  
  levels(annotation)[levels(annotation) %in% old] <- new
  
  return(annotation)
}

#' @export
dotSize <- function(size, alpha=1) {
  ggplot2::guides(colour = guide_legend(override.aes = list(size=size,
                                                            alpha=alpha)))
}

#' @export
checkDims <- function(cm, con) {
  cat("Dimensions of cm : ",paste((dim(cm)), collapse=" "),"\n")
  cat("Dimensions of con: ",paste((dim(con$embedding)), collapse=" "),"\n")
  
  if(dim(cm)[2]!=dim(con$embedding)[1])
    stop("Dimensions don't match.")
  
  message("All OK!")
}

#' @export
grepl.replace <- function(x, patterns, result = NULL) {
  if(is.null(result)) result <- patterns
  if(length(patterns) != length(result)) stop("'patterns' and 'result' must have equal lengths.")
  
  for(i in 1:length(patterns)) {
    x[grepl(patterns[i], x)] <- result[i]
  }
  return(x)
}

#' @export
dotPlot2 <- function (markers, count.matrix, cell.groups, marker.colour = "black", 
                      cluster.colour = "black", xlab = "Marker", ylab = "Cluster", 
                      n.cores = 1, text.angle = 45, gene.order = NULL, cols = c("blue", 
                                                                                "red"), col.min = -2.5, col.max = 2.5, dot.min = 0, 
                      dot.scale = 6, scale.by = "radius", scale.min = NA, scale.max = NA, 
                      verbose = TRUE, ...) 
{
  scale.func <- switch(scale.by, size = scale_size, radius = scale_radius, 
                       stop("'scale.by' must be either 'size' or 'radius'"))
  if (!is.character(markers)) {
    stop("'markers' must be a character vector.")
  }
  missing.markers <- setdiff(markers, colnames(count.matrix))
  if (length(missing.markers) > 0) {
    message("Not all markers are in 'count.matrix'. The following are missing: ", 
            paste(missing.markers, collapse = " "))
    stop("Please update 'markers'.")
  }
  marker.table <- table(markers)
  if (sum(marker.table > 1) != 0) {
    message("The following genes are present more than once in 'markers': ", 
            paste(names(marker.table[marker.table > 1]), collapse = " "), 
            " These genes will only be plotted at first instace. Consider revising. ")
  }
  if (verbose) {
    message("Extracting gene expression... ")
  }
  if (inherits(cell.groups, "factor")) {
    tryCatch({
      if (verbose) {
        message("Treating 'cell.groups' as a factor.")
      }
      cell.groups %<>% as.factor()
    }, error = function(e) stop("Could not convert 'cell.groups' to a factor\n", 
                                e))
  }
  p.df <- plapply(markers, function(g) data.frame(Expr = count.matrix[names(cell.groups), 
                                                                      g], Type = cell.groups, Gene = g), n.cores = n.cores, 
                  progress = verbose, ...) %>% Reduce(rbind, .)
  if (is.logical(gene.order) && gene.order) {
    gene.order <- unique(markers)
  }
  
  if (!is.null(gene.order)) {
    p.df %<>% dplyr::mutate(Gene = factor(as.character(.data$Gene), 
                                          levels = gene.order))
  }
  if (verbose) {
    message("Calculating expression distributions... ")
  }
  data.plot <- levels(cell.groups) %>% plapply(function(t) {
    markers %>% lapply(function(g) {
      df <- p.df %>% dplyr::filter(.data$Type == t, .data$Gene == 
                                     g)
      pct.exp <- sum(df$Expr > 0)/dim(df)[1] * 100
      avg.exp <- mean(df$Expr[df$Expr > 0])
      res <- data.frame(gene = g, pct.exp = pct.exp, avg.exp = avg.exp)
      return(res)
    }) %>% Reduce(rbind, .)
  }, n.cores = n.cores, progress = verbose, ...) %>% stats::setNames(levels(cell.groups)) %>% 
    dplyr::bind_rows(., .id = "cluster")
  data.plot$cluster %<>% factor(., levels = rev(unique(.)))
  data.plot %<>% dplyr::arrange(.data$gene)
  data.plot$avg.exp.scaled <- data.plot$gene %>% unique %>% 
    sapply(function(g) {
      data.plot %>% .[.$gene == g, "avg.exp"] %>% scale %>% 
        setMinMax(min = col.min, max = col.max)
    }) %>% unlist %>% as.numeric
  data.plot$pct.exp[data.plot$pct.exp < dot.min] <- NA
  cluster.colour %<>% rev
  if (!is.null(gene.order)) data.plot %<>% mutate(gene = gene %>% factor(levels = gene.order))
  plot <- ggplot(data.plot, aes(gene, cluster)) + 
    geom_point(aes_string(size = "pct.exp", color = "avg.exp.scaled")) + 
    scale.func(range = c(0, dot.scale), limits = c(scale.min, 
                                                   scale.max)) + theme(axis.text.x = element_text(angle = text.angle, 
                                                                                                  hjust = 1, colour = marker.colour), axis.text.y = element_text(colour = cluster.colour), 
                                                                       panel.background = element_rect(fill = "white", colour = "black", 
                                                                                                       size = 1, linetype = "solid"), panel.grid.major = element_blank(), 
                                                                       panel.grid.minor = element_blank()) + guides(size = guide_legend(title = "Percent expressed"), 
                                                                                                                    color = guide_colorbar(title = "Average expression")) + 
    labs(x = xlab, y = ylab) + scale_color_gradient(low = cols[1], 
                                                    high = cols[2])
  return(plot)
}

#' @export
createEmbeddings <- function(con, n.iterations = 1, min.group.size = 10, ncomps = 100) {
  con$embedding <- NULL
  con$embeddings <- NULL
  
  # PCA space
  
  ## Build graph
  con$buildGraph(space = "PCA", ncomps = ncomps)
  
  ## Identify clusters
  if (length(con$clusters) == 0) con$findCommunities(n.iterations = n.iterations, min.group.size = min.group.size)
  
  ## Embed
  con$embedGraph(method = "UMAP", embedding.name = "UMAP_PCA")
  con$embedGraph(method = "largeVis", embedding.name = "largeVis_PCA")
  
  # PCA space, AS01
  con$buildGraph(space = "PCA", alignment.strength = 0.1, ncomps = ncomps)
  
  ## Embed
  con$embedGraph(method = "UMAP", embedding.name = "UMAP_PCA_AS01")
  con$embedGraph(method = "largeVis", embedding.name = "largeVis_PCA_AS01")
  
  # CPCA space
  
  ## Build graph
  con$buildGraph(space = "CPCA", ncomps = ncomps)
  
  ## Embed
  con$embedGraph(method = "UMAP", embedding.name = "UMAP_CPCA")
  con$embedGraph(method = "largeVis", embedding.name = "largeVis_CPCA")
  
  # CPCA space, AS01
  con$buildGraph(space = "CPCA", alignment.strength = 0.1, ncomps = ncomps)
  
  ## Build graph
  
  ## Embed
  con$embedGraph(method = "UMAP", embedding.name = "UMAP_CPCA_AS01")
  con$embedGraph(method = "largeVis", embedding.name = "largeVis_CPCA_AS01")
  
  return(con)
}

getConditionPerCell <- function(con, conditions = c("CTRL","PD","MSA")) {
  spc <- con$getDatasetPerCell() 
  out <- spc %>% 
    as.character() %>% 
    grepl.replace(conditions) %>% 
    as.factor() %>% 
    `names<-`(spc %>% names())
  
  return(out)
}