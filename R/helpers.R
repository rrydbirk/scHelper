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
createEmbeddings <- function(con, n.iterations = 2, min.group.size = 10, ncomps = 40, alignment.strength = TRUE) {
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
  
  if (alignment.strength) {
    # PCA space, AS01
    con$buildGraph(space = "PCA", alignment.strength = 0.1, ncomps = ncomps)
    
    ## Embed
    con$embedGraph(method = "UMAP", embedding.name = "UMAP_PCA_AS01")
    con$embedGraph(method = "largeVis", embedding.name = "largeVis_PCA_AS01")
  }
  
  # CPCA space
  
  ## Build graph
  con$buildGraph(space = "CPCA", ncomps = ncomps)
  
  ## Embed
  con$embedGraph(method = "UMAP", embedding.name = "UMAP_CPCA")
  con$embedGraph(method = "largeVis", embedding.name = "largeVis_CPCA")
  
  if (alignment.strength) {
    # CPCA space, AS01
    con$buildGraph(space = "CPCA", alignment.strength = 0.1, ncomps = ncomps)
    
    ## Build graph
    
    ## Embed
    con$embedGraph(method = "UMAP", embedding.name = "UMAP_CPCA_AS01")
    con$embedGraph(method = "largeVis", embedding.name = "largeVis_CPCA_AS01")
  }
  
  return(con)
}

#' @export
getConditionPerCell <- function(con, conditions = c("CTRL","PD","MSA")) {
  spc <- con$getDatasetPerCell() 
  out <- spc %>% 
    as.character() %>% 
    grepl.replace(conditions) %>% 
    as.factor() %>% 
    `names<-`(spc %>% names())
  
  return(out)
}

#' @export
sget <- function(x, list.object) {
  sapply(x, '[[', list.object)
}

#' @export
lget <- function(x, list.object) {
  lapply(x, '[[', list.object)
}

#' @export
prepareObjectsForPython <- function(con, annotation, out.dir, prefix, embedding, genes.to.omit = NULL, verbose = T) {
  # Checks
  requireNamespace("sccore")
  sccore::checkPackageInstalled(c("qs","fastMatMR","conos","dplyr","magrittr")) # Need extension that's Seurat/Conos specific
  anno <- annotation[!is.na(annotation)]
  
  # Create objects
  cm.merge <- con$getJointCountMatrix(raw = T) # Consider creating parameter for this. Don't use F for Conos!
  
  anno.df <- anno %>% 
    {data.frame(cellid = names(.), annotation = unname(.))}
  
  emb <- con$embeddings[[embedding]]
  if (is.null(emb)) stop("Embedding doesn't exist.")
  
  # Create index
  idx = Reduce(intersect, list(rownames(cm.merge), anno.df$cellid, rownames(emb)))
  if (verbose) message(paste0("Using index of ",length(idx)," cells"))
  
  # Sort
  cm.merge %<>% .[match(idx, rownames(.)), ]
  anno.df %<>% .[match(idx, .$cellid), ]
  emb %<>% .[match(idx, rownames(.)), ]
  
  # Omit genes
  if (!is.null(genes.to.omit)) {
    if (all(genes.to.omit %in% colnames(cm.merge))) {
      cm.merge %<>% .[, !colnames(.) %in% genes.to.omit]
    } else {
      stop(paste0("These genes are not present: ",paste(setdiff(genes.to.omit, colnames(cm.merge)), collapse = " ")))
    }
  }
  
  # Save objects
  fastMatMR::write_fmm(cm.merge, 
                       paste0(paste(out.dir, prefix, sep = "/"), ".mtx"))
  if (verbose) message(paste0("Sparse matrix saved as ",paste0(paste(out.dir, prefix, sep = "/"), ".mtx")))
  
  object.list <- list(rownames(cm.merge), colnames(cm.merge), anno.df, emb)
  suffices <- c(".cells", ".genes", ".annotation", ".embedding")
  
  for(x in seq(4)) {
    write.table(object.list[[x]], 
                paste0(paste(out.dir, prefix, sep = "/"), suffices[[x]]), 
                sep = ",", 
                dec = ".", 
                row.names = F, 
                col.names = F)
    
    if (verbose) message(paste0("Created ", paste0(paste(out.dir, prefix, sep = "/"), suffices[[x]])))
  }
}