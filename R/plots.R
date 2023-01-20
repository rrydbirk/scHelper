#' @import magrittr conos
#' @importFrom cowplot plot_grid
NULL

#' @export
plotEmbeddingOverview <- function(con, dot.size = 0.5) {
  # Generate factors
  spc <- con$getDatasetPerCell()
  
  cpc <- getConditionPerCell(con)
  
  # Get embedding plots
  embeddings <- con$embeddings %>% 
    names() %>% 
    lapply(\(emb) con$plotGraph(embedding = emb, 
                                title = paste0(emb,", clusters"), 
                                size = dot.size)
    )
  
  # Get sample plots
  samples <- con$embeddings %>% 
    names() %>% 
    lapply(\(emb) {
      con$plotGraph(embedding = emb,
                    groups = spc, 
                    mark.groups = FALSE, 
                    show.legend = FALSE,
                    size = dot.size,
                    title = paste0(emb,", samples"))
    })
  
  # Get condition plots
  conditions <- con$embeddings %>% 
    names() %>% 
    lapply(\(emb) {
      con$plotGraph(embedding = emb,
                    groups = cpc, 
                    mark.groups = FALSE, 
                    show.legend = FALSE,
                    size = dot.size,
                    title = paste0(emb,", conditions"))
    })
  
  # Plot
  seq(1:length(con$embeddings)) %>% 
    lapply(\(x) {
      list(embeddings[[x]], 
           samples[[x]], 
           conditions[[x]]) %>% 
        plot_grid(plotlist = ., ncol = 3)
    }) %>% 
    plot_grid(plotlist = ., ncol = 1)
}