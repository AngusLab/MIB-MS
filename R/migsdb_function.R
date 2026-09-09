#' Run MSigDB enrichment (generic collection/subcollection)
#'
#' Runs clusterProfiler::enricher() against a chosen MSigDB collection, clusters
#' redundant terms by Jaccard similarity, and generates dotplot, cnet, and STRING
#' network visualizations - same pattern as enrich_kegg/reactome_function.
#'
#' @param df Stats output, significant proteins/kinases. One column has gene symbols.
#' @param universe Background universe to use for enrichment.
#' @param gene_col Name of column with gene symbols (matches msigdbr's gene_symbol).
#' @param collection MSigDB collection, e.g. "H" (Hallmark), "C2", "C6".
#' @param subcollection MSigDB subcollection, e.g. "CP:PID", "CP:WIKIPATHWAYS". NULL for whole collection.
#' @param type String used to identify the type of comparison.
#' @param min_count Minimum number of genes needed in a term.
#' @param padj_fallback Max padj to keep if filtering leaves nothing at 0.05.
#' @param top_n_pathways Max number of pathways to keep after filtering.
#' @keywords MSigDB
#' @examples
#' msigdb_function()
#' @import dplyr stringr tidyverse tidyr ggsci ggplot2 svglite clusterProfiler igraph msigdbr enrichplot
#' @export

msigdb_function <- function(df, universe, gene_col = "Protein",
                            collection = "H", subcollection = NULL,
                            type = "Upregulated", min_count = 3,
                            padj_fallback = 0.2, top_n_pathways = 50) {

  original_wd <- getwd()
  on.exit(setwd(original_wd))
  sub_tag <- if (is.null(subcollection)) "" else paste0("_", gsub("[:/]", "-", subcollection))
  dir_name <- paste0("MSigDB_", collection, sub_tag, "_", type)
  dir.create(file.path(dir_name), showWarnings = FALSE)
  setwd(file.path(dir_name))
  cat("Running MSigDB (", collection, sub_tag, ") for: ", type, " | ", getwd(), "\n")

  msig <- msigdbr::msigdbr(species = "Homo sapiens", collection = collection,
                           subcollection = subcollection)

  if (is.null(msig) || nrow(msig) == 0) {
    cat("No MSigDB gene sets found for collection:", collection, sub_tag, ". SO sad. skipping.\n")
    return()
  }

  term2gene <- msig %>% dplyr::select(gs_name, gene_symbol)
  term2name <- msig %>% dplyr::distinct(gs_name, gs_description)

  x <- clusterProfiler::enricher(
    gene       = unique(df[[gene_col]]),
    universe   = unique(universe[[gene_col]]),
    TERM2GENE  = term2gene,
    TERM2NAME  = term2name,
    minGSSize  = 3,
    pvalueCutoff = 0.05,
    pAdjustMethod = "BH"
  )

  if (is.null(x) || nrow(x@result) == 0) {
    cat("No MSigDB results for:", type, "-", collection, sub_tag, ". It's giving low energy.\n")
    return()
  }

  msig_results <- as.data.frame(x@result) %>%
    filter(Count >= min_count) %>%
    arrange(p.adjust)

  if (nrow(msig_results) == 0) {
    cat("No MSigDB results after Count filter for:", type, "-", collection, sub_tag, ". Ope.\n")
    return()
  }
  cat("Finished running MSigDB for: ", type, "\n")

  x_eligible <- x
  x_eligible@result <- msig_results

  if (nrow(msig_results) >= 2) {
    x_simplified <- tryCatch({
      x_sim <- pairwise_termsim(x_eligible, method = "JC")
      sim <- x_sim@termsim
      d <- as.dist(1 - sim)
      hc <- hclust(d)
      cluster <- cutree(hc, h = 0.5)
      x_sim@result %>%mutate(cluster = cluster[match(Description, names(cluster))]) %>%
        group_by(cluster) %>%
        slice_min(p.adjust, n = 1, with_ties = FALSE) %>%
        ungroup()
    }, error = function(e) {
      cat("pairwise_termsim failed for ", type, ": ", e$message, ":(. \n")
      msig_results
    })
  } else {
    cat("Only", nrow(msig_results), "MSigDB term(s) for:", type, "; skipping similarity clustering.\n")
    x_simplified <- msig_results
  }
  cat("Ran MSigDB and conducted Jaccard similarity: ", type, ".\n")

  keep_ids <- head(x_simplified$ID[order(x_simplified$p.adjust)], top_n_pathways)
  keep_ids <- keep_ids[keep_ids %in% (x_simplified %>% filter(p.adjust <= padj_fallback))$ID]

  if (length(keep_ids) == 0) {
    cat("Nothing passed padj_fallback for:", type, "-", collection, sub_tag, "- skipping.\n")
    return()
  }

  x_filtered <- x
  x_filtered@result <- x@result[x@result$ID %in% keep_ids, ] %>% arrange(p.adjust)

  b <- dotplot(x_filtered, showCategory = 15,
               title = paste0(gsub("_", " ", type), " MSigDB ", collection, sub_tag, " Terms")) +
    theme_pub() +
    scale_fill_gradient(low = "red2", high = "navy")

  b$data$Description <- factor(str_wrap(b$data$Description, width = 40),
                               levels = str_wrap(levels(factor(b$data$Description)), width = 40))
  n_categories <- length(unique(b$data$Description))
  plot_height <- max(5, 3 + n_categories * 0.35)
  max_label <- max(nchar(x_filtered@result$Description), na.rm = TRUE)
  plot_width <- max(7, min(14, 5 + max_label * 0.08))

  ggsave(paste0(type, " MSigDB_", collection, sub_tag, " results.svg"), b,
         device = "svg", height = plot_height, width = plot_width)
  ggsave(paste0(type, " MSigDB_", collection, sub_tag, " results.png"), b,
         device = "png", height = plot_height, width = plot_width, dpi = 300)
  write.csv(x_filtered@result, paste0(type, "-MSigDB_", collection, sub_tag, " terms.csv"))
  cat("Made dotplot for MSigDB for: ", type, ".\n")

  if (nrow(x_filtered@result) > 1) {
    ego3 <- make_cnet_plot(x_filtered, go_type = paste0("MSigDB_", collection),
                           type = type, bio_type = "MSigDB")
    cat("Finished network analysis MSigDB for: ", type, ".\n")

    string_analysis(ego3, go_type = paste0("MSigDB_", collection), type = type)
    cat("Finished STRING network analysis MSigDB for: ", type, ".\n")
  } else {
    cat("Only", nrow(x_filtered@result), "filtered MSigDB term for:", type, "; skipping network analysis.\n")
    ego3 <- x_filtered
  }

  cat("Finished MSigDB analysis for: ", type, "-", collection, sub_tag, ".\n")
  return(ego3)
}
