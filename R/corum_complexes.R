#' CORUM protein complex enrichment
#'
#' Runs enrichment against curated CORUM protein complexes. Accepts the raw CORUM
#' download (one row per complex, semicolon-delimited subunits) and does the
#' organism filtering + long-format conversion internally.
#'
#' @param df Stats output, significant proteins/kinases. One column has gene symbols.
#' @param universe Background universe - detected proteome/kinome.
#' @param gene_col Name of column with gene symbols.
#' @param corum_df Loaded into package
#' @param organism Value to filter CORUM's "organism" column on. Default "Human".
#' @param type String used to identify the type of comparison.
#' @param min_count Minimum number of proteins needed in a complex.
#' @param padj_fallback Max padj to keep if nothing passes at 0.05.
#' @param top_n_pathways Max number of complexes to keep after filtering.
#' @keywords CORUM
#' @examples
#' corum_function()
#' @import dplyr stringr tidyr ggplot2 svglite clusterProfiler enrichplot
#' @export

corum_function <- function(df, universe, gene_col = "Protein",
                           corum,
                           type = "Upregulated",
                           min_count = 3,
                           padj_fallback = 0.2,
                           top_n_pathways = 50) {

  original_wd <- getwd()
  on.exit(setwd(original_wd))
  dir_name <- paste0("CORUM_", type)
  dir.create(file.path(dir_name), showWarnings = FALSE)
  setwd(file.path(dir_name))
  cat("Running CORUM complex enrichment for: ", type, " | ", getwd(), "\n")


  corum_long <- corum %>%
    dplyr::select(complex_name, subunits_gene_name) %>%
    tidyr::separate_rows(subunits_gene_name, sep = ";") %>%
    dplyr::mutate(subunits_gene_name = str_trim(subunits_gene_name)) %>%
    dplyr::filter(subunits_gene_name != "") %>%
    dplyr::distinct()

  term2gene <- corum_long %>%
    dplyr::rename(ComplexName = complex_name, gene_symbol = subunits_gene_name)

  gene_hits <- unique(df[[gene_col]])
  gene_universe <- unique(universe[[gene_col]])

  x <- clusterProfiler::enricher(
    gene          = gene_hits,
    universe      = gene_universe,
    TERM2GENE     = term2gene,
    minGSSize     = 2,
    pvalueCutoff  = 0.05,
    pAdjustMethod = "BH"
  )

  if (is.null(x) || nrow(x@result) == 0) {
    cat("No CORUM results for:", type, ".Ope .\n")
    return()
  }

  corum_results <- as.data.frame(x@result) %>%
    filter(Count >= min_count) %>%
    arrange(p.adjust)

  if (nrow(corum_results) == 0) {
    cat("No CORUM results after Count filter for:", type, ". Sad, you have no data.\n")
    return()
  }
  cat("Finished running CORUM for: ", type, "\n")

  x_eligible <- x
  x_eligible@result <- corum_results

  if (nrow(corum_results) >= 2) {
    x_simplified <- tryCatch({
      x_sim <- pairwise_termsim(x_eligible, method = "JC")
      sim <- x_sim@termsim
      d <- as.dist(1 - sim)
      hc <- hclust(d)
      cluster <- cutree(hc, h = 0.5)
      x_sim@result %>%
        mutate(cluster = cluster[match(Description, names(cluster))]) %>%
        arrange(cluster, p.adjust)
    }, error = function(e) {
      cat("pairwise_termsim failed for ", type, ": ", e$message, "\n")
      corum_results %>% mutate(cluster = NA)
    })
  } else {
    cat("Only", nrow(corum_results), "CORUM complex(es) for:", type, "; skipping similarity clustering.\n")
    x_simplified <- corum_results %>% mutate(cluster = NA)
  }
  cat("Ran CORUM and conducted Jaccard similarity: ", type, ".\n")

  x_ranked <- x_simplified %>% filter(p.adjust <= padj_fallback) %>% arrange(p.adjust)
  keep_ids <- head(x_ranked$ID, top_n_pathways)

  if (length(keep_ids) == 0) {
    cat("Nothing passed padj_fallback for CORUM:", type, ". This is REALLY not significnat.\n")
    return()
  }

  x_filtered <- x
  x_filtered@result <- x@result[x@result$ID %in% keep_ids, ] %>% arrange(p.adjust)

  b <- dotplot(x_filtered, showCategory = 15,
               title = paste0(gsub("_", " ", type), " CORUM Complex Enrichment")) +
    theme_pub() +
    scale_fill_gradient(low = "red2", high = "navy")

  b$data$Description <- factor(str_wrap(b$data$Description, width = 40),
                               levels = str_wrap(levels(factor(b$data$Description)), width = 40))
  n_categories <- length(unique(b$data$Description))
  plot_height <- max(5, 3 + n_categories * 0.35)
  max_label <- max(nchar(x_filtered@result$Description), na.rm = TRUE)
  plot_width <- max(7, min(14, 5 + max_label * 0.08))

  ggsave(paste0(type, "_CORUM_results.svg"), b, device = "svg",
         height = plot_height, width = plot_width)
  ggsave(paste0(type, "_CORUM_results.png"), b, device = "png",
         height = plot_height, width = plot_width, dpi = 300)
  write.csv(x_filtered@result, paste0(type, "_CORUM_terms.csv"), row.names = FALSE)
  cat("Made dotplot for CORUM for: ", type, ".\n")

  if (nrow(x_filtered@result) > 1) {
    ego3 <- make_cnet_plot(x_filtered, go_type = "CORUM", type = type, bio_type = "CORUM")
    cat("Finished network analysis CORUM for: ", type, ".\n")

    string_analysis(ego3, go_type = "CORUM", type = type)
    cat("Finished STRING network analysis CORUM for: ", type, ".\n")
  } else {
    cat("Only", nrow(x_filtered@result), "filtered CORUM complex for:", type,
        "; skipping network analysis.\n")
    ego3 <- x_filtered
  }

  cat("Finished CORUM analysis for: ", type, ".\n")
  return(ego3)
}
