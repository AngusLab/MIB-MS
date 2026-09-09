#' Kinase group/family enrichment
#'
#' Tests whether hit kinases are overrepresented in specific kinome groups/families
#' relative to the detected kinome background. Unlike GO/KEGG/MSigDB, groups are
#' mutually exclusive by construction, so no Jaccard/redundancy-reduction step is used.
#'
#' @param df Stats output, significant kinases. One column has gene symbols.
#' @param universe Background universe - should be your FULL DETECTED KINOME, not the
#'   whole proteome/genome, since this tests overrepresentation within the kinome.
#' @param gene_col Name of column with gene symbols.
#' @param kinome_annotation Lookup table with columns matching gene_col's symbol type
#'   plus a "group" column (and optionally "family" for finer resolution). This can be found jsut by loading the human.kinome
#' @param type String used to identify the type of comparison.
#' @param level Either "group" (8 major kinome groups) or "family" (finer resolution).
#' @param min_count Minimum number of kinases needed in a group/family to keep it.
#' @keywords Kinase Group Enrichment
#' @examples
#' kinase_group()
#' @import dplyr stringr ggplot2 svglite clusterProfiler enrichplot
#' @export

kinase_group <- function(df, universe, gene_col = "Protein",
                         kinome_annotation,
                                  type = "Upregulated",
                                  level = "Group",
                                  min_count = 2) {

  original_wd <- getwd()
  on.exit(setwd(original_wd))
  dir_name <- paste0("KinaseGroup_", level, "_", type)
  dir.create(file.path(dir_name), showWarnings = FALSE)
  setwd(file.path(dir_name))
  cat("Running Kinase", level, "enrichment for: ", type, " | ", getwd(), "\n")


  term2gene <- kinome_annotation %>%
    dplyr::select(all_of(c(level, "Gene"))) %>%
    dplyr::filter(!is.na(.data[[level]]), .data[[level]] != "") %>%
    dplyr::distinct()

  gene_hits <- unique(df[[gene_col]])
  gene_universe <- unique(universe[[gene_col]])

  x <- clusterProfiler::enricher(
    gene          = gene_hits,
    universe      = gene_universe,
    TERM2GENE     = term2gene,
    minGSSize     = 3,
    pvalueCutoff  = 1,
    pAdjustMethod = "BH"
  )

  if (is.null(x) || nrow(x@result) == 0) {
    cat("No kinase", level, "enrichment results for:", type, ". No data for you: BOOM.TOASTED.\n")
    return()
  }

  group_results <- as.data.frame(x@result) %>%
    filter(Count >= min_count) %>%
    arrange(p.adjust)

  if (nrow(group_results) == 0) {
    cat("No kinase", level, "results after Count filter for:", type, ". Ope.\n")
    return()
  }
  cat("Finished kinase", level, "enrichment for: ", type, "\n")

  x_filtered <- x
  x_filtered@result <- group_results

  write.csv(group_results, paste0(type, "_kinase_", level, "_enrichment.csv"), row.names = FALSE)

  # Bar chart
  group_results_plot <- group_results %>%
    tidyr::separate_wider_delim(GeneRatio, delim = "/", names = c("num", "den")) %>%
    mutate(GeneRatio = as.numeric(num) / as.numeric(den)) %>%
    dplyr::select(-num, -den) %>%
    arrange(GeneRatio) %>%
    mutate(Description = factor(Description, levels = unique(Description)))

  bar <- ggplot(group_results_plot, aes(x = GeneRatio, y = Description, fill = p.adjust)) +
    geom_col() +
    geom_text(aes(label = Count), hjust = -0.2, size = 3.2) +
    scale_fill_gradient(low = "red2", high = "navy") +
    labs(x = "Gene Ratio", y = NULL,
         title = paste0(gsub("_", " ", type), " Kinase ", tools::toTitleCase(level), " Enrichment")) +
    theme_pub() +
    expand_limits(x = max(group_results_plot$GeneRatio) * 1.15)

  n_categories <- nrow(group_results_plot)
  plot_height <- max(3, 1.5 + n_categories * 0.4)

  ggsave(paste0(type, "_kinase_", level, "_barplot.svg"), bar,
         device = "svg", width = 8, height = plot_height)
  ggsave(paste0(type, "_kinase_", level, "_barplot.png"), bar,
         device = "png", width = 8, height = plot_height, dpi = 300)
  cat("Made kinase", level, "barplot for: ", type, ".\n")

  # cnetplot: which specific kinases map to which group
  if (nrow(x_filtered@result) > 1) {
    ego3 <- tryCatch({
      make_cnet_plot(x_filtered, go_type = paste0("Kinase_", level), type = type,
                     bio_type = "KinaseGroup")
    }, error = function(e) {
      cat("Cnetplot failed for kinase", level, ":", type, "-", e$message, "\n")
      x_filtered
    })
    cat("Finished network plot for kinase", level, ": ", type, ".\n")
  } else {
    cat("Only", nrow(x_filtered@result), level, "passed filtering for:", type,
        "; skipping cnetplot.\n")
    ego3 <- x_filtered
  }

  cat("Finished kinase", level, "enrichment for: ", type, ".\n")
  return(ego3)
}
