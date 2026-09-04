#' Reactome pathway enrichment for a set of genes
#' Runs ReactomePA's enrichPathway against a background universe, clusters similar terms,
#' and generates dotplot, cnet, STRING, and pathway view visualizations.
#'
#' @param df Stats output, significant proteins. One column is gene names, and one column includes the stats
#' @param universe Background universe to use fo enrichemnt
#' @param gene_col Name of column with Gene names
#' @param lfc_col Name of column with ranking stats
#' @param type String used to identify the type of comparison
#' @param top_n_pathways Number of pathways to display in dotplot
#' @param viewpath_n Number of pathways from Pathview
#' @param min_count Minimum number of proteins needed in a GO term
#' @param padj_fallback Max padj if nothing significant @ padj<=0.05
#' @keywords Reactome
#' @examples
#' reactome_function()
#' @import dplyr stringr tidyverse tidyr ggsci ggplot2 svglite forcats clusterProfiler igraph org.Hs.eg.db ReactomePA enrichplot GOSemSim ggraph STRINGdb pathview
#' @export

reactome_function <- function(df,universe_df,type = "Upregulated",gene_col = "Protein",
                              lfc_col  = "logFC", top_n_pathways = 50,lfc_cutoff = 0.5,
                              min_count = 3,viewpath_n = 20,padj_fallback=0.2) {

  original_wd <- getwd()
  on.exit(setwd(original_wd))
  dir_name <- paste0("Reactome_", type)
  dir.create(file.path(dir_name), showWarnings = FALSE)
  setwd(file.path(dir_name))
  message("Running Reactome for: ", type, " | ", getwd())

  #Maping Gene symbol back to EntrezID
  genes<-mapIds(org.Hs.eg.db, df[[gene_col]], "ENTREZID", "SYMBOL")
  df$entrez<- genes
  #If there are duplicate entrez mapped genes, keeps the highest ranking
  df <- df %>% filter(!is.na(entrez), !is.na(.data[[lfc_col]])) %>%
    arrange(desc(abs(.data[[lfc_col]]))) %>%
    distinct(entrez, .keep_all = TRUE)

  universe_genes<-mapIds(org.Hs.eg.db, universe_df[[gene_col]], "ENTREZID", "SYMBOL")
  universe_df$entrez<- universe_genes
  #If there are duplicate entrez mapped genes, keeps the highest ranking
  universe_df <- universe_df %>% filter(!is.na(entrez)) %>%
    distinct(entrez, .keep_all = TRUE)


  cut_off_df<- df%>%filter(abs(.data[[lfc_col]]) > lfc_cutoff)

  gene_list <- cut_off_df[[lfc_col]]
  names(gene_list) <- cut_off_df$entrez
  gene_list <- gene_list[!duplicated(names(gene_list))]
  gene_list <- sort(gene_list, decreasing = TRUE)
  de <- names(gene_list)

  if (length(de) == 0) {
    message("No genes to test for: ", type, " — skipping.")
    return()
  }

  x <- enrichPathway(gene= de,
                     universe= universe_df$entrez,
                     pvalueCutoff =0.05,
                     pAdjustMethod = "BH",
                     readable = TRUE,
                     organism = "human")


  if (is.null(x) || nrow(x@result) == 0) {
    message("No Reactome results for: ", type, " — skipping.")
    return()
  }

  react <- as.data.frame(x@result)

  react_simplified  <- tryCatch({
    x_sim<-pairwise_termsim(x, method = "JC")
    sim<-x_sim@termsim
    d<-as.dist(1-sim)
    hc<-hclust(d)
    cluster<-cutree(hc,h=0.5)# JC >= 0.5 means distance <= 0.5
    react_simplifed<-x_sim@result%>%mutate(cluster=cluster[match(Description,names(cluster))])%>%
      group_by(cluster)%>%
      slice_min(p.adjust,n=1,with_ties = F)%>%
      ungroup()

  },error = function(e) {
    message("pairwise_termsim failed for ", type, ": ", e$message)
    return(NULL)
  }
  )

  react_filtered <- react_simplified  %>% filter(p.adjust <= padj_fallback, Count >= min_count) %>%
    arrange(p.adjust) %>%
    head(top_n_pathways) %>%
    mutate(Description = factor(Description, levels = rev(unique(Description))))%>%
    separate_wider_delim(GeneRatio,delim = "/",names=c("GeneRatio_num","GeneRatio_den"))%>%
    mutate(GeneRatio=as.numeric(GeneRatio_num)/as.numeric(GeneRatio_den))%>%
    dplyr::select(-GeneRatio_num, -GeneRatio_den)%>%
    arrange(GeneRatio)%>%
    as.data.frame()


  if (nrow(react_filtered) == 0) {
    message("No Reactome results after filtering for: ", type, " — skipping.")
    return()
  }

  write.csv(react_filtered, paste0(type, "_REACTOME_terms_filtered.csv"), row.names = FALSE)
  write.csv(react, paste0(type, "_REACTOME_terms_full.csv"), row.names = FALSE)

  react_dot <- ggplot(react_filtered, aes(x = GeneRatio, y = reorder(Description, GeneRatio))) +
    geom_point(aes(size = Count, color = pvalue)) +
    scale_color_gradient(low= "red", high = "blue") +
    scale_size_continuous(name="Gene Count") +
    labs(x = "Gene Ratio",y= NULL, title = paste0(gsub("_"," ",type), " Reactome Pathways")) +
    theme_pub()

  react_dot$data$Description <- factor(str_wrap(react_dot$data$Description, width = 40),
                                       levels = str_wrap(levels(factor(react_dot$data$Description)),
                                                         width = 40))
  ## Height based on how many categories there are
  n_categories<-length(unique(react_dot$data$Description))
  plot_height <- max(5, 3 + n_categories * 0.35)

  ## Setting plot width based on size of label
  max_label <- max(nchar(as.vector(react_filtered$Description)), na.rm = TRUE)
  plot_width <- max(5, min(10, 5 + max_label * 0.08))

  ggsave(paste0(type, "_Reactome_dotplot.svg"), react_dot,
         device = "svg", width = plot_width, height = plot_height)
  ggsave(paste0(type, "_Reactome_dotplot.png"), react_dot,
         device = "png", width = plot_width, height = plot_height)


  keep_ids <- react_filtered$ID
  x_filtered <- x[x$ID %in% keep_ids, asis = TRUE]
  x_filtered@pvalueCutoff <- 0.2

  ego3<-make_cnet_plot(x_filtered,go_type="Reactome", type=type,bio_type = "Reactome")
  string_analysis(ego3,go_type="Reactome",type=type)


  pathways_to_plot <- as.vector(head(react_filtered$Description, viewpath_n))
  message("Running viewPathway on top ", length(pathways_to_plot), " pathways")
  gene_list_nodup <- gene_list[!duplicated(names(gene_list))]
  old_overlaps <- options(ggrepel.max.overlaps = Inf)
  on.exit(options(old_overlaps), add = TRUE)

  for (i in seq_along(pathways_to_plot)) {
    tryCatch({
      a <- viewPathway(pathways_to_plot[i],
                       readable   = TRUE,
                       foldChange = gene_list_nodup,
                       organism   = "human",
                       layout     = "kk") +
        theme(panel.background = element_rect(fill = "white", colour = "white"),
              plot.background  = element_rect(fill = "white", colour = "white")) +
        scale_color_gradient2(low = "blue", mid = "white", high     = "red", midpoint = 0) +
        ggtitle(paste0(pathways_to_plot[i], " — ", type)) +
        labs(color = "Log2FC")

      built <- ggplot_build(a)
      is_unobserved<- built$data[[2]]$colour == "grey50"
      is_observed <- !is_unobserved

      built$data[[2]]$colour[is_unobserved]<- scales::alpha("grey50", 0.3)
      built$data[[2]]$size[is_unobserved]<- 2
      built$data[[2]]$alpha[is_unobserved] <- 0.2
      built$data[[3]]$size[is_unobserved]<- 1
      built$data[[3]]$segment.size[is_unobserved]<- 0
      built$data[[2]]$size[is_observed]<- 3
      built$data[[3]]$fontface[is_observed]<- 2
      built$data[[3]]$point.padding[is_observed]<- 0.1
      built$data[[3]]$box.padding[is_observed]<- 0.2
      built$data[[3]]$segment.curvature[is_observed]<- 0
      built$data[[3]]$min.segment.length[is_observed]<- 0.5
      built$data[[3]]$nudge_x <- 0
      built$data[[3]]$nudge_y<- 0

      clean_name <- gsub("[[:punct:]]", "_", pathways_to_plot[i])
      ggsave(paste0(clean_name, "_", type, "_pathway.png"),
             plot   = ggplot_gtable(built),
             width  = 8, height = 8, dpi = 300)
      message("Done: ", pathways_to_plot[i])

    }, error = function(e) {
      message("viewPathway failed for: ", pathways_to_plot[i], " — ", e$message)
    })
  }

  message("Reactome complete for: ", type)
  return(invisible(react_filtered))
}
