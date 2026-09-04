#' Enrich GO terms for a set of genes
#' Runs clusterProfiler's enrichGO against a background universe, filters and simplifies
#' the results, and generates dotplot, cnet, and STRING network visualizations.
#'
#' @param df Stats output, significant proteins. One column is gene names, and one column includes the stats
#' @param universe Background universe to use fo enrichemnt
#' @param gene_col Name of column with Gene names
#' @param lfc_col Name of column with ranking stats
#' @param type String used to identify the type of comparison
#' @param go_type Ontology being tested. Either BP, CC, or MF
#' @param min_count Minimum number of proteins needed in a GO term
#' @keywords Enrich GO
#' @examples
#' enrich_go()
#' @import dplyr stringr tidyverse tidyr ggsci ggplot2 svglite forcats clusterProfiler igraph org.Hs.eg.db enrichplot pathview
#' @export

enrich_go<-function(df,universe,gene_col="Protein",lfc_col="logFC" ,
                            type="Upregulated_all",go_type="BP", min_count=3){

  # save starting directory and make new directory
  original_wd <- getwd()
  on.exit(setwd(original_wd))
  dir_name <- paste0("EnrichGO_",go_type,"_" ,type)
  dir.create(file.path(dir_name), showWarnings = FALSE)
  setwd(file.path(dir_name))
  message("Running Enrich GO for: ", type,": ",go_type, " | ", getwd())


  genes<-mapIds(org.Hs.eg.db, df[[gene_col]], "ENTREZID", "SYMBOL")
  df$entrez<- genes
  #If there are duplicate entrez mapped genes, outputs vector of entrez IDs
  df <- df %>% filter(!is.na(entrez)) %>%
    arrange(desc(abs(.data[[lfc_col]]))) %>%
    distinct(entrez, .keep_all = TRUE)

  universe_genes<-mapIds(org.Hs.eg.db, universe[[gene_col]], "ENTREZID", "SYMBOL")
  universe$entrez<- universe_genes
  #If there are duplicate entrez mapped genes, outputs vector of entrez IDs
  universegenes <- universe %>% filter(!is.na(entrez)) %>%
    distinct(entrez)%>%pull(entrez)

  go_terms_results<- enrichGO(gene= df$entrez,
                              universe = universegenes,
                              OrgDb = org.Hs.eg.db,
                              ont = go_type,
                              pAdjustMethod = "BH",
                              keyType = "ENTREZID",
                              pvalueCutoff = 0.05,
                              readable = TRUE)

  if (is.null(go_terms_results) || nrow(go_terms_results@result) == 0) {
    message("No GO results for: ", type," ",go_type ," — skipping.")
    return()
  }

  up_go_term_df<- as.data.frame(go_terms_results@result)
  boring_terms <- c("kinase activity", "phosphorylation",
                    "protein phosphorylation", "ATP binding",
                    "transferase activity", "phosphotransferase activity")
  up_go_term_df <- up_go_term_df %>%
    filter(!grepl(paste(boring_terms, collapse = "|"),
                  Description, ignore.case = TRUE))%>%
    filter(Count>=min_count)

  go_terms_results@result<-up_go_term_df

  if (nrow(up_go_term_df) == 0) {
    message(paste("No enrichgo results found for:", type, "Skipping visualization."))
    return()
  }
  b<-dotplot(go_terms_results, showCategory = 15,
             title = paste0(gsub("_"," ",type), " GO ", go_type, " Terms"))+
    theme_pub()+
    scale_fill_gradient (low = "red2", high = "navy")


  b$data$Description <- factor(str_wrap(b$data$Description, width = 40),
                               levels = str_wrap(levels(factor(b$data$Description)),
                                                 width = 40))
  ## Height based on how many categories there are
  n_categories<-length(unique(b$data$Description))
  plot_height <- max(5, 3 + n_categories * 0.35)

  ## Setting plot width based on size of label
  max_label <- max(nchar(up_go_term_df$Description), na.rm = TRUE)
  plot_width <- max(7, min(14, 5 + max_label * 0.08))

  ggsave(paste0(type," EnrichGO results.svg"), b, device = "svg", height = plot_height, width=plot_width)
  ggsave(paste0(type," EnrichGO results.png"), b, device = "png",
         height = plot_height, width=plot_width, dpi=300)
  write.csv(up_go_term_df, paste0(type, "-",go_type ,"GO Terms.csv"))


  # Simplify GO terms using clusterprofiler - removes redundat GO terms
  ego2 <- clusterProfiler::simplify(go_terms_results,
                                    cutoff = 0.5, by = "p.adjust", select_fun = min)


  ego3<-make_cnet_plot(ego2,go_type=go_type, type=type,bio_type="GO")

  string_analysis(ego3,go_type=go_type,type=type)

  print(paste0("Done with enrichGO for: ", type," | ", go_type))

  return(ego3)



}
