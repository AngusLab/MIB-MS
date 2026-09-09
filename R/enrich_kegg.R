#' Enrich KEGG pathways for a set of genes
#'
#' Runs clusterProfiler's enrichKEGG against a background universe, clusters similar terms
#' by Jaccard similarity, and generates dotplot, cnet, STRING, and pathview visualizations.
#'
#' @param df Stats output, only the significant ones. One column is gene names, and one column includes the stats
#' @param universe Background universe to use fo enrichemnt
#' @param gene_col Name of column with Gene names
#' @param lfc_col Name of column with ranking stats
#' @param type String used to identify the type of comparison
#' @keywords Enrich KEGG
#' @examples
#' enrich_kegg()
#' @import dplyr stringr tidyverse tidyr ggsci ggplot2 svglite forcats clusterProfiler igraph org.Hs.eg.db enrichplot DOSE pathview
#' @export

enrich_kegg<-function(df,universe, gene_col="Protein" ,lfc_col="logFC",
                              type = paste0(unique(df$Contrast),"-all-up") ){
  # save starting directory and make new directory
  original_wd <- getwd()
  on.exit(setwd(original_wd))
  dir_name <- paste0("EnrichKEGG_", type)
  dir.create(file.path(dir_name), showWarnings = FALSE)
  setwd(file.path(dir_name))
  cat("Running Enrich KEGG for: ", type, " | ", getwd(),"\n")

  genes<-mapIds(org.Hs.eg.db, df[[gene_col]], "ENTREZID", "SYMBOL")
  df$entrez<- genes
  #If there are duplicate entrez mapped genes, keeps the highest ranking
  df <- df %>% filter(!is.na(entrez), !is.na(.data[[lfc_col]])) %>%
    arrange(desc(abs(.data[[lfc_col]]))) %>%
    distinct(entrez, .keep_all = TRUE)

  #converting universe to entrez
  universe_genes<-mapIds(org.Hs.eg.db, universe[[gene_col]], "ENTREZID", "SYMBOL")
  universe$entrez<- universe_genes
  #If there are duplicate entrez mapped genes, outputs vector of entrez IDs
  universegenes <- universe %>% filter(!is.na(entrez)) %>%
    distinct(entrez)%>%pull(entrez)

  enrichkegg <- enrichKEGG(df$entrez,
                           keyType = "kegg",
                           minGSSize = 3,
                           organism = "hsa",
                           universe=universegenes,
                           pvalueCutoff = 0.05,
                           pAdjustMethod = "BH")

  # Check if kegg object has results before proceeding
  if (is.null(enrichkegg) || !("result" %in% slotNames(enrichkegg)) || nrow(enrichkegg@result) == 0) {
    cat(paste("No enrichkegg results found for module:", type, ". No data for you: BOOM.TOASTED.\n"))
    return()
  }
  enrichkegg_results<- as.data.frame(enrichkegg@result)%>%
    filter(Count >=3)%>%
    arrange(p.adjust)

  if (nrow(enrichkegg_results) == 0) {
    cat(paste("No enrichkegg results found for module:", type, "Skipping pathway visualization.\n"))
    return()
  }
  cat("Finished running Enrich KEGG for: ", type, ". Now onto Jaccard. ", "\n")

  kegg_eligible<-enrichkegg
  kegg_eligible@result<-enrichkegg_results

  #Jaccard similarity

  if (nrow(enrichkegg_results) >= 2) {
    kegg_simplified <- tryCatch({
      x_sim <- pairwise_termsim(kegg_eligible, method = "JC")
      sim <- x_sim@termsim
      d <- as.dist(1 - sim)
      hc <- hclust(d)
      cluster <- cutree(hc, h = 0.5)  # JC >= 0.5 means distance <= 0.5
      x_sim@result %>% mutate(cluster = cluster[match(Description, names(cluster))]) %>%
        arrange(cluster, p.adjust)
    }, error = function(e) {
      cat("pairwise_termsim failed for ", type, ": ", e$message, "\n")
      enrichkegg_results %>% mutate(cluster = NA)
    })
  } else {
    cat("Only", nrow(enrichkegg_results), "enriched KEGG term for:", type, "; skipping similarity clustering.\n")
    kegg_simplified <- enrichkegg_results
  }

  keep_ids <- kegg_simplified$ID
  kegg_raw<-enrichkegg
  kegg_raw@result<-enrichkegg@result[enrichkegg@result$ID %in% keep_ids,]

  kegg_readable <- DOSE::setReadable(kegg_raw,OrgDb = org.Hs.eg.db,keyType = "ENTREZID")
  kegg_results <- as.data.frame(kegg_readable@result) %>%
    arrange(p.adjust)

  b<-dotplot(kegg_readable, showCategory = 10,
             title = paste0(gsub("_"," ",type)," Enriched Pathways"))+
    theme_pub()+
    scale_fill_gradient (low = "red2", high = "navy")

  b$data$Description <- factor(str_wrap(b$data$Description, width = 40),
                               levels = str_wrap(levels(factor(b$data$Description)),
                                                 width = 40))
  ## Height based on how many categories there are
  n_categories<-length(unique(b$data$Description))
  plot_height <- max(5, 3 + n_categories * 0.35)

  ## Setting plot width based on size of label
  max_label <- max(nchar(kegg_results$Description), na.rm = TRUE)
  plot_width <- max(7, min(14, 5 + max_label * 0.08))

  ggsave(paste0(type," EnrichKegg results.svg"), b,
         device = "svg", height = plot_height, width = plot_width)

  ggsave(paste0(type," EnrichKegg results.png"), b,
         device = "png", height = plot_height, width = plot_width)

  write.csv(kegg_results, paste0(type," Enriched Kegg results.csv"))
  cat("Made dotplot for Enrich Kegg for: ", type, ".\n")

  #Cnet and string analysis
  if(nrow(kegg_results)>1){
  kegg_network<-make_cnet_plot(kegg_readable,go_type="KEGG", type=type,bio_type = "KEGG")
  cat("Made network plots for Enrich Kegg for: ", type, ".\n")

  string_analysis(kegg_network,go_type="KEGG", type=type)
  cat("STRING analysis finished for Enrich Kegg for: ", type, ".\n")

  }

  pathview_results <- as.data.frame(kegg_raw@result)%>% arrange(p.adjust) %>%
    head(20)

  cat("Pathways going to Pathview:", nrow(pathview_results), "\n")
  cat("Pathview output directory:", getwd(), "\n")

  for(k in 1:nrow(pathview_results)){
    tryCatch({
      genes.2<- as.vector(pathview_results[k,"geneID"])
      genes.2<- strsplit(genes.2, split = "/", fixed = T)
      genes.2<- unlist(genes.2)
      #Grabbing lfc from initial df
      gene.list<- df%>%filter(entrez %in% genes.2)%>%
        dplyr::select(all_of(c("entrez", lfc_col)))
      gene_list<-setNames(gene.list[[lfc_col]],gene.list[["entrez"]])
      suffix<-pathview_results[k,"Description"] |> str_replace_all("[^A-Za-z0-9 _-]", "") |>
        str_squish()|>
        str_trunc(20,ellipsis = "")
      cat( "Running Pathview",k, "/", nrow(pathview_results), ":", pathview_results$ID[k],"\n")
      pathview(gene.data = gene_list,
               pathway.id = pathview_results[k,"ID"],
               kegg.dir = getwd(),
               species = "hsa",
               out.suffix = paste0(suffix," enrichKEGG results"),
               limit = list(gene = 2),low   = list(gene = "#2166AC"),
               mid   = list(gene = "white"),high  = list(gene = "#B2182B"))
      cat("Finished Pathview:", pathview_results$ID[k], "\n")
    }, error = function(e) {
      cat(paste("Error processing enrichkegg ID:",
                    pathview_results[k,"ID"], "at row", k),".\n")
      cat(paste("Error message:", e$message),"\n")
    },finally = {
      cat("Done with this kegg analysis.\n")
    })
  }
  cat("Done with enrichKEGG.\n")
}
