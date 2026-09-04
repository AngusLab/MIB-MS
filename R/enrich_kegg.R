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

enrich_kegg<-function(df,universe=universe, gene_col="Protein" ,lfc_col="logFC",
                              type = paste0(unique(df$Contrast),"-all-up") ){
  # save starting directory and make new directory
  original_wd <- getwd()
  on.exit(setwd(original_wd))
  dir_name <- paste0("EnrichKEGG_", type)
  dir.create(file.path(dir_name), showWarnings = FALSE)
  setwd(file.path(dir_name))
  message("Running Enrich KEGG for: ", type, " | ", getwd())

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

  enrichgenes<- unlist(as.vector(df$entrez))

  enrichkegg <- enrichKEGG(enrichgenes,
                           keyType = "kegg",
                           minGSSize = 3,
                           organism = "hsa",
                           universe=universegenes,
                           pvalueCutoff = 0.05,
                           pAdjustMethod = "BH")

  # Check if kegg object has results before proceeding
  if (is.null(enrichkegg) || !("result" %in% slotNames(enrichkegg)) || nrow(enrichkegg@result) == 0) {
    message(paste("No enrichkegg results found for module:", type, "Skipping pathway visualization."))
    return()
  }
  enrichkegg_results<- as.data.frame(enrichkegg@result)%>%
    filter(Count >=3)%>%
    arrange(p.adjust)

  if (nrow(enrichkegg_results) == 0) {
    message(paste("No enrichkegg results found for module:", type, "Skipping pathway visualization."))
    return()
  }

  kegg_eligible<-enrichkegg[enrichkegg$ID %in% enrichkegg_results,asis=T]

  #Jaccard similarity
  keep_ids <- enrichkegg_results$ID

  kegg_simplified  <- tryCatch({
    x_sim<-pairwise_termsim(kegg_eligible, method = "JC")
    sim<-x_sim@termsim
    d<-as.dist(1-sim)
    hc<-hclust(d)
    cluster<-cutree(hc,h=0.5)# JC >= 0.5 means distance <= 0.5
    x_sim@result%>%mutate(cluster=cluster[match(Description,names(cluster))])%>%
      group_by(cluster)%>%
      slice_min(p.adjust,n=1,with_ties = F)%>%
      ungroup()

  },error = function(e) {
    message("pairwise_termsim failed for ", type, ": ", e$message)
    return(NULL)
  }
  )
  keep_ids <- kegg_simplified$ID
  kegg_raw<-enrichkegg[enrichkegg$ID %in% keep_ids,asis=T]
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

  ggsave(paste0(type," Upregulated EnrichKegg results.svg"), b,
         device = "svg", height = plot_height, width = plot_width)

  ggsave(paste0(type," Upregulated EnrichKegg results.png"), b,
         device = "png", height = plot_height, width = plot_width)

  write.csv(kegg_results, paste0(type," Enriched Kegg results.csv"))

  #Cnet and string analysis
  kegg_network<-make_cnet_plot(kegg_readable,go_type="KEGG", type=type,bio_type = "KEGG")
  string_analysis(kegg_network,go_type="KEGG", type=type)


  pathview_results <- as.data.frame(kegg_raw@result)%>% arrange(p.adjust) %>%
    head(20)

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
        str_trunc(60,ellipsis = "")

      pathview(gene.data = gene_list,
               pathway.id = pathview_results[k,"ID"],
               species = "hsa",
               out.suffix = paste0(suffix," enrichKEGG results"),
               limit = list(gene = 2),low   = list(gene = "#2166AC"),
               mid   = list(gene = "white"),high  = list(gene = "#B2182B"))
    }, error = function(e) {
      message(paste("Error processing enrichkegg ID:",
                    pathview_results[k,"ID"], "at row", k))
      message(paste("Error message:", e$message))
    },finally = {
      message("Done with this kegg analysis")
    })
  }
  print(paste0("Done with enrichKEGG."))
}
