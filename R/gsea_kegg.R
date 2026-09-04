#' This function performs GSEA KEGG from clsuter profiler and outputs graphs
#'
#'
#' @param df Stats output. One column is gene names, and one column includes the stats
#' @param gene_col Name of column with Gene names
#' @param lfc_col Name of column with ranking stats
#' @param type String used to identify the type of comparison
#' @param n_perm Number of permutations when running GSEA
#' @keywords GSEAKegg
#' @examples
#' gsea_kegg()
#' @import dplyr stringr tidyverse tidyr ggsci ggplot2 svglite forcats clusterProfiler igraph org.Hs.eg.db enrichplot pathview
#' @export

gsea_kegg<-function(df,gene_col="Protein" ,lfc_col="logFC",
                            type =paste0(unique(df$Contrast),"-all"),
                            n_perm=10000){
  # save starting directory and make new directory
  original_wd <- getwd()
  on.exit(setwd(original_wd))
  dir_name <- paste0("GSEAKEGG_", type)
  dir.create(file.path(dir_name), showWarnings = FALSE)
  setwd(file.path(dir_name))
  message("Running GSEA KEGG for: ", type, " | ", getwd())

  #Maping Gene symbol back to EntrezID
  genes<-mapIds(org.Hs.eg.db, df[[gene_col]], "ENTREZID", "SYMBOL")
  df$entrez<- genes
  #If there are duplicate entrez mapped genes, keeps the highest ranking
  df <- df %>% filter(!is.na(entrez), !is.na(.data[[lfc_col]])) %>%
    arrange(desc(abs(.data[[lfc_col]]))) %>%
    distinct(entrez, .keep_all = TRUE)

  #Using ranking column
  gene_list<- df[[lfc_col]]
  names(gene_list)<- df$entrez
  gene_list<- sort(gene_list, decreasing = T)

  kegg <- gseKEGG(geneList = gene_list,
                  keyType = "kegg",
                  pAdjustMethod = "BH",
                  minGSSize = 3,
                  organism = "hsa",
                  nPermSimple=n_perm)

  if (is.null(kegg) || !("result" %in% slotNames(kegg)) || nrow(kegg@result) == 0) {
    message(paste("No significant gseKEGG results found for module:", type, ". Skipping gseKEGG-related steps."))
    return()
  }

  #Make Kegg results a df
  kegg_results<- as.data.frame(kegg@result)
  kegg_results$hsa.id <- rownames(kegg_results)
  #Make dotplot

  a<-dotplot(kegg, showCategory = 15,
             title = paste0(gsub("_"," ",type)," GSEA KEGG Enriched Pathways") ,
             split=".sign")+
    facet_grid(.~.sign)+
    theme_pub() +
    scale_fill_gradient(low = "red2", high = "navy")

  a$data$Description <- factor(str_wrap(a$data$Description, width = 40),
                               levels = str_wrap(levels(factor(a$data$Description)),
                                                 width = 40))
  #Width based on how many categories there are
  n_categories<-length(unique(a$data$Description))
  plot_height <- max(5, 3 + n_categories * 0.35)
  ## Setting plot width based on size of label
  max_label <- max(nchar(kegg_results$Description), na.rm = TRUE)
  plot_width <- max(7, min(14, 5 + max_label * 0.08))


  ggsave(paste0(type," GSEKegg results.svg"), a, device = "svg",
         width=plot_width, height =plot_height)
  ggsave(paste0(type," GSEKegg results.png"), a, device = "png",
         width=plot_width, height =plot_height)

  pathview_results <- as.data.frame(kegg_results@result)%>% arrange(p.adjust) %>%
    head(20)

  #Make pathview obj
  for(j in 1:nrow(kegg_results)){
    tryCatch({
      genes.2<- as.vector(pathview_results[j,c("core_enrichment")])
      genes.2<- strsplit(genes.2, split = "/", fixed = T)
      genes.2<- unlist(genes.2)
      #Grabbing lfc from initial df
      gene.list<- df%>%filter(entrez %in% genes.2)%>%
        dplyr::select(all_of(c("entrez", lfc_col)))
      gene_list<-setNames(gene.list[[lfc_col]],gene.list[["entrez"]])
      suffix<-pathview_results[j,"Description"] |> str_replace_all("[^A-Za-z0-9 _-]", "") |>
        str_squish()|>
        str_trunc(60,ellipsis = "")

      pathview(gene.data = gene_list, pathway.id = pathview_results[j,"ID"],
               species = "hsa",
               out.suffix = paste0(suffix," gseKEGG results"),
               limit = list(gene = 2),low   = list(gene = "#2166AC"),
               mid   = list(gene = "white"),high  = list(gene = "#B2182B"))
    }, error = function(e) {
      message(paste("Error processing pathway ID:", pathview_results[j,"ID"], "at row"))
      message(paste("Error message:", e$message))
    }, finally = {
      message("Done with this kegg analysis")
    })
  }


  write.csv(kegg_results,paste0(type, " GSEAKegg results.csv"))
  return(kegg_results)
}
