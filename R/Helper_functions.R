#' Internal helper for enrichment network plots
#'
#' @keywords internal

make_cnet_plot<-function(bio_result,go_type="BP", type="helper",bio_type="GO"){

  # term-gene network
  p_cnet <- cnetplot(bio_result,
                     showCategory = 10,
                     layout = igraph::layout_with_fr,
                     color_category = "#2C7FB8",
                     color_item = "grey10",
                     color_edge = "grey75",
                     size_category = 1.5,
                     size_item = 1,
                     size_edge = 0.2) +
    theme_void() +
    theme(legend.position = "right",
          plot.title = element_text(face = "bold"),
          text = element_text(size = 11),
          panel.background = element_rect(fill = "white", colour = "white"),
          plot.background = element_rect(fill = "white", colour = "white"))+
    ggtitle(paste0(go_type," ",gsub("_"," ",type)," cnetplot"))

  for (i in seq_along(p_cnet$layers)) {
    p_cnet$layers[[i]]$aes_params$alpha <- 0.6
  }

  p_cnet
  ggsave(paste0(go_type,"_","Cnetplot_of_",type,".png"), p_cnet, width = 8, height = 6, dpi = 300)
  ggsave(paste0(go_type,"_","Cnetplot_of_",type,".svg"), p_cnet, width = 8, height = 6, device = "svg")

  # term-term similarity network - takes simiplified GO terms and calculates similarity
  ## Makes plots showing how much they're connected
  if(bio_type=="GO"){
    d <- GOSemSim::godata('org.Hs.eg.db', ont=go_type)
    ego2 <- pairwise_termsim(bio_result, method="Wang", semData = d)

  }else{
    ego2<-pairwise_termsim(bio_result,method="JC")
  }

  if (nrow(ego2@result) < 2) {
    message("Too few terms for emapplot. Skipping.")
  } else {
    #Emap- essentially clusters of the common terms
    p_emap <- emapplot(ego2, showCategory = 10, layout = "kk") +
      theme_void()+
      theme(legend.position = "right",
            plot.title = element_text(size = 12, face = "bold", hjust = 0),
            text = element_text(size = 11),
            plot.background = element_rect(fill = "white", colour = "white")) +
      ggtitle(paste0(go_type, " ", gsub("_"," ",type), " emapplot"))+
      scale_color_gradient (low = "red2", high = "navy")

    ggsave(paste0(go_type, "_", type, "_emapplot.png"), p_emap,
           width = 10, height = 8, dpi = 300)
    ggsave(paste0(go_type, "_", type, "_emapplot.svg"), p_emap, device = "svg",
           width = 10, height = 8, dpi = 300)
    #Tree plot with key words at the end.
    p_tree <- treeplot(ego2, showCategory = 30,
                       group_color=ggsci::pal_jco("default")(10),
                       cluster.params = list(n = 5,label_words_n = 2,label_format = 25)) +
      scale_color_gradient(low  = "red2",
                           high = "navy",
                           name = "p.adjust",
                           oob = scales::squish)+
      theme_void()+
      theme(plot.title= element_text(face = "bold", hjust = 0.5),
            plot.background  = element_rect(fill = "white", colour = "white"),
            panel.background = element_rect(fill = "white", colour = "white"))+
      ggtitle(paste0(go_type," ",gsub("_"," ",type), " treeplot"))

    ggsave(paste0(go_type, "_", type, "_treeplot.png"), p_tree,
           width = 14, height = 8, dpi = 300)
    ggsave(paste0(go_type, "_", type, "_treeplot.svg"), p_tree, device = "svg",
           width = 14, height = 8, dpi = 300)

    #Simplified dotplot
    b<-dotplot(ego2, showCategory = 10,
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
    max_label <- max(nchar(ego2@result$Description), na.rm = TRUE)
    plot_width <- max(7, min(14, 5 + max_label * 0.08))

    ggsave(paste0(type," Simplified",go_type, "results.svg"), b,
           device = "svg", height = plot_height, width = plot_width)

    ggsave(paste0(type," Simplified",go_type, "results.png"), b,
           device = "png", height = plot_height, width = plot_width)

  }

  return(ego2)
}


string_analysis<-function(simplifed_term,go_type="BP", type="helper"){
  # STRING interaction network for hit genes
  genes <- unique(unlist(strsplit(simplifed_term@result$geneID, "/")))

  string_db <- STRINGdb$new(version = "12", species = 9606, score_threshold = 400)
  mapped <- string_db$map(data.frame(gene = genes), "gene", removeUnmappedRows = TRUE)
  hits <- mapped$STRING_id
  ppi <- string_db$get_interactions(hits)

  # Map STRING IDs back to gene symbols
  id_map <- mapped[, c("STRING_id", "gene")]
  ppi$from_symbol<- id_map$gene[match(ppi$from, id_map$STRING_id)]
  ppi$to_symbol <- id_map$gene[match(ppi$to, id_map$STRING_id)]
  ppi2 <- ppi[!is.na(ppi$from_symbol) & !is.na(ppi$to_symbol), ]

  if(nrow(ppi2)<2){
    message("Too few mapped PPIs")
    return()
  }

  # Build network graph from ppi df
  g <- igraph::graph_from_data_frame(ppi2[, c("from_symbol", "to_symbol")],directed = FALSE)

  # Basic STRING-like plot
  a<-ggraph(g, layout = "fr") +
    geom_edge_link(alpha = 0.2, color = "grey60") +
    geom_node_point(size = 4, color = "#2C7FB8") +
    geom_node_text(aes(label = name), repel = TRUE, size = 4) +
    theme_void() +
    theme(plot.title= element_text(face = "bold", hjust = 0.5),
          panel.background = element_rect(fill = "white", colour = "white"),
          plot.background = element_rect(fill = "white", colour = "white"))+
    ggtitle(paste0(go_type," ",gsub("_"," ",type)," STRING interactions"))
  a

  ggsave(paste0(go_type,"_",type,"_basic String plot.png"),
         a, width = 10, height = 8, dpi = 300)
  ggsave(paste0(go_type,"_",type,"_basic String plot.svg"),
         a, width = 10, height = 8, device="svg")

  # Gene -> GO term mapping
  gene2term <- stack(setNames(strsplit(simplifed_term@result$geneID, "/"),simplifed_term@result$Description))
  colnames(gene2term) <- c("gene", "term")
  # Assign one term per gene by lowest p adjustment
  gene2term_full <- merge(gene2term, simplifed_term@result[, c("Description", "p.adjust")],by.x = "term",
                          by.y = "Description")
  gene_term <- gene2term_full%>%group_by(gene)%>%
    slice_min(p.adjust,n=1, with_ties = F)%>%
    ungroup()

  #Adding Complex info to the String network
  igraph::V(g)$Complex <- gene_term$term[match(igraph::V(g)$name, gene_term$gene)]
  igraph::V(g)$Complex <- as.character(igraph::V(g)$Complex)

  # fill missing term from neighbors
  na_nodes <- igraph::V(g)[is.na(igraph::V(g)$Complex)]
  print(na_nodes)

  for (v in na_nodes) {
    nb <- igraph::neighbors(g, v)
    nb_terms <- as.character(V(g)$Complex[nb])
    nb_terms <- na.omit(nb_terms)

    if (length(nb_terms) == 0) next

    tab<-sort(table(nb_terms), decreasing=TRUE)
    if (tab[1]/sum(tab)>=0.7){
      igraph::V(g)$Complex[v]<-names(tab)[1]
    }
  }
  igraph::V(g)$Complex[is.na(igraph::V(g)$Complex)] <- "Unassigned"
  levs <- unique(igraph::V(g)$Complex)
  levs_no_unassigned <- setdiff(levs, "Unassigned")

  # Plotting complex info ontop of basic STRING graph
  extended_cols <- colorRampPalette(ggsci::pal_jco("default")(10))(length(levs_no_unassigned))
  colors <- setNames(extended_cols, levs_no_unassigned)
  colors["Unassigned"] <- "grey70"
  igraph::V(g)$Complex <- factor(igraph::V(g)$Complex, levels = names(colors))

  # Plot colored by term
  b_string<-ggraph(g, layout = "fr") +
    geom_edge_link(alpha = 0.2, color = "grey85") +
    geom_node_point(aes(color = Complex), size = 4) +
    geom_node_text(aes(label = name), repel = TRUE, size = 4) +
    scale_color_manual(values = colors) +
    theme_void() +
    theme(plot.title= element_text(face = "bold", hjust = 0.5),
          panel.background = element_rect(fill = "white", colour = "white"),
          plot.background = element_rect(fill = "white", colour = "white"))+
    ggtitle(paste0(go_type," ",gsub("_"," ",type)," Common protein Function"))

  ggsave(paste0(go_type,"_",type,"_String plot with complexes.png"), b_string,
         width = 10, height = 8, dpi = 300)

  ggsave(paste0(go_type,"_",type,"_String plot with complexes.svg"), b_string,
         width = 10, height = 8, device="svg")

  # Cluster modules
  set.seed(123)
  cl <- igraph::cluster_leiden(g, objective_function = "modularity",resolution = 1)
  igraph::V(g)$module <- factor(cl$membership)


  extended_cols <- colorRampPalette(ggsci::pal_jco("default")(10))(length(unique(V(g)$module)))
  c<-ggraph(g, layout = "fr") +
    geom_edge_link(alpha = 0.15, color = "grey80") +
    geom_node_point(aes(color = module), size = 4) +
    geom_node_text(aes(label = name), repel = TRUE, size = 4) +
    scale_color_manual(values = extended_cols)+
    theme_void() +
    theme(plot.title= element_text(face = "bold", hjust = 0.5),
          panel.background = element_rect(fill = "white", colour = "white"),
          plot.background = element_rect(fill = "white", colour = "white"))+
    ggtitle(paste0(go_type," ",gsub("_"," ",type)," proteins with leiden clustering"))

  ggsave(paste0(go_type,"_",type,"_String plot with leiden clustered modules.png"),c)
  ggsave(paste0(go_type,"_",type,"_String plot with leiden clustered modules.svg"),c, device = "svg")

}

theme_pub <- function(base_size = 11) {
  theme_bw(base_size = base_size) +
    theme(
      axis.line         = element_line(linewidth = 0.4, color = "black"),
      axis.ticks        = element_line(linewidth = 0.4, color = "black"),
      axis.text         = element_text(size = base_size - 1, color = "black"),
      axis.title        = element_text(size = base_size, face = "bold"),
      strip.background  = element_blank(),
      strip.text        = element_text(size = base_size, face = "bold"),
      legend.key.size   = unit(0.4, "cm"),
      legend.text       = element_text(size = base_size - 2),
      legend.title      = element_text(size = base_size - 1, face = "bold"),
      plot.title        = element_text(size = base_size + 1, face = "bold", hjust = 0),
      plot.subtitle     = element_text(size = base_size - 1, color = "grey40", hjust = 0),
      panel.spacing     = unit(0.6, "lines"),
      plot.margin       = margin(8, 8, 8, 8)
    )
}
