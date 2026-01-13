#' This function cleans up and analyzes the DIA output from DIANN, adding a pseudocount
#'
#'
#' @param df Unique protein groups matrix from DIANN
#' @param kinases Either the human or mouse kinome spreadsheet
#' @param peptides report.parquet from DIANN
#' @param metadata File with column titled "Sample.ID" with the column names of the abundance values in the protein groups file, and a column titled "Treatment" with the corresponding treatment replicate
#' @param directory Output folder
#' @param unique_df Either "protein" for pg.matrix or "gene" for unique_gene.matrix
#' @keywords DIANN
#' @examples
#' diann.cleanup()
#' @import dplyr fuzzyjoin stringr tibble arsenal tidyverse tidyr data.table sjmisc ggpubr ggsci ggplot2 svglite rstatix pheatmap arrow
#' @export

diann.cleanup<- function(df, kinases, metadata, directory, unique_df, peptide ){
  ifelse(!dir.exists(directory), dir.create(directory), paste0(directory," folder exists already"))
  setwd(directory)
  ifelse(!dir.exists("Heatmaps"), dir.create("Heatmaps"), "Heatmaps folder exists already")
  ifelse(!dir.exists("Volcano_plots"), dir.create("Volcano_plots"), "Volcano_plots folder exists already")
  ifelse(!dir.exists("PCA"), dir.create("PCA"), "PCA folder exists already")
  ifelse(!dir.exists("Results"), dir.create("Results"), "Results folder exists already")
  ifelse(!dir.exists("QC"), dir.create("QC"), "QC folder exists already")

  #Post-DIANN quant Cleanup

  log_message <- function(message, file = "Processing_log.txt") {
    timestamp <- format(Sys.time(), "[%Y-%m-%d %H:%M:%S]")
    cat(timestamp, message, "\n", file = file, append = TRUE)
    print(message)
  }
  #filtering out precursor rows that map back to multiple protein IDs

  unique<- peptide%>%
    mutate(cleaned_genes=str_trim(Genes))%>%
    filter(!is.na(cleaned_genes), cleaned_genes!="")%>%
    filter(!str_detect(cleaned_genes,";"))

  # subsetting based on input file, either unique genes matrix or protein groups
  if(unique_df=="protein"){
    log_message(paste0("Input: ", unique_df))
    #making a table with number of peptide that map back to a gene
    unique_counts<- unique%>%
      distinct(cleaned_genes, Stripped.Sequence)%>%
      count(cleaned_genes, name="unique_gene_peptides")
    #filtering out proteins that have only 1 mapped peptide

    keep<- unique_counts%>%
      filter(unique_gene_peptides>1)%>%
      pull(cleaned_genes)

    protein<- df[df$Genes %in% keep,]
    protein<- protein%>%
      filter(!str_detect(Protein.Group, ";"))

    rownames(protein)<-protein$Genes; protein<- protein[,-c(1,2,4)]
  }else{
    print("Not the Protein DF")
  }

  if(unique_df=="gene"){
    log_message(paste0("Input: ", unique_df))
    #making a table with number of peptide that map back to a gene
    unique_counts<- unique%>%
      distinct(cleaned_genes, Stripped.Sequence)%>%
      count(cleaned_genes, name="unique_gene_peptides")

    #filtering out proteins that have only 1 mapped peptide

    keep<- unique_counts%>%
      filter(unique_gene_peptides>1)%>%
      pull(cleaned_genes)

    protein<- df[df$Genes %in% keep,]

    protein<-protein%>%
      group_by(Genes)%>%
      summarise(across(everything(), median, na.rm=TRUE))
    rownames(protein)<-protein$Genes
  }else{
    print("Not the Gene matrix")
  }


  #Printing stats
  y<- as.character(length(unique(df$Genes)))
  x<- as.character(length(unique(protein$Genes)))
  log_message(paste0("Number of proteins kept: ", x, "/",y))

  kinase.peptides1<- df[df$Genes%in% human.kinome$Gene,]
  kinase.peptides<- protein[protein$Genes%in% human.kinome$Gene,]
  y<- as.character(length(unique(kinase.peptides1$Genes)))
  x<- as.character(length(unique(kinase.peptides$Genes)))
  log_message(paste0("Number of kinases kept: ", x, "/",y))

  #Match Sample ID to treatment
  protein<-protein %>%
    rename_with(~deframe(metadata)[.x], .cols = metadata$Sample.ID) %>%
    dplyr::select(Genes, any_of(metadata$Treatment))

  log_message("Matched to sample ID")

  histo<- protein%>%
    pivot_longer(!Genes)%>%
    mutate(value=log2(value))%>%
    mutate(value= tidyr::replace_na(value, 0))%>%
    mutate(type ="Pre-median")%>%
    mutate(name=sub("Z.", "", name))%>%
    mutate(name=sub("ZZ,","",name))

  #Log2 transform
  x<- protein[,-c(1)]
  x<- log2(x)
  x[sapply(x, is.infinite)] <- NA
  df.2<- protein[,c(1)]
  df.3<- cbind(df.2, x)
  colnames(df.3)[1]<- "Genes"
  log_message("Log transformed")


  # Median center
  medians<- c()
  for(i in 2:ncol(df.3)){
    x<- median(df.3[,i], na.rm = T)
    medians<- append(medians, x)
  }
  max.med<- max(medians)

  med.norm<- list()

  for(i in 2:ncol(df.3)){
    x<- median(df.3[,i], na.rm = T)
    u<- ((df.3[,i]/x)*max.med)
    med.norm[[i-1]]<- flatten(as.data.frame(u))
  }

  df.4 <- do.call(rbind, med.norm)
  df.4<- as.data.frame(t(df.4))
  df.4<- sapply(df.4, function(x) as.numeric(as.character(x)))
  df.4<- as.data.frame(df.4)
  df.4$Genes<- df.3$Genes
  df.4<- df.4[,c(ncol(df.4), 1:ncol(df.4)-1)]; colnames(df.4)<- colnames(df.3)
  log_message("Median normalized")

  histo2<- df.4%>%
    pivot_longer(!Genes)%>%
    mutate(value= tidyr::replace_na(value, 0))%>%
    mutate(type ="Post-median")%>%
    mutate(name=sub("Z.", "", name))%>%
    mutate(name=sub("ZZ,","",name))

  histogram<- rbind(histo, histo2)

  b<- ggplot(histogram, aes(x=(value),fill = type))+
    geom_histogram(color="#e9ecef", alpha=0.6, position = "identity")+
    scale_fill_manual(values=c("#69b3a2", "#404080")) +
    ggtitle("Log2(Relative Protein Abundances Median Normalization)")+
    xlab("log2(Protein Abundances)")+
    theme_bw()+
    facet_wrap(~name, ncol = 3)

  b

  ggsave(filename = file.path("QC","Global abundances histogram.svg"), b, device= "svg", width = 10, height = 6)
  ggsave(filename = file.path("QC","Global pre-abundances histogram.pdf"), b, device= "pdf", width = 10, height = 6)


  c<- ggplot(histogram%>%filter(type=="Pre-median")%>%filter(value!=0), aes(x=(value),y = name))+
    geom_jitter(color="darkgray", alpha=0.9)+
    ggtitle("Log2(Relative Protein Abundances Pre-Median Normalization)")+
    xlab("log2(Protein Abundances)")+
    theme_bw()+
    stat_summary(fun.x = median, fun.xmin = median, fun.xmax = median,
                 geom = "crossbar", width = 0.5)

  c

  ggsave(filename = file.path("QC","Global abundances pre median normalization.svg"), c, device= "svg", width = 10, height = 6)
  ggsave(filename = file.path("QC","Global abundances pre median normalization.pdf"),c, device= "pdf", width = 10, height = 6)

  d<- ggplot(histogram%>%filter(type=="Post-median")%>%filter(value!=0), aes(x=(value),y = name))+
    geom_jitter(color="darkgray", alpha=0.9)+
    ggtitle("Log2(Relative Protein Abundances Post-Median Normalization)")+
    xlab("log2(Protein Abundances)")+
    theme_bw()+
    stat_summary(fun.x = median, fun.xmin = median, fun.xmax = median,
                 geom = "crossbar", width = 0.5)


  d
  ggsave(filename = file.path("QC","Global abundances post median normalization.svg"), d, device= "svg", width = 10, height = 6)
  ggsave(filename = file.path("QC","Global abundances post median normalization.pdf"),d, device= "pdf", width = 10, height = 6)

  #Filter for columns where at least one group has 3 values
  rownames(df.4)<- df.4[,1]
  df.4<- as.data.frame(t(df.4))
  df.4<- df.4[-c(1),]
  df.4$type<- rownames(df.4)
  df.4<- df.4[,c(ncol(df.4), 1:ncol(df.4)-1)]
  df.4$type<- sub("-.*", '', df.4$type)
  df.4$type<- sub("_.*", '', df.4$type)

  x<- as.vector(NULL)
  g<- as.data.frame(table(df.4$type))
  g<-max(g$Freq)


  if(g>2){
    for(i in 2:ncol(df.4)){
      a<- df.4[(1)]
      b<- df.4[i]
      df.5<- cbind(a,b)
      y<- colnames(df.5)[2]
      colnames(df.5)[2]<- "Protein"
      df.6<- df.5 %>%
        group_by(type)%>%
        summarise(total_non_na = sum(!is.na(Protein)))

      if(any(df.6$total_non_na >2)){
        x<- append(x, y)
      }
    }
  }else if(g == 2){
    for(i in 2:ncol(df.4)){
      a<- df.4[(1)]
      b<- df.4[i]
      df.5<- cbind(a,b)
      y<- colnames(df.5)[2]
      colnames(df.5)[2]<- "Protein"
      df.6<- df.5 %>%
        group_by(type)%>%
        summarise(total_non_na = sum(!is.na(Protein)))

      if(any(df.6$total_non_na >1)){
        x<- append(x, y)
      }
    }
  } else{
    for(i in 2:ncol(df.4)){
      a<- df.4[(1)]
      b<- df.4[i]
      df.5<- cbind(a,b)
      y<- colnames(df.5)[2]
      colnames(df.5)[2]<- "Protein"
      df.6<- df.5 %>%
        group_by(type)%>%
        summarise(total_non_na = sum(!is.na(Protein)))

      if(any(df.6$total_non_na==1)){
        x<- append(x, y)
      }
    }}


  log_message("Filtered for proteins with enouch sample representation")
  print("Filtered for proteins with enouch sample representation")



  #Impute 0s with LOD -1
  df.4<- df.4[colnames(df.4)%in% x]
  df.5<-data.frame(sapply(df.4, function(x) as.numeric(as.character(x))))
  log_message(paste0("The minimum value before adding pseudocount:", min(df.5, na.rm = T)))
  rownames(df.5)<- rownames(df.4)

  missing_mat <- is.na(df.5) * 1
  svg(filename = file.path("QC","Heatmap of missing values.svg"), width = 10, height = 10)
  pheatmap(missing_mat,
           cluster_rows=TRUE,
           cluster_cols=TRUE,
           border_color = "black",
           show_colnames = FALSE,)

  dev.off()
  min_value<- min(df.5, na.rm = T)-1
  df.5[is.na(df.5)]<-min_value
  log_message(paste0("The minimum value after adding pseudocount:", min(df.5, na.rm = T)))

  histo3<- as.data.frame(t(df.5)); histo3$Genes<- rownames(histo3)
  histo3<- histo3%>%
    pivot_longer(!Genes)%>%
    mutate(value= tidyr::replace_na(value, 0))%>%
    mutate(type ="Post-imputation")%>%
    mutate(name = sub("Z.","",name))%>%
    mutate(name = sub("ZZ.","",name))

  histogram<- rbind(histogram, histo3)

  e<- ggplot(histogram, aes(x=(value),fill = type))+
    geom_histogram(color="#e9ecef", alpha=0.6)+
    scale_fill_manual(values=c("red","#69b3a2", "#404080")) +
    ggtitle("Log2(Relative Protein Abundances) Post Filtering and Imputation")+
    xlab("log2(Protein Abundances)")+
    theme_bw()+
    facet_wrap(~name, ncol = 3)

  e
  ggsave(filename = file.path("QC","Global abundances post imputation and filtering.svg"), e, device= "svg", width = 10, height = 6)
  ggsave(filename = file.path("QC","Global abundances post imputation and filtering.pdf"),e, device= "pdf", width = 10, height = 6)

  imputed_df.2<- df.5
  imputed_df.2$Treatment<- rownames(imputed_df.2)
  imputed_df.2$Treatment<- sub("-.*", '', imputed_df.2$Treatment)
  imputed_df.2$Treatment<- sub("_.*", '', imputed_df.2$Treatment)
  imputed_df.2<- imputed_df.2[,c(ncol(imputed_df.2), 1:ncol(imputed_df.2)-1)]

  log_message("Imputed")

  new_df<- imputed_df.2
  new_df<- new_df[,c(-1)]
  cols.num<- colnames(new_df)
  new_df[cols.num]<- sapply(new_df[cols.num], as.numeric)
  df_pca<- prcomp(new_df)


  bc.pca.var<- df_pca$sdev^2
  bc.pca.var<- round(bc.pca.var/sum(bc.pca.var)*100,1)
  df.pca_data<- data.frame(sample= rownames(df_pca$x),
                           x=df_pca$x[,1],
                           y = df_pca$x[,2])

  df.pca_data$type<- sub("_.*", "", df.pca_data$sample)
  df.pca_data$type<- sub("^ZZ.","", df.pca_data$type)
  df.pca_data$type<- sub("^Z.", "", df.pca_data$type)
  df.pca_data<-df.pca_data %>%
    group_by(type) %>%
    mutate(Rep = row_number()) %>%
    ungroup()

  df.pca_data$Rep<- as.character(as.numeric(df.pca_data$Rep))

  a<-ggplot(data=df.pca_data,aes(label=sample,x=x,y=y, color = type, shape = Rep)) +
    ggtitle("Global PCA Plot") +
    geom_point(aes( size = 4)) +
    xlab(paste("PC1: ",bc.pca.var[1],"%",sep="")) +
    ylab(paste("PC2: ",bc.pca.var[2],"%",sep=""))+
    scale_color_jco()+
    theme_classic()+
    theme(legend.text = element_text(size = 10))+
    guides(color=guide_legend(override.aes = list(size = 5)),
           shape = guide_legend(override.aes = list(size=5)))
  a

  ggsave(filename = file.path("PCA","All proteins PCA plot.svg"), a, device= "svg", width = 10, height = 6)
  ggsave(filename = file.path("PCA","All proteins plot.pdf"), a, device= "pdf", width = 10, height = 6)

  write.csv(imputed_df.2, file = file.path("Results","Relative Global Protein Abundance.csv"))


  #Matching to either the human or mouse kinome
  imputed_df.2<- as.data.frame(t(imputed_df.2))
  imputed_df.3<- imputed_df.2[-c(1),]
  subset_kinase<- imputed_df.3[rownames(imputed_df.3) %in% kinases$Gene,]
  imputed_df.2<- as.data.frame(t(subset_kinase))

  log_message("Matched to kinome")
  imputed_df.3<-data.frame(sapply(imputed_df.2[,2:ncol(imputed_df.2)], function(x) as.numeric(as.character(x))))
  imputed_df.3$Treatment<- rownames(imputed_df.2)
  rownames(imputed_df.3)<- rownames(imputed_df.2)
  imputed_df.3$Treatment<- sub("-.*", '', imputed_df.3$Treatment)
  imputed_df.3$Treatment<- sub("_.*", '', imputed_df.3$Treatment)
  imputed_df.3<- imputed_df.3[,c(ncol(imputed_df.3), 1:ncol(imputed_df.3)-1)]


  write.csv(imputed_df.3, file = file.path("Results","Relative Kinase Protein Abundance.csv"))

  #ANOVA
  DF<- imputed_df.3
  kinases<-colnames(DF)
  #choose which columns to drop, typically this is the first column
  kinases<-kinases[-c(1)]
  #collapse to all kinases separated by + sign
  kinases<-paste(kinases,collapse="+")

  ##Comparative Stats - All Groups ANOVA
  df<-summary(tableby(as.formula(paste('Treatment~',kinases)),
                      numeric.stats=c("mean"),
                      stats.labels=list(mean="Mean"),
                      numeric.simplify=TRUE,
                      cat.simplify=TRUE,
                      numeric.test="anova",
                      total=FALSE,
                      digits.p=20,
                      data=DF), text = NULL)
  df<-as.data.frame(df)
  write.csv(df, file=file.path("Results","ANOVA.csv"))

  log_message("Completed Summary statistics")
  #PCA plots

  new_df<- imputed_df.3
  new_df<- new_df[,c(-1)]
  cols.num<- colnames(new_df)
  new_df[cols.num]<- sapply(new_df[cols.num], as.numeric)
  df_pca<- prcomp(new_df)


  bc.pca.var<- df_pca$sdev^2
  bc.pca.var<- round(bc.pca.var/sum(bc.pca.var)*100,1)
  df.pca_data<- data.frame(sample= rownames(df_pca$x),
                           x=df_pca$x[,1],
                           y = df_pca$x[,2])

  df.pca_data$type<- sub("_.*", "", df.pca_data$sample)
  df.pca_data$type<- sub("^ZZ.","", df.pca_data$type)
  df.pca_data$type<- sub("^Z.", "", df.pca_data$type)
  df.pca_data<-df.pca_data %>%
    group_by(type) %>%
    mutate(Rep = row_number()) %>%
    ungroup()

  df.pca_data$Rep<- as.character(as.numeric(df.pca_data$Rep))

  a<-ggplot(data=df.pca_data,aes(label=sample,x=x,y=y, color = type, shape = Rep)) +
    ggtitle("Kinome PCA Plot") +
    geom_point(aes( size = 4)) +
    xlab(paste("PC1: ",bc.pca.var[1],"%",sep="")) +
    ylab(paste("PC2: ",bc.pca.var[2],"%",sep=""))+
    scale_color_jco()+
    theme_classic()+
    theme(legend.text = element_text(size = 10))+
    guides(color=guide_legend(override.aes = list(size = 5)),
           shape = guide_legend(override.aes = list(size=5)))


  ggsave(filename = file.path("PCA","PCA plot.svg"), a, device= "svg", width = 10, height = 6)
  ggsave(filename = file.path("PCA","PCA plot.pdf"), a, device= "pdf", width = 10, height = 6)

  log_message("Made PCA plot")
  #Generating a heatmap and clustering using Euclidian distance

  new_df$type<- rownames(new_df)
  new_df$type<- sub("_.*", "", new_df$type)
  new_df$type<- sub("^ZZ.","", new_df$type)
  new_df$type<- sub("^Z.", "", new_df$type)

  new_df<-new_df %>%
    group_by(type) %>%
    mutate(Rep = row_number()) %>%
    ungroup()

  new_df$Rep<- as.character(as.numeric(new_df$Rep))
  new_df$type<- paste0(new_df$type,"-Rep",new_df$Rep)
  new_df2<- new_df[,-((ncol(new_df)-1):ncol(new_df))]
  rownames(new_df2)<- new_df$type
  ## Z-score of the replicates

  df.matrix<- as.matrix(new_df2)
  df.matrix<- as.data.frame(t(df.matrix))

  cal_z_score <- function(f){
    (f - mean(f)) / sd(f)
  }

  df.matrix.z.score <- t(apply(df.matrix, 1, cal_z_score))


  svg(filename = file.path("Heatmaps","Z-score of the Log2 LFQ kinome intensities by replicate.svg"), width = 10, height = 10)

  pheatmap(df.matrix.z.score,
           cluster_rows = T,
           cluster_cols = T,
           clustering_distance_rows = 'euclidean',
           clustering_distance_cols = "euclidean",
           fontsize_row = 3,
           cellwidth = 20,
           colorRampPalette(c("#000080", "white", "#DC143C"))(100),
           angle_col = 45,
           main="Z-score of the Log2 of the Kinase LFQ intensity (By Replicate)")

  dev.off()
  log_message("Made heatmap of each replicate")
  ## Z-score of the average
  new_df3<-new_df[,-ncol(new_df)]
  new_df3$type<- sub("-.*", "", new_df3$type)
  new_df3<- new_df3%>%group_by(type)%>%
    summarise(across(everything(), mean))
  new_df3<- as.data.frame(new_df3)
  rownames(new_df3)<- new_df3$type
  new_df4<- new_df3[,-1]
  df.matrix<- as.matrix(new_df4)
  df.matrix<- as.data.frame(t(new_df4))
  df.matrix.z.score2 <- t(apply(df.matrix, 1, cal_z_score))

  svg(filename = file.path("Heatmaps","Z-score of the averaged Log2 LFQ kinome intensities.svg"), width = 10, height = 10)

  pheatmap(df.matrix.z.score2,
           cluster_rows = T,
           cluster_cols = T,
           clustering_distance_rows = 'euclidean',
           clustering_distance_cols = "euclidean",
           fontsize_row = 3,
           cellwidth = 20,
           colorRampPalette(c("#000080", "white", "#DC143C"))(100),
           angle_col = 45,
           main="Z-score of the Log2 of the Averaged Kinase LFQ intensity")

  dev.off()
  log_message("Made heatmap of averaged replicates")

  # Heatmap of average of significant kinases
  colnames(df)[1]<- "Kinases"
  sig.output<- df %>% filter(`p value` <= 0.05)
  y<- as.character(length(unique(sig.output$Kinases)))
  print(paste0("Number of significant kinases: ", y))
  df.matrix.z.score2<- as.data.frame(df.matrix.z.score2)
  df.matrix.z.score2$Kinases<- rownames(df.matrix.z.score2)
  df.matrix.z.score2<- df.matrix.z.score2[df.matrix.z.score2$Kinases %in% sig.output$Kinases,]
  df.matrix.z.score2<- df.matrix.z.score2[,-c(ncol(df.matrix.z.score2))]

  svg(filename = file.path("Heatmaps","Z-score of the averaged Log2 LFQ significant kinases.svg"), width = 10, height = 10)
  pheatmap(df.matrix.z.score2,
           cluster_rows = T,
           cluster_cols = T,
           clustering_distance_rows = 'euclidean',
           clustering_distance_cols = "euclidean",
           fontsize_row = 3.5,
           cellwidth = 20,
           colorRampPalette(c("#000080", "white", "#DC143C"))(100),
           angle_col = 45,
           main="Z-score of the Log2(Kinase LFQ intensity)")

  dev.off()
  log_message("Made heatmap of averaged replicates - signficant")
  log_message("Finished script")

  return(imputed_df.3)

}

