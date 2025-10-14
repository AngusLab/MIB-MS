#' This function cleans up and analyzes the DIA output from DIANN
#'
#'
#' @param df Unique protein groups matrix from DIANN
#' @param kinases Either the human or mouse kinome spreadsheet
#' @param peptides report.parquet from DIANN
#' @keywords DIANN
#' @export
#' @examples
#' diann.cleanup()

diann.cleanup<- function(df, kinases ,peptides){
  
  #Matching to either the human or mouse kinome
  df.3<- kinases %>% inner_join(df, by=c("Gene"= "Genes"))
  df.3<- df.3[,-c(2)]
  print(("Matched to kinome"))
  
  #filtering out <2 Razor + unique peptides and rows where the gene name is blank
  
  dia.kinase.names<- as.vector(df.3$Gene)
  
  peptides.kinases<- peptide[peptide$Genes %in% dia.kinase.names,]
  ## Get rid of duplicate precursor IDs and getting rid of peptides that aren't proteotypic
  peptides.kinases.2<- peptides.kinases[!duplicated(peptides.kinases$Precursor.Id),]
  peptides.kinases.2<- peptides.kinases.2[,c(1,4,10,15,30)]
  peptides.kinases.2<- peptides.kinases.2%>%filter(Proteotypic != 0)
  peptides.kinases.3<- as.data.frame(table(peptides.kinases.2$Genes))
  kinase.1.peptide<- peptides.kinases.3%>% filter(Freq==1)
  y<- as.character(length(unique(kinase.1.peptide$Var1)))
  print(paste0("Number of kinases with only 1 peptide: ", y))
  df.3<- df.3[!df.3$Gene %in% kinase.1.peptide$Var1,]
  #Get rid of those kinases
  
  #Log2 transform
  x<- df.3[,-c(1)]
  x<- log2(x)
  x[sapply(x, is.infinite)] <- NA
  df.3<- df.3[,c(1)]
  df.3<- cbind(df.3, x)
  colnames(df.3)[1]<- "Genes"
  print("Log transformed")
  #Match Sample ID to treatment 
  df.4<-df.3 %>% 
    rename_with(~deframe(sample)[.x], .cols = sample$Sample.ID) %>% 
    dplyr::select(Genes, any_of(sample$Treatment))
  print("Matched to sample ID")
  
  df.4<- as.data.frame(t(df.4))
  colnames(df.4)<- df.4[1,]
  df.4<- df.4[-c(1),]
  df.4$type<- rownames(df.4)
  df.4<- df.4[,c(ncol(df.4), 1:ncol(df.4)-1)]
  df.4$type<- sub("-.*", '', df.4$type)
  df.4$type<- sub("_.*", '', df.4$type)
  
  x<- as.vector(NULL)
  g<- as.data.frame(table(df.4$type))
  g<-max(g$Freq)
  
  
  if(g>2){#Filter for columns where at least one group has 3 values
    for(i in 2:ncol(df.4)){
      a<- df.4[(1)]
      b<- df.4[i]
      df.5<- cbind(a,b)
      y<- colnames(df.5)[2]
      colnames(df.5)[2]<- "Kinase"
      df.6<- df.5 %>%
        group_by(type)%>%
        summarise(total_non_na = sum(!is.na(Kinase)))
      
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
      colnames(df.5)[2]<- "Kinase"
      df.6<- df.5 %>%
        group_by(type)%>%
        summarise(total_non_na = sum(!is.na(Kinase)))
      
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
      colnames(df.5)[2]<- "Kinase"
      df.6<- df.5 %>%
        group_by(type)%>%
        summarise(total_non_na = sum(!is.na(Kinase)))
      
      if(any(df.6$total_non_na==1)){
        x<- append(x, y)
        
      }
      
      
      
    }}
  
  
  print("Filtered for Kinases with enouch sample representation")
  
  df.4<- df.4[colnames(df.4)%in% x]
  df.5<-data.frame(sapply(df.4, function(x) as.numeric(as.character(x))))
  rownames(df.5)<- rownames(df.4)
  
  #Impute 
  impute_normal <- function(object, width=0.3, downshift=1.8, seed=100) {
    
    if (!is.matrix(object)) object <- as.matrix(object)
    mx <- max(object, na.rm=TRUE)
    mn <- min(object, na.rm=TRUE)
    if (mx - mn > 20) warning("Please make sure the values are log-transformed")
    
    set.seed(seed)
    object <- apply(object, 2, function(temp) {
      temp[!is.finite(temp)] <- NA
      temp_sd <- stats::sd(temp, na.rm=TRUE)
      temp_mean <- mean(temp, na.rm=TRUE)
      shrinked_sd <- width * temp_sd   # shrink sd width
      downshifted_mean <- temp_mean - downshift * temp_sd   # shift mean of imputed values
      n_missing <- sum(is.na(temp))
      temp[is.na(temp)] <- stats::rnorm(n_missing, mean=downshifted_mean, sd=shrinked_sd)
      temp
    })
    return(object)
  }
  
  
  imputed_df<- impute_normal(df.5)
  imputed_df.2<- as.data.frame(imputed_df)
  
  imputed_df.2$Treatment<- rownames(imputed_df.2)
  imputed_df.2$Treatment<- sub("-.*", '', imputed_df.2$Treatment)
  imputed_df.2$Treatment<- sub("_.*", '', imputed_df.2$Treatment)
  imputed_df.2<- imputed_df.2[,c(ncol(imputed_df.2), 1:ncol(imputed_df.2)-1)]
  
  print("Imputed")
  
  write.table(imputed_df.2, "Relative Kinase Protein Abundance.txt")
  write.csv(imputed_df.2, "Relative Kinase Protein Abundance.csv")
  
  DF<- imputed_df.2
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
                      data=DF), text=NULL)
  df<-as.data.frame(df)
  write.csv(df, file="ANOVA.csv")
  
  print("Completed Summary statistics")
  #PCA plots
  
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
    ggtitle("Kinome PCA Plot") +
    geom_point(aes( size = 4)) +
    xlab(paste("PC1: ",bc.pca.var[1],"%",sep="")) +
    ylab(paste("PC2: ",bc.pca.var[2],"%",sep=""))+
    scale_color_jco()+
    theme_classic()+
    theme(legend.text = element_text(size = 10))+
    guides(color=guide_legend(override.aes = list(size = 5)),
           shape = guide_legend(override.aes = list(size=5)))
  
  
  ggsave("PCA plot.svg", a, device= "svg", width = 10, height = 6)
  ggsave("PCA plot.pdf", a, device= "pdf", width = 10, height = 6)
  
  print("Made PCA plot")
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
  
  
  svg("Heatmap of the z-score of the Log2 LFQ kinome intensities by replicate.svg", width = 10, height = 10)
  
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
  print("Made heatmap of each replicate")
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
  
  svg("Heatmap of the z-score of the averaged Log2 LFQ kinome intensities.svg", width = 10, height = 10)
  
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
  print("Made heatmap of averaged replicates")
  
  # Heatmap of average of significant kinases
  colnames(df)[1]<- "Kinases"
  sig.output<- df %>% filter(`p value` <= 0.05)
  y<- as.character(length(unique(sig.output$Kinases)))
  print(paste0("Number of significant kinases: ", y))
  df.matrix.z.score2<- as.data.frame(df.matrix.z.score2)
  df.matrix.z.score2$Kinases<- rownames(df.matrix.z.score2)
  df.matrix.z.score2<- df.matrix.z.score2[df.matrix.z.score2$Kinases %in% sig.output$Kinases,]
  df.matrix.z.score2<- df.matrix.z.score2[,-c(ncol(df.matrix.z.score2))]
  
  svg("Heatmap of the z-score of the averaged Log2 LFQ significant kinases.svg", width = 10, height = 10)
  
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
  print("Made heatmap of averaged replicates - signficant")
  
  
  return(imputed_df.2)
}
