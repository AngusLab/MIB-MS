#' This function cleans up and analyzes the DDA output from MaxQuant LFQ
#'
#'
#' @param df proteinGroups.csv from MaxQuant
#' @param kinases Either the human or mouse kinome spreadsheet
#' @param metadata File with column titled "Sample.ID" with the column names of the abundance values in the protein groups file, and a column titled "Treatment" with the corresponding treatment replicate
#' @keywords DIANNE
#' @examples
#' maxquant.cleanup()
#' @import dplyr fuzzyjoin stringr tibble arsenal tidyverse tidyr data.table sjmisc ggpubr ggsci ggplot2 svglite rstatix pheatmap
#' @export

maxquant.cleanup<- function(df, kinases, metadata){

  #grabbing only the necessary columns
  df.2<- df[,c(7,11)]
  df.3<- select(df, contains("LFQ"))
  df.2<- cbind(df.2, df.3)

  #filtering out <2 Razor + unique peptides and rows where the gene name is blank
  df.2<- df.2%>%filter(`Razor + unique peptides` > 1)
  df.2<- df.2[!df.2$`Gene names` == "",]
  print("Filtered for >1 unique peptide")
  #Matching to either the human or mouse kinome
  df.3<- kinases %>% inner_join(df.2, by=c("Gene"= "Gene names"))
  df.3<- df.3[,-c(3)]
  print("Matched to kinome")
  #Keeping only the sample names
  x<- df.3[,-c(1:2)]
  #y<- colnames(x)
  #a<-sub('.*\\_',"",y)

  #colnames(x)<- a
  #Log2 transform
  x<- log2(x)
  x[sapply(x, is.infinite)] <- NA
  df.3<- df.3[,c(1:2)]
  df.3<- cbind(df.3, x)
  write.csv(df.3, "Log2LFQforPerseus.txt")
  print("Log transformed")
  #Match Sample ID to treatment
  df.4<-df.3 %>%
    rename_with(~deframe(metadata)[.x], .cols = metadata$Sample.ID)%>%  select(Gene, Group, any_of(metadata$Treatment))
  print("Matched to sample ID")

  df.4<- as.data.frame(t(df.4))
  colnames(df.4)<- df.4[1,]
  df.4<- df.4[-c(1,2),]
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


  print("Filtered for Kinases with enough sample representation")

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
  imputed_df.2<- inner_join(metadata, imputed_df.2, by = "Treatment")

  imputed_df.2$Treatment<- sub("-.*", '', imputed_df.2$Treatment)
  imputed_df.2$Treatment<- sub("_.*", '', imputed_df.2$Treatment)
  imputed_df.2<- imputed_df.2[,-c(1)]
  imputed.df3<- imputed_df.2
  imputed.df3$Treatment<- sub("^ZZ.","", imputed.df3$Treatment)
  imputed.df3$Treatment<- sub("^Z.", "", imputed.df3$Treatment)


  print("Imputed")

  write.table(imputed.df3, "Kinase_protein_abundances.txt")
  write.csv(imputed.df3, "Kinase_protein_abundances.csv")
  write.csv(imputed_df.2, "Start.here.for.stats.function.csv")

  DF<- imputed.df3
  kinases<-colnames(DF)
  #choose which columns to drop, typically this is the first 2 columns but should check each time
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
  write.csv(df, file="stats_allgroups.csv")

  print("Completed Summary statistics")
  #PCA plots
  new_df<- imputed.df3
  new_df$Treatment<- paste0(new_df$Treatment, "-", rownames(new_df))

  rownames(new_df)<- new_df[,1]
  new_df<- new_df[,-c(1)]

  cols.num<- colnames(new_df)
  new_df[cols.num]<- sapply(new_df[cols.num], as.numeric)
  df_pca<- prcomp(new_df)


  bc.pca.var<- df_pca$sdev^2
  bc.pca.var<- round(bc.pca.var/sum(bc.pca.var)*100,1)
  df.pca_data<- data.frame(sample= rownames(df_pca$x),
                           x=df_pca$x[,1],
                           y = df_pca$x[,2])

  df.pca_data$type<- sub("-.*", "", df.pca_data$sample)
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

  a
  ggsave("PCA plot.svg", a, device= "svg", width = 10, height = 6)
  print("Made PCA plot")
  #Generating a heatmap and clustering using Euclidian distance

  new_df$type<- rownames(new_df)
  new_df$type<- sub("-.*", "", new_df$type)
  new_df<-new_df %>%
    group_by(type) %>%
    mutate(Rep = row_number()) %>%
    ungroup()

  new_df$Rep<- as.character(as.numeric(new_df$Rep))
  new_df$type<- paste0(new_df$type,"-Rep",new_df$Rep)
  rownames(new_df)<- new_df$type
  new_df2<- new_df[,-((ncol(new_df)-1):ncol(new_df))]
  rownames(new_df2)<- rownames(new_df)

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
  colnames(df)[ncol(df)]<- "p.value"
  colnames(df)[1]<- "Kinases"
  sig.output<- df %>% filter(p.value <= 0.05)
  df.matrix.z.score3<- as.data.frame(df.matrix.z.score2)
  df.matrix.z.score3$Kinases<- rownames(df.matrix.z.score3)
  df.matrix.z.score3<- df.matrix.z.score3[df.matrix.z.score3$Kinases %in% sig.output$Kinases,]
  df.matrix.z.score3<- df.matrix.z.score3[,-c(ncol(df.matrix.z.score3))]

  svg("Heatmap of the z-score of the averaged Log2 LFQ kinome intensities significant kinases.svg", width = 10, height = 10)

  pheatmap(df.matrix.z.score3,
           cluster_rows = T,
           cluster_cols = T,
           clustering_distance_rows = 'euclidean',
           clustering_distance_cols = "euclidean",
           fontsize_row = 8,
           cellwidth = 20,
           colorRampPalette(c("#000080", "white", "#DC143C"))(100),
           angle_col = 45,
           main="Z-score of the Log2(Kinase LFQ intensity)")

  dev.off()
  print("Made heatmap of averaged replicates - signficant")
  return(imputed_df.2)
}

