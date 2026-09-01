#' This function cleans up and analyzes the DIA output from DIANN, adding a pseudocount and imputing
#'
#'
#' @param df Unique protein groups matrix from DIANN
#' @param kinases Either the human or mouse kinome spreadsheet
#' @param peptides report.pr.matrix from DIANN
#' @param metadata File with column titled "Sample.ID" with the column names of the abundance values in the protein groups file, and a column titled "Treatment" with the corresponding treatment replicate
#' @param directory Output folder
#' @param unique_df Either "protein" for pg.matrix or "gene" for unique_gene.matrix
#' @keywords DIANN
#' @examples
#' diann.cleanup()
#' @import dplyr fuzzyjoin stringr tibble arsenal tidyverse tidyr data.table sjmisc ggpubr ggsci ggplot2 svglite rstatix pheatmap arrow reshape2 viridis ggpointdensity purrr
#' @export


diann.cleanup<- function(df, kinases, metadata, directory, unique_df, peptide ){
  ifelse(!dir.exists(directory), dir.create(directory), paste0(directory," folder exists already"))
  setwd(directory)
  ifelse(!dir.exists("Heatmaps"), dir.create("Heatmaps"), "Heatmaps folder exists already")
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

  if(unique_df=="genes"){
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

  kinase.peptides1<- df[df$Genes%in% kinases$Gene,]
  kinase.peptides<- protein[protein$Genes%in% human.kinome$Gene,]
  y<- as.character(length(unique(kinase.peptides1$Genes)))
  x<- as.character(length(unique(kinase.peptides$Genes)))
  log_message(paste0("Number of kinases kept: ", x, "/",y))

  #Match Sample ID to treatment
  protein<-protein %>%
    rename_with(~deframe(metadata)[.x], .cols = metadata$Sample.ID) %>%
    dplyr::select(Genes, any_of(metadata$Treatment))

  #starting stats file
  stats<-as.data.frame()
  stats<-protein%>%

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
    u<- ((df.3[,i]-x)+max.med)
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


  ggsave(filename = file.path("QC","Global abundances histogram.svg"), b, device= "svg", width = 10, height = 6)
  ggsave(filename = file.path("QC","Global pre-abundances histogram.pdf"), b, device= "pdf", width = 10, height = 6)


  c<- ggplot(histogram%>%filter(type=="Pre-median")%>%filter(value!=0), aes(x=(value),y = name))+
    geom_jitter(color="#6287AF", alpha=0.9)+
    ggtitle("Log2(Relative Protein Abundances Pre-Median Normalization)")+
    xlab("log2(Protein Abundances)")+
    theme_bw()+
    stat_summary(fun.x = median, fun.xmin = median, fun.xmax = median,
                 geom = "crossbar", width = 0.5)


  ggsave(filename = file.path("QC","Global abundances pre median normalization.svg"), c, device= "svg", width = 10, height = 6)
  ggsave(filename = file.path("QC","Global abundances pre median normalization.pdf"),c, device= "pdf", width = 10, height = 6)

  d<- ggplot(histogram%>%filter(type=="Post-median")%>%filter(value!=0), aes(x=(value),y = name))+
    geom_jitter(color="#6287AF", alpha=0.9)+
    ggtitle("Log2(Relative Protein Abundances Post-Median Normalization)")+
    xlab("log2(Protein Abundances)")+
    theme_bw()+
    stat_summary(fun.x = median, fun.xmin = median, fun.xmax = median,
                 geom = "crossbar", width = 0.5)


  ggsave(filename = file.path("QC","Global abundances post median normalization.svg"), d, device= "svg", width = 10, height = 6)
  ggsave(filename = file.path("QC","Global abundances post median normalization.pdf"),d, device= "pdf", width = 10, height = 6)

  # Pre-imputation PCA plot
  new_df<- df.4
  new_df<- as.data.frame(t(df.4))
  colnames(new_df)<- new_df[1,]; new_df<- new_df[-c(1),]
  cols.num<- colnames(new_df)
  new_df[cols.num]<- sapply(new_df[cols.num], as.numeric)
  keep_cols <- colMeans(is.na(new_df)) < 0.5
  df_filt <- new_df[, keep_cols]

  # Step 2 — impute remaining NAs with column means (PCA-only)
  for (i in seq_len(ncol(df_filt))) {
    col_vals <- df_filt[, i]
    df_filt[is.na(col_vals), i] <- mean(col_vals, na.rm = TRUE)
  }

  df_pca<- prcomp(df_filt, scale. = T, center = T)


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
    ggtitle("Global PCA Plot pre filtering and imputation") +
    geom_point(aes( size = 4)) +
    xlab(paste("PC1: ",bc.pca.var[1],"%",sep="")) +
    ylab(paste("PC2: ",bc.pca.var[2],"%",sep=""))+
    scale_color_jco()+
    theme_classic()+
    theme(legend.text = element_text(size = 10))+
    guides(color=guide_legend(override.aes = list(size = 5)),
           shape = guide_legend(override.aes = list(size=5)))

  ggsave(filename = file.path("QC","PCA plot pre filter and imputation.svg"), a, device= "svg", width = 10, height = 6)
  ggsave(filename = file.path("QC","PCA plot pre filter and imputation.pdf"), a, device= "pdf", width = 10, height = 6)

  #Filter for columns where at least one group n-1 observed
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

      if(any(df.6$total_non_na >(g-2))){
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
    }
  }


  log_message("Filtered for proteins with enouch sample representation")
  print("Filtered for proteins with enouch sample representation")

  #Multi-model imputation
  df.4<- df.4[colnames(df.4)%in% x]
  df.5<-data.frame(sapply(df.4, function(x) as.numeric(as.character(x))))
  rownames(df.5)<- rownames(df.4)
  print(paste0("The minimum value before adding pseudocount:", min(df.5, na.rm = T)))
  log_message(paste0("The minimum value before adding pseudocount:", min(df.5, na.rm = T)))


  hybrid_impute <- function(data,
                            left_shift = 1.8,
                            width = 0.3,
                            seed = 123) {

    set.seed(seed)
    #floor_method <- "global_min"
    mat <- as.matrix(data)
    #if (ncol(mat) != length(cond)) stop("cond must have same length as ncol(data)")
    cond <- as.factor(sub("_.*","",metadata$Treatment))
    protein_ids <- colnames(mat)
    sample_ids <- rownames(mat)

    #Get min values for both the gloabl df, each sample type, and each protein type
    min_value<- min(data, na.rm = T)
    sample_min <- apply(mat, 1, function(x) if(all(is.na(x))) NA else min(x, na.rm=TRUE))
    protein_min <- apply(mat, 2, function(x) if(all(is.na(x))) NA else min(x, na.rm=TRUE))

    # Outputs
    imputed_mat <- mat
    impute_method <- matrix(NA_character_, nrow=nrow(mat), ncol=ncol(mat),
                            dimnames=list(sample_ids, protein_ids)) #make matrix with same col and row names, same dimensions

    # For each condition compute distribution parameters (aka mean and sd) from observed values
    cond_stats <- lapply(levels(cond), function(clev) {
      rows <- which(cond == clev)
      vals <- mat[rows, ,drop=FALSE]
      obs_vals <- vals[!is.na(vals)]
      list(mean = mean(obs_vals, na.rm=TRUE),
           sd   = sd(obs_vals, na.rm=TRUE),
           row = rows)
    })
    names(cond_stats) <- levels(cond)

    # Loop: for each protein (col) and each condition, impute missing according to rule
    for (i in seq_len(ncol(mat))) {
      for (rname in levels(cond)) {
        rows <- cond_stats[[rname]]$row
        vals <- mat[rows,i]

        n_rep <- length(rows)
        n_obs <- sum(!is.na(vals))
        missing_idx <- which(is.na(vals)) #check for missing values
        if (length(missing_idx) == 0) next  # nothing to impute for this protein/condition

        # observed in >= either 50% or at least 2/3 replicates
        if (n_obs >= max(2,ceiling(n_rep*0.5))) {
          # left-shifted Gaussian using condition-level stats
          cond_mean <- cond_stats[[rname]]$mean
          cond_sd   <- cond_stats[[rname]]$sd
          if (is.na(cond_mean) || is.na(cond_sd) || cond_sd == 0) {
            #use global min for safety
            imputed_vals <- rep(min_value, length(missing_idx))
            method_label <- "global_min"
          } else {
            mu <- cond_mean - left_shift * cond_sd
            sigma <- cond_sd * width
            imputed_vals <- rnorm(length(missing_idx), mean = mu, sd = sigma)
            method_label <- "leftshift_gaussian"
          }
          # write back to the correct columns
          imputed_mat[rows[missing_idx],i] <- imputed_vals
          impute_method[rows[missing_idx],i] <- method_label

        } else{
          # Observed <50% of replicates or is 0: Global LOD
          for (mi in seq_along(missing_idx)) {
            row_idx <- rows[missing_idx[mi]]
            imputed_mat[row_idx,i] <- min_value
            impute_method[row_idx,i] <- "global_min"
          }
        }
      } # end per-condition
    } # end per-protein

    return(list(imputed = imputed_mat, method = impute_method))
  }

  res<- hybrid_impute(df.5)
  df_imputed <- as.data.frame(res$imputed)
  df_method  <- as.data.frame(res$method)

  # imputation checks
  #Looking at number of missing values
  missing_mat <- is.na(df.5) * 1
  svg(filename = file.path("QC","Heatmap of missing values.svg"), width = 10, height = 10)
  pheatmap(missing_mat,
           cluster_rows=TRUE,
           cluster_cols=TRUE,
           border_color = "black",
           show_colnames = FALSE,)

  dev.off()
  # Looking at distribution of intensities of imputated values compared to non-missing
  test.2<- df.5
  test.2$Sample<- rownames(df.5)
  orig_df <- reshape2::melt(test.2)
  colnames(orig_df) <- c("Sample","Protein","Intensity")
  orig_df$Type <- ifelse(is.na(orig_df$Intensity),"Missing","Observed")
  df_imputed.2<-df_imputed
  df_imputed.2$Sample<- rownames(df_imputed)
  imp_df <- reshape2::melt(df_imputed.2); colnames(imp_df) <- c("Sample","Protein","Intensity")
  imp_df$Type <- "Imputed"
  # only keep imputed cells where original was missing
  imputed_only <- merge(imp_df, orig_df[orig_df$Type=="Missing", c("Sample","Protein")], by=c("Sample","Protein"))

  impute<-ggplot2::ggplot() +
    geom_density(data = orig_df[orig_df$Type=="Observed",], aes(Intensity), alpha=0.3) +
    geom_density(data = imputed_only, aes(Intensity), color="red") +
    facet_wrap(~sub("_.*","", Sample), scales="free") +
    ggtitle("Observed (grey) and Imputed (red) densities by Condition")+
    theme_bw()
  ggsave(filename = file.path("QC","Imputation_check.svg"), impute, device= "svg", width = 10, height = 6)
  ggsave(filename = file.path("QC","Imputation_check.svg.pdf"), impute, device= "pdf", width = 10, height = 6)

  #Export csv of where and what kind of imputaiton was done
  write.csv(df_method, file = file.path("QC","Imputation_method.csv"))

  ##Done checking


  histo3<- as.data.frame(t(df_imputed)); histo3$Genes<- rownames(histo3)
  histo3<- histo3%>%
    pivot_longer(!Genes)%>%
    mutate(value= tidyr::replace_na(value, 0))%>%
    mutate(type ="Post-imputation")%>%
    mutate(name = sub("Z.","",name))%>%
    mutate(name = sub("ZZ.","",name))

  histogram<- rbind(histogram, histo3)
  histogram$value<- as.numeric(as.character(histogram$value))
  e<- ggplot(histogram, aes(x=(value),fill = type))+
    geom_histogram(color="#e9ecef", alpha=0.6, position= "identity")+
    scale_fill_manual(values=c("red","#69b3a2", "#404080")) +
    ggtitle("Log2(Relative Protein Abundances) Post Filtering and Imputation")+
    xlab("log2(Protein Abundances)")+
    theme_bw()+
    facet_wrap(~name, ncol = 3)

  e
  ggsave(filename = file.path("QC","Global abundances post imputation and filtering.svg"), e, device= "svg", width = 10, height = 6)
  ggsave(filename = file.path("QC","Global abundances post imputation and filtering.pdf"),e, device= "pdf", width = 10, height = 6)

  imputed_df.2<- df_imputed
  imputed_df.2$Treatment<- rownames(imputed_df.2)
  imputed_df.2$Treatment<- sub("-.*", '', imputed_df.2$Treatment)
  imputed_df.2$Treatment<- sub("_.*", '', imputed_df.2$Treatment)
  imputed_df.2<- imputed_df.2[,c(ncol(imputed_df.2), 1:ncol(imputed_df.2)-1)]



  log_message("Imputed")

  density<- imputed_df.3%>%
    tidyverse::rownames_to_columns("Treatment")%>%pivot_longer(cols = !Treatment)
  p2 <- ggplot(density, aes(x=value, group=Treatment, fill=Treatment)) +
    geom_density(adjust=1.5, alpha=.2) +
    theme_bw()
  ggsave(filename = file.path("QC","Final Global abundances density.svg"), p2, device= "svg", width = 10, height = 6)
  ggsave(filename = file.path("QC","Final Global abundances density.pdf"),p2, device= "pdf", width = 10, height = 6)



  new_df<- imputed_df.2
  new_df<- new_df[,c(-1)]
  cols.num<- colnames(new_df)
  new_df[cols.num]<- sapply(new_df[cols.num], as.numeric)
  df_pca<- prcomp(new_df, scale. = T, center = T)


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

  ggsave(filename = file.path("PCA","All proteins PCA plot scaled.svg"), a, device= "svg", width = 10, height = 6)
  ggsave(filename = file.path("PCA","All proteins plot.pdf"), a, device= "pdf", width = 10, height = 6)

  cleaned_up_global<- df_imputed
  cleaned_up_global$Treatment<- sub("Z.","",rownames(cleaned_up_global))
  cleaned_up_global$Treatment<- sub("ZZ.","",cleaned_up_global$Treatment)

  write.csv(cleaned_up_global, file = file.path("Results","Relative Global Protein Abundance.csv"))

  #QC stats



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

  density<- imputed_df.3%>%pivot_longer(cols = !Treatment)
  p2 <- ggplot(density, aes(x=value, group=Treatment, fill=Treatment)) +
    geom_density(adjust=1.5, alpha=.2) +
    theme_bw()
  ggsave(filename = file.path("QC","Final Kinase abundances density.svg"), p2, device= "svg", width = 10, height = 6)
  ggsave(filename = file.path("QC","Final Kinase abundances density.pdf"),p2, device= "pdf", width = 10, height = 6)



  cleaned_up_kinases<- imputed_df.3
  cleaned_up_kinases$Treatment<- sub("Z.","",rownames(cleaned_up_kinases))
  cleaned_up_kinases$Treatment<- sub("ZZ.","",cleaned_up_kinases$Treatment)
  write.csv(cleaned_up_kinases, file = file.path("Results","Relative Kinase Protein Abundance.csv"))

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
  df_pca<- prcomp(new_df, scale=T, center = T)


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


  heatmap_variables<-function(mat){
    n_rows <- nrow(mat)
    n_cols <- ncol(mat)
    if (n_rows <= 100) {
      cell_h <- 12
      font_r <- 10
    } else if (n_rows <= 250) {
      cell_h <- 8
      font_r <- 7
    } else if (n_rows <= 450) {   # covers your max
      cell_h <- 5
      font_r <- 4
    } else {
      cell_h <- 3
      font_r <- 2
    }
    if (n_cols <= 8) {
      cell_w <- 25
    } else if (n_cols <= 20) {
      cell_w <- 20
    } else if (n_cols <= 40) {
      cell_w <- 15
    } else {
      cell_w <- 10
    }

    show_rows <- n_rows <= 450
    svg_height_in <- max(6, (n_rows * cell_h) / 72 + 2)
    svg_width_in <- max(8, min(40, (n_cols * cell_w) / 72 + 2))

    list(
      n_rows = n_rows,
      n_cols = n_cols,
      cell_h = cell_h,
      font_r = font_r,
      cell_w = cell_w,
      show_rows = show_rows,
      svg_height_in = svg_height_in,
      svg_width_in = svg_width_in
    )
  }
  vars<- heatmap_variables(df.matrix.z.score)


  svg(filename = file.path("Heatmaps","Z-score of the Log2 LFQ kinome intensities by replicate.svg"),
      width =vars$svg_width_in, height = vars$svg_height_in)

  pheatmap(df.matrix.z.score,
           cluster_rows = T,
           cluster_cols = T,
           clustering_distance_rows = 'euclidean',
           clustering_distance_cols = "euclidean",
           fontsize_row = vars$font_r,
           cellheight = vars$cell_h,
           cellwidth = vars$cell_w,
           show_rownames = vars$show_rows,
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

  vars<- heatmap_variables(df.matrix.z.score2)

  svg(filename = file.path("Heatmaps","Z-score of the averaged Log2 LFQ kinome intensities.svg"),
      width = vars$svg_width_in, height = vars$svg_height_in)

  pheatmap(df.matrix.z.score2,
           cluster_rows = T,
           cluster_cols = T,
           clustering_distance_rows = 'euclidean',
           clustering_distance_cols = "euclidean",
           fontsize_row = vars$font_r,
           cellheight = vars$cell_h,
           cellwidth = vars$cell_w,
           show_rownames = vars$show_rows,
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

  ##Heatmap of each replicate
  df.matrix.z.score3<- as.data.frame(df.matrix.z.score)
  df.matrix.z.score3$Kinases<- rownames(df.matrix.z.score3)
  df.matrix.z.score3<- df.matrix.z.score3[df.matrix.z.score3$Kinases %in% sig.output$Kinases,]
  df.matrix.z.score3<- df.matrix.z.score3[,-c(ncol(df.matrix.z.score3))]

  vars<- heatmap_variables(df.matrix.z.score3)
  svg(filename = file.path("Heatmaps","Z-score of the Log2 LFQ significant kinases reps.svg"),
      width = vars$svg_width_in, height = vars$svg_height_in)
  pheatmap(df.matrix.z.score3,
           cluster_rows = T,
           cluster_cols = T,
           clustering_distance_rows = 'euclidean',
           clustering_distance_cols = "euclidean",
           fontsize_row = vars$font_r,
           cellheight = vars$cell_h,
           cellwidth = vars$cell_w,
           show_rownames = vars$show_rows,
           colorRampPalette(c("#000080", "white", "#DC143C"))(100),
           angle_col = 45,
           main="Z-score of the Significant Log2(Kinase LFQ intensity)- By Rep")

  dev.off()
  log_message("Made heatmap of replicates - signficant")
  #Heatmap of average
  df.matrix.z.score2<- as.data.frame(df.matrix.z.score2)
  df.matrix.z.score2$Kinases<- rownames(df.matrix.z.score2)
  df.matrix.z.score2<- df.matrix.z.score2[df.matrix.z.score2$Kinases %in% sig.output$Kinases,]
  df.matrix.z.score2<- df.matrix.z.score2[,-c(ncol(df.matrix.z.score2))]

  vars<- heatmap_variables(df.matrix.z.score2)

  svg(filename = file.path("Heatmaps","Z-score of the averaged Log2 LFQ significant kinases.svg"),
      width = vars$svg_width_in, height = vars$svg_height_in)
  pheatmap(df.matrix.z.score2,
           cluster_rows = T,
           cluster_cols = T,
           clustering_distance_rows = 'euclidean',
           clustering_distance_cols = "euclidean",
           fontsize_row = vars$font_r,
           cellheight = vars$cell_h,
           cellwidth = vars$cell_w,
           show_rownames = vars$show_rows,
           colorRampPalette(c("#000080", "white", "#DC143C"))(100),
           angle_col = 45,
           main="Z-score of the ANOVA Significant \nLog2(Kinase LFQ intensity)")

  dev.off()
  log_message("Made heatmap of averaged replicates - signficant")
  log_message("Finished script")


  sig.output<- df %>% filter(`p value` <= 0.05)%>% arrange(`p value`) %>% head(100)
  ##Heatmap of each replicate top 100
  df.matrix.z.score4<- df.matrix.z.score3[rownames(df.matrix.z.score3) %in% sig.output$Kinases,]
  #df.matrix.z.score4<- df.matrix.z.score4[,-c(ncol(df.matrix.z.score4))]

  ##Heatmap of each replicate

  print(head(df.matrix.z.score4))

  vars<- heatmap_variables(df.matrix.z.score4)
  svg(filename = file.path("Heatmaps","Top 100 significant proteins-rep.svg"),
      width = vars$svg_width_in, height = vars$svg_height_in)
  pheatmap(df.matrix.z.score4,
           cluster_rows = T,
           cluster_cols = T,
           clustering_distance_rows = 'euclidean',
           clustering_distance_cols = "euclidean",
           fontsize_row = vars$font_r,
           cellheight = vars$cell_h,
           cellwidth = vars$cell_w,
           show_rownames = vars$show_rows,
           colorRampPalette(c("#000080", "white", "#DC143C"))(100),
           angle_col = 45,
           main="Z-score of top 100 Significant \nLog2(Kinase LFQ intensity)- By Rep")

  dev.off()
  log_message("Made heatmap of replicates - signficant")

  ##Heatmap of each average top 100
  df.matrix.z.score5<- df.matrix.z.score2[rownames(df.matrix.z.score2) %in% sig.output$Kinases,]
  #df.matrix.z.score5<- df.matrix.z.score5[,-c(ncol(df.matrix.z.score5))]

  ##Heatmap of each replicate

  print(head(df.matrix.z.score5))

  vars<- heatmap_variables(df.matrix.z.score5)
  svg(filename = file.path("Heatmaps","Top 100 significant proteins-average.svg"),
      width = vars$svg_width_in, height = vars$svg_height_in)
  pheatmap(df.matrix.z.score5,
           cluster_rows = T,
           cluster_cols = T,
           clustering_distance_rows = 'euclidean',
           clustering_distance_cols = "euclidean",
           fontsize_row = vars$font_r,
           cellheight = vars$cell_h,
           cellwidth = vars$cell_w,
           show_rownames = vars$show_rows,
           colorRampPalette(c("#000080", "white", "#DC143C"))(100),
           angle_col = 45,
           main="Z-score of top 100 Significant \nLog2(Kinase LFQ intensity)")

  dev.off()
  log_message("Made heatmap of average - signficant")

  return(imputed_df.3)

}
