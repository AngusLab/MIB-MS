#' This function conducts t-tests across all parameters
#'
#'
#' @param DF output from either diann.cleanup or maxquant.cleanup
#' @param directory Output directory
#' @keywords stats
#'
#' @examples
#' statistical.testing.test()
#' @import dplyr fuzzyjoin stringr tibble arsenal tidyverse tidyr data.table sjmisc ggpubr ggsci ggplot2 svglite rstatix pheatmap
#' @export

statistical.testing.test<-function(DF, directory){

  ifelse(!dir.exists(directory), dir.create(directory), paste0(directory," folder exists already"))
  setwd(directory)
  ifelse(!dir.exists("Volcano_plots"), dir.create("Volcano_plots"), "Volcano_plots folder exists already")
  ifelse(!dir.exists("Results"), dir.create("Results"), "Results folder exists already")


  log_message <- function(message, file = "Stats_log.txt") {
    timestamp <- format(Sys.time(), "[%Y-%m-%d %H:%M:%S]")
    cat(timestamp, message, "\n", file = file, append = TRUE)
    print(message)
  }
  all_kinase<- as.data.frame(NULL)
  for (i in 3:ncol(kinome.analysis)) {

    kinase_name<-colnames(kinome.analysis)[i]
    test<-DF[,c(1,i)]
    colnames(test)<- c("Treatment", "Kinase")
    treatment<- as.vector(unique(test$Treatment))
    ttest<- as.data.frame((NULL))

    for (j in 1:(length(treatment)-1)){
      for (k in (j+1):(length(treatment))){
        A<-treatment[j]
        B<- treatment[k]
        if(is.na(B)) next

        test_sub<-test[test$Treatment %in% c(A,B),]
        tryCatch({
          result<-test_sub%>%
            rstatix::t_test(Kinase~Treatment, detailed = T)
          result$Kinase<- kinase_name
          result<- result[,c(5,6,10,2,3,1,16)]
          ttest<- rbind(ttest, result)
        }, error=function(e){
          print(paste0(e$message, "at Kinase:", kinase_name, " ", A, " vs. ", B))
        })
      }
    }
    all_kinase<- rbind(all_kinase, ttest)
  }
  all_kinase <- all_kinase %>%
    rstatix::adjust_pvalue(method = "BH")%>%
    as.data.frame()

  ##compute - log10p.value
  all_kinase$neglog10.p.value = -log10(all_kinase$p)
  ##compute -log10BH.p.adj
  all_kinase$neglog10.BH.adj.p.value = -log10(all_kinase$p.adj)


  colnames(all_kinase)<- c("group1", "group2", "p.value",  "group1.mean", "group2.mean", "LFC", "Kinase","BH.adj.p.value", "neglog10.p.value", "neglog10.BH.adj.p.value")

  y<- as.character(length(unique(all_kinase$Kinase)))
  log_message(paste0("Number of unique kinases: ", y))

  ###Set significance threshold for colors (blue, grey, red)
  ####raw p value
  all_kinase$significant <- "ns"
  all_kinase$significant[all_kinase$p.value < 0.05 & all_kinase$LFC <= -0.5] <- "down"
  all_kinase$significant[all_kinase$p.value < 0.05 & all_kinase$LFC >=  0.5] <- "up"
  ####adjusted p value
  all_kinase$significant_adj <- "ns"
  all_kinase$significant_adj[all_kinase$BH.adj.p.value < 0.05 & all_kinase$LFC <= -0.5] <- "down"
  all_kinase$significant_adj[all_kinase$BH.adj.p.value < 0.05 & all_kinase$LFC >=  0.5] <- "up"



  all_kinase$group1 <- sub("^ZZ\\.", "", sub("^Z\\.", "", all_kinase$group1))
  all_kinase$group2 <- sub("^ZZ\\.", "", sub("^Z\\.", "", all_kinase$group2))


  write.csv(all_kinase,file=file.path("Results","T.test comparisons all kinases.csv"))


  all_kinase$ID <- paste0(all_kinase$group1, " v. ", all_kinase$group2)

  log_message("Making pdf volcano plots")

  make_volcano <- function(df, x_var, y_var, y_lab, significant_lab ,out_dir, n_labels=25) {
    A <- unique(df$group1)
    B <- unique(df$group2)

    sig_df   <- df[df[[significant_lab]] != "ns", ]
    top_hits <- sig_df[order(sig_df[[y_var]], decreasing = TRUE), ]
    top_hits <- head(top_hits, n_labels)

    p <- ggplot2::ggplot(df, aes(x = .data[[x_var]], y = .data[[y_var]], color = .data[[significant_lab]])) +
      geom_point(size = 1.5, alpha = 0.7) +
      theme_bw(base_size = 14) +
      geom_hline(yintercept = 1.3,  linetype = "dashed", color = "grey40", linewidth = 0.5) +
      geom_vline(xintercept =  0.5, linetype = "dashed", color = "grey40", linewidth = 0.5) +
      geom_vline(xintercept = -0.5, linetype = "dashed", color = "grey40", linewidth = 0.5) +
      scale_color_manual(values = c(up = "#A73030FF", down = "#003C67FF", ns = "#868686FF")) +
      ggrepel::geom_label_repel(data = top_hits,
                                aes(label = Kinase),
                                size = 3.5,
                                max.overlaps = Inf,
                                box.padding = 0.5,
                                force=10,
                                force_pull = 0.5,
                                ylim=c(0,NA),
                                fill=alpha("white",0.7),
                                label.size = NA,
                                segment.color = "grey50",
                                segment.size = 0.3,
                                color="black",
                                show.legend = FALSE) +
      labs(x = expression(log[2] ~ "Fold Change MIB Binding"),
           y = y_lab,
           title= paste(A, "vs", B)) +
      theme(plot.title= element_text(hjust = 0.5, face = "bold", size = 16),
            axis.text= element_text(color = "black"),
            legend.position= "none",
            panel.grid.minor = element_blank(),
            panel.grid.major = element_line(color = "grey90"),
            plot.margin = margin(10, 10, 10, 10))+
      coord_cartesian(clip="off")+
      expand_limits(y=max(df[[y_var]])*1.2)

    ggsave(p, file = file.path(out_dir, paste0(A, " vs ", B, ".svg")), device = "svg", width = 8, height = 8)
    ggsave(p, file = file.path(out_dir, paste0(A, " vs ", B, ".pdf")), device = "pdf", width = 8, height = 8)
    #log_message(paste0("Saved volcano plot: ", A, " vs. ", B))
  }

  for (group_id in unique(all_kinase$ID)) {
    df <- all_kinase[all_kinase[["ID"]] == group_id, ]
    make_volcano(df, "LFC", "neglog10.p.value","-log10 p-value","significant","Volcano_plots",n_labels = 25)
    make_volcano(df, "LFC", "neglog10.BH.adj.p.value", "-log10 B.H. adj p-value","significant_adj","Adjusted Volcano_plots",n_labels = 25)
  }

  return(all_kinase)
}

