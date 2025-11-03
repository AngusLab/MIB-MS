#' This function conducts t-tests across all parameters
#'
#'
#' @param DF output from either diann.cleanup or maxquant.cleanup
#' @keywords stats
#'
#' @examples
#' statistical.testing()
#' @import dplyr fuzzyjoin stringr tibble arsenal tidyverse tidyr data.table sjmisc ggpubr ggsci ggplot2 svglite rstatix pheatmap
#' @export

statistical.testing<-function(DF){

  all_kinase<- as.data.frame(NULL)
  for (i in 2:ncol(DF)) {
    treatment<- DF[1]
    test<- DF[i]
    x<- colnames(test)
    test<- cbind(treatment, test)
    colnames(test)<- c("Treatment", "Kinase")

    e<-test%>%
      t_test(Kinase~Treatment, detailed = T)%>%
      adjust_pvalue(method = "BH") %>%
      add_significance("p.adj")
    #print("done with t.test")
    e$Kinase<- x
    final<- e[,c(5,6,10,16,2,3,1,18)]


    all_kinase<- rbind(all_kinase, final)
  }

  ##compute - log10p.value
  all_kinase$neglog10.p.value = -log10(all_kinase$p)
  ##compute -log10BH.p.adj
  all_kinase$neglog10.BH.adj.p.value = -log10(all_kinase$p.adj)


  colnames(all_kinase)<- c("group1", "group2", "p.value", "BH.adj.p.value", "group1.mean", "group2.mean", "LFC", "Kinase", "neglog10.p.value", "neglog10.BH.adj.p.value")

  y<- as.character(length(unique(all_kinase$Kinase)))
  print(paste0("Number of unique kinases: ", y))

  ###Set significance threshold for colors (blue, grey, red)
  all_kinase$significant <- "ns"
  all_kinase[all_kinase$p.value<0.05 & all_kinase$LFC <= (-1) , "significant"] <- "down"
  all_kinase[all_kinase$p.value<0.05 & all_kinase$LFC >= (1) , "significant"] <- "up"


  ####Set threshold of what kinases to label
  all_kinase$labels<-"FALSE"
  all_kinase[all_kinase$p.value<0.05 & all_kinase$LFC <= (-1) , "labels"] <- "TRUE" ###Adjust based on density of data
  all_kinase[all_kinase$p.value<0.05 & all_kinase$LFC >= (1) , "labels"] <- "TRUE" ###Adjust based on density of data

  all_kinase$group1 <- sub("^ZZ.", "", all_kinase$group1)
  all_kinase$group2 <- sub("^ZZ.", "", all_kinase$group2)
  all_kinase$group1 <- sub("^Z.", "", all_kinase$group1)
  all_kinase$group2 <- sub("^Z.", "", all_kinase$group2)


  write.csv(all_kinase,"T.test comparisons all kinases.csv")


  all_kinase$ID<- paste0(all_kinase$group1, " v. ", all_kinase$group2)
  groups<- as.vector(unique(all_kinase$ID))
  print("Making pdf volcano plots")
  for (i in groups) {
    df<- all_kinase%>% filter(ID == i)
    A<- unique(df$group1)
    B<- unique(df$group2)
    p<-ggscatter(df, x="LFC", y="neglog10.p.value",
                 title=paste(A,"vs", B),
                 label=ifelse(df$labels, df$Kinase, ""),
                 font.label = list(size=8, color ="black"),
                 repel=TRUE, label.rectangle = FALSE,  color="significant",
                 palette=c(down = "#003C67FF", ns = "#868686FF", up = "#A73030FF"),
                 xlab="log2 fold change MIB binding",
                 ylab="-log10 p-value",
                 #xticks.by=4, ###Set x-axis tick intervals
                 #xlim=c(-7.5,7.5),        ###Set limits of x axis
                 #ylim=c(0,4),         ###set limits of y axis
                 legend="none")

    p<-p+geom_hline(yintercept = 1.3, linetype = 2) +
      geom_vline(xintercept = 1, linetype=2) +
      geom_vline(xintercept = -1, linetype=2) +
      grids(linetype="dashed")+border()

    ggsave(p, file=paste0(A," vs ", B,".svg"),device = "svg", width=5, height=5)
    ggsave(p, file=paste0(A," vs ", B,".pdf"), width=5, height=5)

  }

  return(all_kinase)
}
