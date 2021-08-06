library(DESeq2)
library(pheatmap)
library(CEMiTool)
library(tidyverse)
library(data.table)
library(STRINGdb)
library(ggrepel)

cem.gmt<-gmt_input
cem.metadata<-samples_metadata
colnames(cem.metadata[,3:12])<-Tissues[1:10]
cem.ppi<-ppi1.hc700[, 4:5]
cem.count.data<-cem.data2.df[, c(2,11:30)]


#Creating sample annotation from metadata:
sample_annot.cem <- data.frame(cem.metadata$Samples, cem.metadata$TISSUES)
colnames(sample_annot.cem)<-  c("SampleName", "Class")

#Creating the data.frame with counts 
count.cem <-data.frame(data.table(counts[,cem.metadata$Samples]))
rownames(count.cem) <- make.names(cem.count.data$GeneID, unique=TRUE)


foreach::registerDoSEQ()



cem <- CEMiTool::cemitool(count.cem, sample_annot.cem, 
                          interactions = cem.ppi,cem.gmt,
                          filter=TRUE, plot=TRUE, 
                          verbose=TRUE, force_beta = TRUE, min_ngen = 6)

cem <- mod_ora(cem, cem.gmt)
generate_report(cem, force=T, output_format=c("pdf_document", "html_document"),
                title = "Differential Co-expression analysis in Pikeperch")


#Output is a csv dataframe containing all genes that are associated with a co-expression module
fwrite(CEMiTool::module_genes(cem)[order(CEMiTool::module_genes(cem)$modules),], file="Julien_genesallmodulesDiagnose.csv", quote=F, sep="\t")

#Gene Set Enrichment Analysis in a more fashion way with pheatmap:
NESheatmap <- function(CEMiToolfile=cem){
  NES <- t(as.matrix(CEMiToolfile@enrichment$nes[,-1]))
  NES[is.na(NES)] <- 0
  colnames(NES) <- CEMiToolfile@enrichment$nes[,1]
  #ordering:
  NES <- NES[,order(nchar(colnames(NES)))]
  heatmap <- pheatmap::pheatmap(
    mat=NES,
    color= viridis::viridis(length(seq(min(NES, na.rm = T), max(NES,na.rm = T), length.out = 10)) -1),
    show_colnames     = TRUE,
    show_rownames     = TRUE,
    cluster_cols = F,
    fontsize = 18,
    main="GSEA of co-expression modules",
    legend_breaks = c(-4, -2, 0, 2, 4, max(NES)),
    legend_labels= c("-4","-2","0","2","4", "NES\n")                    
  )
  return(heatmap)
}
cemheatmap <- NESheatmap(cem)
ggsave(plot=cemheatmap, filename=paste0("GSEA","_","heatmap",".","png"),
       dpi="retina",width = 8, height =8, units = "in", device="png")

cem <- plot_interactions(cem)
cem.plots <- show_plot(cem, "interaction")
cem.plots[6]

