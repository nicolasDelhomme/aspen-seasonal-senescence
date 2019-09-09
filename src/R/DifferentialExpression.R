#' ---
#' title: "Differential Expression"
#' author: "Nicolas Delhomme & Domenique Andre"
#' date: "`r Sys.Date()`"
#' output:
#'  html_document:
#'    toc: true
#'    number_sections: true
#' ---
#' # Setup
#' * Working directory

#' * Libraries
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(DESeq2))
suppressPackageStartupMessages(library(gplots))
suppressPackageStartupMessages(library(here))
suppressPackageStartupMessages(library(hyperSpec))
suppressPackageStartupMessages(library(RColorBrewer))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(VennDiagram))

#' * Helper files
suppressMessages(source("~/Git/UPSCb/src/R/plotMA.R"))
suppressMessages(source("~/Git/UPSCb/src/R/volcanoPlot.R"))

#' * Graphics
pal=brewer.pal(8,"Dark2")
hpal <- colorRampPalette(c("blue","white","red"))(100)
mar <- par("mar")

#' * Functions
#' 1. plot specific gene expression
"line_plot" <- function(dds,vst,gene_id){
    sel <- grepl(gene_id,rownames(vst))
    stopifnot(sum(sel)==1)

    return(
        ggplot(bind_cols(as.data.frame(colData(dds)),
                         melt(vst[sel,])),
               aes(x=Time,y=value,col=Experiment,group=Experiment)) +
            geom_point() + geom_smooth() +
            scale_y_continuous(name="VST expression") + 
            ggtitle(label=paste("Expression for: ",gene_id))
    )
}

#' 2. extract the DE results. Default cutoffs are
#' from Schurch _et al._, RNA, 2016
"extract_results" <- function(dds,vst,contrast,
                              padj=0.01,lfc=0.5,
                              plot=TRUE,verbose=TRUE,
                              export=TRUE,default_dir="analysis/DE",
                              default_prefix="DE-",
                              labels=colnames(dds),
                              sample_sel=1:ncol(dds)){
    
    if(length(contrast)==1){
        res <- results(dds,name=contrast)
    } else {
        res <- results(dds,contrast=contrast)
    }
    
    if(plot){
        par(mar=c(5,5,5,5))
        volcanoPlot(res)
        par(mar=mar)
    }
    
    sel <- res$padj <= padj & abs(res$log2FoldChange) >= lfc & ! is.na(res$padj)
    
    if(verbose){
        message(sprintf("There are %s genes that are DE",sum(sel)))
    }
            
    if(export){
        if(!dir.exists(default_dir)){
            dir.create(default_dir,showWarnings=FALSE,recursive=TRUE,mode="0771")
        }
        write.csv(res,file=file.path(default_dir,paste0(default_prefix,"results.csv")))
        write.csv(res[sel,],file.path(default_dir,paste0(default_prefix,"genes.csv")))
    }
    if(plot){
        heatmap.2(t(scale(t(vst[sel,sample_sel]))),
                  distfun = pearson.dist,
                  hclustfun = function(X){hclust(X,method="ward.D2")},
                  trace="none",col=hpal,labRow = FALSE,
                  labCol=labels[sample_sel]
        )
    }
    return(rownames(res[sel,]))
}

#' # _Populus tremula_
#' * Data
load(here("../../analysis/salmon/Potra-104-127-139-removed-dds.rda"))

#' ## Normalisation for visualisation
vsd <- varianceStabilizingTransformation(dds,blind=FALSE)
vst <- assay(vsd)
vst <- vst - min(vst)

#' ## Gene of interest
#' * Potra000962g07909
#' The gene is not expressed
line_plot(dds,vst,"Potra000962g07909")

#' * Potra001661g13641
#' The gene is  medium expressed and show time-lagged expression patterns between
#' the experiments
line_plot(dds,vst,"Potra001661g13641")

#' * Potra003711g22520
#' The gene is medium expressed and show initially opposing expression patterns between experiments
line_plot(dds,vst,"Potra003711g22520")

#' * Potra006413g25676
#' The gene his very lowly expressed and shows somewhat opposite expression patterns
line_plot(dds,vst,"Potra006413g25676")

#' ## Differential Expression
dds <- DESeq(dds)

#' * Dispersion estimation
#' The dispersion estimation is adequate
plotDispEsts(dds)

#' The model used is:
#' 
#' `Experiment * Time` meaning that the `Experiment` and `Time variable` as 
#' well as their interaction `Experiment:Time` is considered. Because we 
#' cannot assume that these two variables explain all the variance in the data,
#' there is also an `Intercept` for the linear model. This also implies that the 
#' model assumes `Cont` at `3` hours to be the baseline; _i.e._ everything is compared 
#' against it.
resultsNames(dds)

#' ## Results
#' In the following we look at the interaction specific genes; _i.e._ genes that 
#' changes at a given time transition in between experiments
#' ### ECM _vs._ Cont at T3
Pa_3 <- extract_results(dds,vst,"Experiment_ECM_vs_Cont",
                        default_prefix="Potra_ECM-vs-Cont_T3_",
                        labels=paste0(colData(dds)$Experiment,
                                      colData(dds)$Time),
                        sample_sel=colData(dds)$Time==3)

#' ### ECM _vs._ Cont at T7
Pa_7 <- extract_results(dds,vst,c(0,1,0,0,0,0,1,0,0,0),
                        default_prefix="Potra_ECM-vs-Cont_T7_",
                        labels=paste0(colData(dds)$Experiment,
                                      colData(dds)$Time),
                        sample_sel=colData(dds)$Time==7)

#' ### ECM _vs._ Cont at T14
Pa_14 <- extract_results(dds,vst,c(0,1,0,0,0,0,0,1,0,0),
                         default_prefix="Potra_ECM-vs-Cont_T14_",
                         labels=paste0(colData(dds)$Experiment,
                                       colData(dds)$Time),
                         sample_sel=colData(dds)$Time==14)

#' ### ECM _vs._ Cont at T21
Pa_21 <- extract_results(dds,vst,c(0,1,0,0,0,0,0,0,1,0),
                         default_prefix="Potra_ECM-vs-Cont_T21_",
                         labels=paste0(colData(dds)$Experiment,
                                       colData(dds)$Time),
                         sample_sel=colData(dds)$Time==21)

#' ### ECM _vs._ Cont at T28
Pa_28 <- extract_results(dds,vst,c(0,1,0,0,0,0,0,0,0,1),
                         default_prefix="Potra_ECM-vs-Cont_T28_",
                         labels=paste0(colData(dds)$Experiment,
                                       colData(dds)$Time),
                         sample_sel=colData(dds)$Time==28)

#' ### Venn Diagram
grid.newpage()
grid.draw(venn.diagram(list(T3=Pa_3,
                            T7=Pa_7,
                            T14=Pa_14,
                            T21=Pa_21,
                            T28=Pa_28),
                       NULL,
                       fill=pal[1:5]))

#' # Session Info 
#'  ```{r session info, echo=FALSE}
#'  sessionInfo()
#'  ```


