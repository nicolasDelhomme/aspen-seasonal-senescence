#' ---
#' title: "Differential Expression"
#' author: "Nicolas Delhomme & Jenna Lihavainen"
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
suppressMessages({
    source(here("UPSCb-common/src/R/volcanoPlot.R"))
    source(here("UPSCb-common/src/R/featureSelection.R"))
})


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
               aes(x=Time,y=value,col=Genotype,group=Genotype)) +
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
load(here("data/analysis/salmon/dds_subsampling_no_outliers.rda"))

#' ## Normalisation for visualisation
vsd <- varianceStabilizingTransformation(dds,blind=FALSE)
vst <- assay(vsd)
vst <- vst - min(vst)
write_tsv(as.data.frame(vst_gen1) %>% rownames_to_column("ID"),here("data/analysis/salmon/variance-stabilised_model-aware_data.tsv"))


#' ## Gene of interest
#' The gene is not expressed
line_plot(dds,vst,"Potra2n4c10046")

#' The gene expression levels differ between the genotypes
line_plot(dds,vst,"Potra2n19c33539")

#' The gene displays enhanced expression during autumn
line_plot(dds,vst,"Potra2n15c28290")

#' The gene displays repressed expression during autumn
line_plot(dds,vst,"Potra2n2c4558")

#' The gene displays different temporal expression pattern between the genotypes
line_plot(dds,vst,"Potra2n14c27101")

#' ## Differential Expression
dds <- DESeq(dds)

#' * Dispersion estimation
#' The dispersion estimation is adequate
plotDispEsts(dds)

#' #' The model used is:
#' #' `Genotype` and `Time`. Time effect is considered within whole data and within individual genotypes. 

ddsTime <- dds
design(ddsTime) = ~ Time
ddsTime$Time <- relevel(ddsTime$Time,"225")
ddsTime <- DESeq(ddsTime)
resultsNames(ddsTime)

# time effect 264 vs 225 DOY (Senescence-associated genes)
dir.create(here("../../../../Git/aspen-seasonal-senescence/data/analysis/DE/"), showWarnings = FALSE)

T264vs225 <- extract_results(ddsTime,vst,"Time_264_vs_225",
                             default_dir = here("data/analysis/DE/Time"),
                             default_prefix = "T264-225_",
                             labels=paste(ddsTime$Time,sep="_"),
                             sample_sel = ddsTime$Time %in% c("225","264"),
                             plot = TRUE,verbose = TRUE)

# time effects (consecutive time points)
dir.create(here("data/analysis/DE/"),showWarnings = FALSE)
T232vs225 <- extract_results(ddsTime,vst,"Time_232_vs_225",
                             default_dir = here("data/analysis/DE/Time"),
                             default_prefix = "T232-225_",
                             labels=paste(ddsTime$Time,sep="_"),
                             sample_sel = ddsTime$Time %in% c("225","232"),
                             plot = TRUE,verbose = TRUE)
    
T237vs232 <- extract_results(ddsTime,vst,list("Time_237_vs_225","Time_232_vs_225"),
                             default_dir = here("data/analysis/DE/Time"),
                             default_prefix = "T237-232_",
                             labels=paste(ddsTime$Time,sep="_"),
                             sample_sel = ddsTime$Time %in% c("232","237"),
                             plot = TRUE,verbose = TRUE)

T241vs237 <- extract_results(ddsTime,vst,list("Time_241_vs_225","Time_237_vs_225"),
                             default_dir = here("data/analysis/DE"),
                             default_prefix = "T241-237_",
                             labels=paste(ddsTime$Time,sep="_"),
                             sample_sel = ddsTime$Time %in% c("241","237"),
                             plot = TRUE,verbose = TRUE)

T246vs241 <- extract_results(ddsTime,vst,list("Time_246_vs_225","Time_241_vs_225"),
                             default_dir = here("data/analysis/DE/Time"),
                             default_prefix = "T246-241_",
                             labels=paste(ddsTime$Time,sep="_"),
                             sample_sel = ddsTime$Time %in% c("246","241"),
                             plot = TRUE,verbose = TRUE)

T250vs246 <- extract_results(ddsTime,vst,list("Time_250_vs_225","Time_246_vs_225"),
                             default_dir = here("data/analysis/DE/Time"),
                             default_prefix = "T250-246_",
                             labels=paste(dds$Genotype,ddsTime$Time,sep="_"),
                             sample_sel = ddsTime$Time %in% c("250","246"),
                             plot = TRUE,verbose = TRUE)

T253vs250 <- extract_results(ddsTime,vst,list("Time_253_vs_225","Time_250_vs_225"),
                             default_dir = here("data/analysis/DE/Time"),
                             default_prefix = "T253-250_",
                             labels=paste(dds$Genotype,ddsTime$Time,sep="_"),
                             sample_sel = ddsTime$Time %in% c("253","250"),
                             plot = TRUE,verbose = TRUE)

T256vs253 <- extract_results(ddsTime,vst,list("Time_256_vs_225","Time_253_vs_225"),
                             default_dir = here("data/analysis/DE/Time"),
                             default_prefix = "T256-253_",
                             labels=paste(ddsTime$Time,sep="_"),
                             sample_sel = ddsTime$Time %in% c("256","253"),
                             plot = TRUE,verbose = TRUE)

T264vs256 <- extract_results(ddsTime,vst,list("Time_264_vs_225","Time_256_vs_225"),
                             default_dir = here("data/analysis/DE/Time"),
                             default_prefix = "T264-256_",
                             labels=paste(ddsTime$Time,sep="_"),
                             sample_sel = ddsTime$Time %in% c("264","256"),
                             plot = TRUE,verbose = TRUE)

# Genotype effect 

ddsGenotype <- dds
design(ddsGenotype) = ~ Genotype
ddsGenotype$Genotype <- relevel(ddsGenotype$Genotype,"1")
ddsGenotype <- DESeq(ddsGenotype)

Gen1vs48 <- extract_results(ddsGenotype,vst,"Genotype_48_vs_1",
                             default_dir = here("data/analysis/DE/Genotype"),
                             default_prefix = "Gen1-gen48_",
                             labels=paste(dds$Genotype,sep="_"),
                             sample_sel = ddsGenotype$Genotype %in% c("48","1"),
                             plot = TRUE,verbose = TRUE)

Gen1vs81 <- extract_results(ddsGenotype,vst,"Genotype_81_vs_1",
                               default_dir = here("data/analysis/DE/Genotype"),
                               default_prefix = "Gen1-gen81_",
                               labels=paste(dds$Genotype,sep="_"),
                               sample_sel = ddsGenotype$Genotype %in% c("81","1"),
                               plot = TRUE,verbose = TRUE)

Gen48vs81 <- extract_results(ddsGenotype,vst,list("Genotype_48_vs_1","Genotype_81_vs_1"),
                             default_dir = here("data/analysis/DE/Genotype"),
                             default_prefix = "Gen48-gen81_",
                             labels=paste(dds$Genotype,sep="_"),
                             sample_sel = ddsGenotype$Genotype %in% c("81","48"),
                             plot = TRUE,verbose = TRUE)


#Time effect in genotype 1

dds_gen1 <- dds[, dds$Genotype %in% c("1")]
dds_gen1$Genotype <- droplevels(dds_gen1$Genotype)
design(dds_gen1) = ~ Time
dds_gen1$Time <- relevel(dds_gen1$Time,"225")
dds_gen1 <- DESeq(dds_gen1)

# time effect 264 vs 225 (Senescence associated genes-SAGs)

T264vs225 <- extract_results(dds_gen1,vst,"Time_264_vs_225",
                             default_dir = here("data/analysis/DE/Time/Gen1"),
                             default_prefix = "T264-225_",
                             sample_sel = dds_gen1$Time %in% c("225","264"),
                             plot = TRUE,verbose = TRUE)

# time effect (compared to 225)

T256vs225 <- extract_results(dds_gen1,vst,"Time_256_vs_225",
                             default_dir = here("data/analysis/DE/Time/Gen1"),
                             default_prefix = "T256-225_",
                             sample_sel = ddsTime$Time %in% c("225","256"),
                             plot = TRUE,verbose = TRUE)

T253vs225 <- extract_results(dds_gen1,vst,"Time_253_vs_225",
                             default_dir = here("data/analysis/DE/Time/Gen1"),
                             default_prefix = "T253-225_",
                             sample_sel = ddsTime$Time %in% c("225","253"),
                             plot = TRUE,verbose = TRUE)

T250vs225 <- extract_results(dds_gen1,vst,"Time_250_vs_225",
                             default_dir = here("data/analysis/DE/Time/Gen1"),
                             default_prefix = "T250-225_",
                             sample_sel = ddsTime$Time %in% c("225","250"),
                             plot = TRUE,verbose = TRUE)


T246vs225 <- extract_results(dds_gen1,vst,"Time_246_vs_225",
                             default_dir = here("data/analysis/DE/Time/Gen1"),
                             default_prefix = "T246-225_",
                             sample_sel = ddsTime$Time %in% c("225","246"),
                             plot = TRUE,verbose = TRUE)

T241vs225 <- extract_results(dds_gen1,vst,"Time_241_vs_225",
                             default_dir = here("data/analysis/DE/Time/Gen1"),
                             default_prefix = "T241-225_",
                             sample_sel = ddsTime$Time %in% c("225","241"),
                             plot = TRUE,verbose = TRUE)

T237vs225 <- extract_results(dds_gen1,vst,"Time_237_vs_225",
                             default_dir = here("data/analysis/DE/Time/Gen1"),
                             default_prefix = "T237-225_",
                             labels=paste(dds$Genotype,ddsTime$Time,sep="_"),
                             sample_sel = ddsTime$Time %in% c("225","237"),
                             plot = TRUE,verbose = TRUE)


#' # Session Info 
#'  ```{r session info, echo=FALSE}
#'  sessionInfo()
#'  ```


