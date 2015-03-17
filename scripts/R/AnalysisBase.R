#{{{
# Analysis with edgeR.
# source("/nfs/omics_scratch/PROJECTS/Chi/Doss/20140923/scripts/R/AnalysisBase.R")
#}}}

#{{{
# Data structure:
#   environmetn
#     --- Fields that should be set by invokin script
#       projdir     : string -- full path to project directory
#       datadir     : string -- <projdir>/processedData
#       htseqdir    : string -- <datadir>/d5_HTSeq
#       gseadir     : string -- <datadir>/d6_GSEAinput
#       resultdir   : string -- <projdir>/results
#       sl_file     : string -- full path to SampleList.txt file
#       gn_file     : string -- full path to gene_names.txt file
#       cg_file     : string -- full path to coding_genes.txt file
#     --- Fields that are created by functions in this script
#       samples     : data.frame from sl_file, "group" column added
#       gene_names  : vector<string>, values==gene-names, names==ens-gene-ids
#       coding_genes: vector<string>, values==ens-gene-ids
#       unique_counts   : ngenes x nsamples matrix of counts, rownames==GIDs, colnames==SampleNames
#       multi_counts    : ngenes x nsamples matrix of counts, rownames==GIDs, colnames==SampleNames

#       raw_counts  : list(data.frame), names==sample-names,
#                       df-colnames==[FeatureId,All,Multi], df-rownames==ens-gene-ids
#     --- Various analyses are build by both invoking and this scripts
#       analyses : environment
#           analysis-name : environment
#             --- Fields prepared by invoking script
#               counts_matrix
#               samples : subset of global samples
#               desing  : model.matrix(...)
#               contrasts  : environment
#                   contrast-name : environment
#                     --- Fields prepared by invoking script
#                       contr : makeContrast(...)
#                       data_idx : vector of indices of columns of 'counts_matrix' that should
#                                  be included in top table
#                     --- Fields that are created by function makeAnalysis() in this script
#                       fit2 : glmLTR(...)
#                       tt : the top table, ready to save
#             --- Fields that are created by function makeAnalysis() in this script
#               dge     : DGEList, with all dispersions
#               norm_counts : matrix, columns==samples, rows==genes, cpm(dge)
#               pca     : PCA data
#               hcdist  : distances for Hierarchical Clustering
#               hc      : Hierarchical Clustering
#               ph      : as.phylo(hc)
#               glmfit  : glmFit(dge,design)
#}}}

library(edgeR)
library(splines)
library(gplots)
library(ape)
source(file.path(Sys.getenv("PROJDIR"), "scripts", "R", "plotphylo.R"))
# source: the patched version of plot.phylo

mkDir <- function(path) { dir.create(path, recursive=T, showWarnings=F) }

loadSampleList <- function(dd) { #{{{
    cat("--- Loading sample list...");  flush(stdout())
    dd$samples <- read.delim( dd$sl_file, 
                              stringsAsFactors=FALSE, check.names=F, colClasses="character" )
    # We need to create 'group' column in samples;
    # the following code hanles R inconsistencies with handling data frames with 0/1/many columns.
    cns <- colnames(dd$samples)
    groups <- rep("g", nrow(dd$samples))
    for (cn in cns[cns!="SampleName" & cns!="Replicate"]) {
        groups <- paste(groups, dd$samples[[cn]], sep="_")
    }
    dd$samples$group <- groups
    cat("Done.\n")
}#}}}

loadGeneNames <- function(dd) { #{{{
    cat("--- Loading gene names...");  flush(stdout())
    gnames_df <- read.delim( dd$gn_file, 
                             stringsAsFactors=FALSE, check.names=F, colClasses="character",
                             header=F, col.names=c("ETID", "EGID", "SYMB") )
    uidx <- ! duplicated(gnames_df$EGID)
    dd$gene_names <- gnames_df$SYMB[uidx]
    names(dd$gene_names) <- gnames_df$EGID[uidx]
    cat("Done.\n")
}#}}}

loadCodingGenes <- function(dd) { #{{{
    cat("--- Loading coding genes...");  flush(stdout())
    cgenes_df <- read.delim( dd$cg_file, 
                             stringsAsFactors=FALSE, check.names=F, colClasses="character",
                             header=F, col.names=c("EGID") )
    dd$coding_genes <- cgenes_df$EGID
    cat("Done.\n")
}#}}}

loadData <- function(dd) { #{{{
    cat("--- Loading counts.")
    counts_list <- list()
    for (sname in dd$samples$SampleName) {
        counts_fn = file.path(dd$htseqdir, paste0(sname,"_htseq.count"))
        counts_list[[sname]] <- read.delim( counts_fn, header=F, stringsAsFactors=F,
                                            col.names=c("FeatureId","All") )
        cat('.');  flush(stdout())
    }
    cat("Done.\n")
    if (length(counts_list)==0) {
        cat("ERROR: no count file is loaded.\n")
        return()
    }
    cat("--- Merging counts.");  flush(stdout())
    count_rownames <- counts_list[[1]]$FeatureId
    count_generows <- which(substr(count_rownames,1,3)=="ENS")
    dd$multi_counts = matrix(nrow=length(count_generows), ncol=nrow(dd$samples))
    rownames(dd$multi_counts) <- counts_list[[1]]$FeatureId[count_generows]
    colnames(dd$multi_counts) <- dd$samples$SampleName
    for (i in 1:length(counts_list)) {
        dd$multi_counts[,i] <- counts_list[[i]]$All[count_generows]
        cat('.');  flush(stdout())
    }
    dd$unique_counts <- dd$multi_counts
    cat("Done.\n")
}#}}}
loadDataPlus <- function(dd) { #{{{
    cat("--- Loading counts.");  flush(stdout())
    counts_list <- list()
    for (sname in dd$samples$SampleName) {
        counts_fn = file.path(dd$htseqdir, paste0(sname,"_htseq.count"))
        counts_list[[sname]] <- read.delim( counts_fn, header=T, stringsAsFactors=F )
        cat('.');  flush(stdout())
    }
    cat("Done.\n")
    if (length(counts_list)==0) {
        cat("ERROR: no count file is loaded.\n")
        return()
    }
    cat("--- Merging counts.");  flush(stdout())
    count_rownames <- counts_list[[1]]$FeatureId
    count_generows <- which(substr(count_rownames,1,3)=="ENS")
    dd$multi_counts = matrix(nrow=length(count_generows), ncol=nrow(dd$samples))
    rownames(dd$multi_counts) <- counts_list[[1]]$FeatureId[count_generows]
    colnames(dd$multi_counts) <- dd$samples$SampleName
    dd$unique_counts = matrix(nrow=length(count_generows), ncol=nrow(dd$samples))
    rownames(dd$unique_counts) <- counts_list[[1]]$FeatureId[count_generows]
    colnames(dd$unique_counts) <- dd$samples$SampleName
    nn <- ncol(counts_list[[1]])
    for (i in 1:length(counts_list)) {
        dd$multi_counts[,i] <- counts_list[[i]]$All[count_generows]
        dd$unique_counts[,i] <- counts_list[[i]]$All[count_generows] -
                                    rowSums(counts_list[[i]][count_generows,3:nn])
        cat('.');  flush(stdout())
    }
    dd$counts_list <- counts_list
    cat("Done.\n")
}#}}}
    

makeAnalysis <- function(dd, aname) { #{{{
    cat(sprintf("Performing analysis %s:\n", aname));  flush(stdout())
    dda <- dd$analyses[[aname]]
    
    cat("--- Creating DGEList.");  flush(stdout())
    counts_matrix <- dda$counts_matrix[ -which(apply(dda$counts_matrix,1,max)<10), ]
    cat(".");  flush(stdout())
    dge <- DGEList(counts=counts_matrix, group=dda$samples$group)
    cat(".");  flush(stdout())
    dge <- calcNormFactors(dge);                        cat(".");  flush(stdout())
    dge <- estimateGLMCommonDisp(dge, dda$design);      cat(".");  flush(stdout())
    dge <- estimateGLMTrendedDisp(dge, dda$design);     cat(".");  flush(stdout())
    dge <- estimateGLMTagwiseDisp(dge, dda$design);     cat(".");  flush(stdout())
    dda$dge <- dge
    dda$norm_counts <- cpm(dge)
    cat("Done.\n");  flush(stdout())
    
    cat("--- Principal component analysis...");  flush(stdout())
    dda$pca <- prcomp(t(dda$norm_counts))
    cat("Done.\n");  flush(stdout())
    cat("--- Hierarhical clustering...");  flush(stdout())
    dda$hcdist <- as.dist( (1 - cor(dda$norm_counts))/2 )
    dda$hc <- hclust(dda$hcdist)
    dda$ph <- as.phylo(dda$hc)
    cat("Done.\n");  flush(stdout())
    
    cat("--- Fitting model...");  flush(stdout())
    dda$glmfit <- glmFit(dda$dge, dda$design)
    cat("Done.\n");  flush(stdout())
    
    cat("--- Building top tables for contrasts:\n");  flush(stdout())
    dda$results <- new.env()
    for (contr in ls(dda$contrasts)) {
        cat(sprintf("      %s", contr));  flush(stdout())
        ddc <- dda$contrasts[[contr]]
        ddc$fit2 <- glmLRT(dda$glmfit, contrast=ddc$contr)
        ddc$tt <- topTags(ddc$fit2, n=nrow(dda$dge), sort.by="none")$table
        ddc$tt$GeneID <- row.names(ddc$tt)
        ddc$tt$GeneName <- dd$gene_names[ddc$tt$GeneID]
        ddc$tt$IsCoding <- "FALSE"
        ddc$tt$IsCoding[ which(ddc$tt$GeneID %in% dd$coding_genes) ] <- "TRUE"
        ddc$tt <- ddc$tt[, c("GeneID","GeneName","IsCoding","logFC","logCPM","LR","PValue","FDR")]
        ddc$tt <- merge(ddc$tt,dda$norm_counts[,ddc$data_idx],by.x="GeneID",by.y="row.names",all=T)
        ddc$tt <- ddc$tt[order(ddc$tt$PValue),]
        nr_top_genes <- sum(ddc$tt$FDR < 0.1)
        cat(sprintf(" - done, number of top genes (FDR<0.1) = %d.\n",nr_top_genes)); flush(stdout())
    }
    cat(sprintf("Done analysis %s.\n", aname));  flush(stdout())
}#}}}

# QC plost {{{
# Multidimensional scaling plot
plotMS <- function(dd, aname, colored="group", components=c(1,2)) { #{{{
    dda <- dd$analyses[[aname]]
    plotMDS(dda$dge, top=nrow(dda$dge),
            col=as.numeric(as.factor(dda$samples[[colored]])), dim.plot=components,
            main=sprintf("Multidimensional scaling plot\n(colored: %s)", colored) )
}#}}}

# PCA plot
plotPCA <- function(dd, aname, colored="group", components=c("PC1","PC2")) { #{{{
    dda <- dd$analyses[[aname]]
    xdata <- as.data.frame(dda$pca$x)[[components[1]]]
    ydata <- as.data.frame(dda$pca$x)[[components[2]]]
    plot( xdata, ydata, col=as.numeric(as.factor(dda$samples[[colored]])),
          pch=19, xlab=components[1], ylab=components[2],
          main=sprintf("Principal components plot\n(colored: %s)", colored) )
    text( xdata, ydata, dda$samples$SampleName, pos=1, cex=0.7)
}#}}}

# HC plot
plotHC <- function(dd, aname, colored="group") { #{{{
    # Here, the pathced version will be invoked -- should have correct font
    dda <- dd$analyses[[aname]]
    plot(dda$ph,
         label.offset=0.0001,       # what are units???
         tip.color=as.numeric(as.factor(dd$samples[[colored]]))
        )
    title(main=sprintf("Sample clusters by correlation distance\n(colored: %s)", colored))
}#}}}

# Dispersion Estimates (biological coefficient of variation) plot
plotDE <- function(dd, aname) { #{{{
    plotBCV(dd$analyses[[aname]]$dge, main="Per-gene biological coefficient of variation")
}#}}}

saveQCPlots <- function(dd, aname, colored="group") { #{{{
    save_dir <- file.path(dd$resultdir, aname)
    mkDir(save_dir)
    fpath <- file.path(save_dir, "QCPlots.pdf")
    pdf(fpath, width=8.5, height=11)
    plotMS(dd,aname,colored)
    plotPCA(dd,aname,colored)
    plotHC(dd,aname,colored)
    plotDE(dd,aname)
    dev.off()
}#}}}
#}}}

saveTopTable <- function(dd, aname, contrast) { #{{{
    save_dir <- file.path(dd$resultdir, aname)
    mkDir(save_dir)
    fpath <- file.path(save_dir, paste0(contrast,".csv"))
    write.table( dd$analyses[[aname]]$contrasts[[contrast]]$tt,
                 file=fpath, quote=FALSE, row.names=FALSE, sep="\t" )
}#}}}

saveGseaInput <- function(dd, aname, contrast) { #{{{
    save_dir <- file.path(dd$resultdir, aname)
    mkDir(save_dir)
    fpath <- file.path(save_dir, paste0(contrast,".rnk"))
    
    gsea_input <- subset( dd$analyses[[aname]]$contrasts[[contrast]]$tt, select=c(GeneName,PValue) )
    gsea_input <- gsea_input[ !is.na(gsea_input$PValue), ]
    gsea_input$PValue <- 1 - gsea_input$PValue
    gsea_input$GeneName <- toupper(gsea_input$GeneName)
    write.table( gsea_input, file=fpath, quote=FALSE, row.names=F, col.names=F, sep="\t" )
}#}}}
