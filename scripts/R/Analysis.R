#{{{
# source("/nfs/omics_scratch/PROJECTS/Chi/Doss/20140923/scripts/R/Analysis.R")
#}}}

source(file.path(Sys.getenv("PROJDIR"), "scripts", "R", "AnalysisBase.R"))

doPrepare <- function() { #{{{
    dd <<- new.env()
    dd$projdir <- Sys.getenv("PROJDIR")
    dd$datadir <- Sys.getenv("procdata_dir")
    dd$htseqdir <- Sys.getenv("htseq_dir")
    dd$gseadir <- Sys.getenv("gsea_dir")
    dd$resultdir <- Sys.getenv("results_dir")
    
    dd$sl_file <- file.path(dd$datadir, "SampleDescription.csv")
    dd$gn_file <- Sys.getenv("genenamesfile")
    dd$cg_file <- Sys.getenv("codinggenesfile")
    
    dd$analyses <- new.env()
}#}}}

doLoad <- function() { #{{{
    cat("Loading data:\n")
    loadSampleList(dd)
    loadGeneNames(dd)
    loadCodingGenes(dd)
    loadDataPlus(dd)
    cat("Done loading data.\n")
}#}}}

doAnalysis <- function(aname, counts_matrix) { #{{{
    dd$analyses[[aname]] <- new.env()
    dda <- dd$analyses[[aname]]
    dda$counts_matrix <- counts_matrix
    dda$samples <- dd$samples
    dda$samples$Phenotype <- relevel( as.factor(dda$samples$Phenotype), "HbAA" )
    dda$design <- model.matrix(~ Phenotype, data=dda$samples)
    colnames(dda$design)[1] <- "Intercept"
    colnames(dda$design) <- sub("Phenotype","", colnames(dda$design))
    
    makeDataIdx <- function(ph1,ph2) {
        isamples <- c(1:nrow(dda$samples))
        idx1 <- isamples[ dda$samples$Phenotype == ph1 ]
        idx2 <- isamples[ dda$samples$Phenotype == ph2 ]
        c(idx1,idx2)
    }
    dda$contrasts <- new.env()
    ddc <- dda$contrasts
    ddc$SSvsAA <- new.env()
    ddc$SSvsAA$contr      <- makeContrasts(SSvsAA=HbSS, levels=dda$design)
    ddc$SSvsAA$data_idx   <- makeDataIdx("HbSS", "HbAA")
    
    makeAnalysis(dd, aname)
}#}}}
doAnalyses <- function() { #{{{
    doAnalysis("Unique", dd$unique_counts)
    doAnalysis("Multi", dd$multi_counts)
}#}}}

doSaveResults <- function() { #{{{
    for (aname in ls(dd$analyses)) {
        saveQCPlots(dd, aname)
        for (contr in ls(dd$analyses[[aname]]$contrasts)) {
            saveTopTable(dd, aname, contr)
            #saveGseaInput(dd, aname, contr)
        }
    }
}#}}}

doAll <- function() { #{{{
    doPrepare()
    doLoad()
    doAnalyses()
    doSaveResults()
}#}}}
