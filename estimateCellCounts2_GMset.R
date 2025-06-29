validationCellType <- function(Y, pheno, modelFix, modelBatch = NULL,
                               L.forFstat = NULL, verbose = FALSE) {
  N <- dim(pheno)[1]
  pheno$y <- rep(0, N)
  xTest <- model.matrix(modelFix, pheno)
  sizeModel <- dim(xTest)[2]
  M <- dim(Y)[1]
  if (is.null(L.forFstat)) {
    L.forFstat <- diag(sizeModel)[-1, ] # All non-intercept coefficients
    colnames(L.forFstat) <- colnames(xTest)
    rownames(L.forFstat) <- colnames(xTest)[-1]
  }
  # Initialize various containers
  sigmaResid <- sigmaIcept <- nObserved <- nClusters <- Fstat <- rep(NA, M)
  coefEsts <- matrix(NA, M, sizeModel)
  coefVcovs <- list()
  if (verbose) {
    cat("[validationCellType] ")
  }
  for (j in seq_len(M)) { # For each CpG
    ## Remove missing methylation values
    ii <- !is.na(Y[j, ])
    nObserved[j] <- sum(ii)
    pheno$y <- Y[j, ]
    
    if (j %% round(M / 10) == 0 && verbose) {
      cat(".")
    } # Report progress
    
    try({ # Try to fit a mixed model to adjust for plate
      if (!is.null(modelBatch)) {
        fit <- try(lme(modelFix, random = modelBatch, data = pheno[ii, ]))
        OLS <- inherits(fit, "try-error")
        # If LME can't be fit, just use OLS
      } else {
        OLS <- TRUE
      }
      
      if (OLS) {
        fit <- lm(modelFix, data = pheno[ii, ])
        fitCoef <- fit$coef
        sigmaResid[j] <- summary(fit)$sigma
        sigmaIcept[j] <- 0
        nClusters[j] <- 0
      } else {
        fitCoef <- fit$coef$fixed
        sigmaResid[j] <- fit$sigma
        sigmaIcept[j] <- sqrt(getVarCov(fit)[1])
        nClusters[j] <- length(fit$coef$random[[1]])
      }
      coefEsts[j, ] <- fitCoef
      coefVcovs[[j]] <- vcov(fit)
      
      useCoef <- L.forFstat %*% fitCoef
      useV <- L.forFstat %*% coefVcovs[[j]] %*% t(L.forFstat)
      Fstat[j] <- (t(useCoef) %*% solve(useV, useCoef)) / sizeModel
    })
  }
  if (verbose) {
    cat(" done\n")
  }
  ## Name the rows so that they can be easily matched to the target data set
  rownames(coefEsts) <- rownames(Y)
  colnames(coefEsts) <- names(fitCoef)
  degFree <- nObserved - nClusters - sizeModel + 1
  
  ## Get P values corresponding to F statistics
  Pval <- 1 - pf(Fstat, sizeModel, degFree)
  
  out <- list(
    coefEsts = coefEsts, coefVcovs = coefVcovs, modelFix = modelFix,
    modelBatch = modelBatch,
    sigmaIcept = sigmaIcept, sigmaResid = sigmaResid,
    L.forFstat = L.forFstat, Pval = Pval,
    orderFstat = order(-Fstat), Fstat = Fstat, nClusters = nClusters,
    nObserved = nObserved,
    degFree = degFree
  )
  
  out
}

estimateCellCounts2_GMset<-function (GenomicMethylSet,referenceset = NULL, CustomCpGs = NULL, returnAll = FALSE, 
    meanPlot = FALSE, verbose = TRUE, lessThanOne = FALSE, cellcounts = NULL, 
    ...) 
{
    ###parameter
    probeSelect = "IDOL"
    cellTypes = c("CD8T", "CD4T", "NK","Bcell", "Mono", "Neu")
    referencePlatform = "IlluminaHumanMethylationEPIC"
    processMethod = "preprocessQuantile"
    compositeCellType = "Blood"
    rgPlatform <- sub("IlluminaHumanMethylation", "", annotation(GenomicMethylSet)[which(names(annotation(GenomicMethylSet)) ==                                                                                       "array")])
    platform <- sub("IlluminaHumanMethylation", "", referencePlatform)
    
    referencePkg <- sprintf("FlowSorted.%s.%s", compositeCellType, 
                            platform)
    subverbose <- max(as.integer(verbose) - 1L, 0L)
    
    
    if (!is.null(referenceset)) {
      referenceRGset <- get(referenceset)
      if (is(referenceRGset, "RGChannelSetExtended")) {
        referenceRGset <- as(referenceRGset, "RGChannelSet")
      }
    }
    else {
      if (!require(referencePkg, character.only = TRUE) && 
          referencePkg != "FlowSorted.BloodExtended.EPIC") {
        stop(strwrap(sprintf("Could not find reference data package for\n                                compositeCellType '%s' and referencePlatform\n                                '%s' (inferred package name is '%s')", 
                             compositeCellType, platform, referencePkg), width = 80, 
                     prefix = " ", initial = ""))
      }
      if (!require(referencePkg, character.only = TRUE) && 
          referencePkg == "FlowSorted.BloodExtended.EPIC") {
        stop(strwrap(sprintf("Could not find reference data package for\n                                compositeCellType '%s' and referencePlatform\n                                '%s' (inferred package name is '%s'),\n                                please contact\n                                Technology.Transfer@dartmouth.edu", 
                             compositeCellType, platform, referencePkg), width = 80, 
                     prefix = " ", initial = ""))
      }
      if ((referencePkg != "FlowSorted.Blood.EPIC") && (referencePkg != 
                                                        "FlowSorted.CordBloodCombined.450k")) {
        referenceRGset <- get(referencePkg)
      }
      else if ((referencePkg == "FlowSorted.Blood.EPIC") | 
               (referencePkg == "FlowSorted.CordBloodCombined.450k")) {
        referenceRGset <- libraryDataGet(referencePkg)
      }
    }
    if (rgPlatform != platform) {
      GenomicMethylSet <- convertArray(GenomicMethylSet, outType = referencePlatform, 
            verbose = TRUE)
    }
    if (!base::all(cellTypes %in% referenceRGset$CellType)) {
        stop(strwrap(sprintf("all elements of argument 'cellTypes' needs to be\n                            part of the reference phenoData columns 'CellType'\n                            (containg the following elements: '%s')", 
            paste(unique(referenceRGset$cellType), collapse = "', '")), 
            width = 80, prefix = " ", initial = ""))
    }
    if (length(unique(cellTypes)) < 2) {
        stop("At least 2 cell types must be provided.")
    }
    
    if ((probeSelect == "IDOL") && (compositeCellType == "Blood")) {
        if ((rgPlatform == "450k")) {
            CustomCpGs <- FlowSorted.Blood.EPIC::IDOLOptimizedCpGs450klegacy
        }
        else {
            CustomCpGs <- FlowSorted.Blood.EPIC::IDOLOptimizedCpGs
        }
    }
    if (verbose) {
        message(strwrap("[estimateCellCounts2_GMset] Combining user data with\n                        reference (flow sorted) data.\n", 
            width = 80, prefix = " ", initial = ""))
    }
    newpd <- DataFrame(sampleNames = base::c(colnames(GenomicMethylSet), colnames(referenceRGset)), 
        studyIndex = rep(base::c("user", "reference"), times = base::c(ncol(GenomicMethylSet), 
            ncol(referenceRGset))))
    referenceRGset$CellType <- as.character(referenceRGset$CellType)
    if (is.null(GenomicMethylSet$CellType)) {
      GenomicMethylSet$CellType <- rep("NA", dim(GenomicMethylSet)[2])
    }
    if (is.null(GenomicMethylSet$Age)) {
      GenomicMethylSet$Age <- rep("NA", dim(GenomicMethylSet)[2])
    }
    if (is.null(GenomicMethylSet$Sex)) {
      GenomicMethylSet$Sex <- rep("NA", dim(GenomicMethylSet)[2])
    }
    if (is.null(referenceRGset$Sex)) {
        referenceRGset$Sex <- rep("NA", dim(referenceRGset)[2])
    }else {
        referenceRGset$Sex <- as.character(referenceRGset$Sex)
    }
    if (is.null(referenceRGset$Age)) {
        referenceRGset$Age <- rep("NA", dim(referenceRGset)[2])
    }else {
        try(referenceRGset$Age <- as.numeric(referenceRGset$Age))
    }
    commoncolumn <- intersect(names(colData(GenomicMethylSet)), names(colData(referenceRGset)))
    restry <- try({
        colData(GenomicMethylSet)[commoncolumn] <- mapply(FUN = as, colData(GenomicMethylSet)[commoncolumn], 
            vapply(colData(referenceRGset)[commoncolumn], class, 
                FUN.VALUE = character(1)), SIMPLIFY = FALSE)
    }, silent = TRUE)
    if ("try-error" %in% class(restry)) {
        commoncolumn <- base::c("CellType", "Sex", "Age")
        colData(GenomicMethylSet)[commoncolumn] <- mapply(FUN = as, colData(GenomicMethylSet)[commoncolumn], 
            vapply(colData(referenceRGset)[commoncolumn], class, 
                FUN.VALUE = character(1)), SIMPLIFY = FALSE)
    }else {
        colData(GenomicMethylSet)[commoncolumn] <- mapply(FUN = as, colData(GenomicMethylSet)[commoncolumn], 
            vapply(colData(referenceRGset)[commoncolumn], class, 
                FUN.VALUE = character(1)), SIMPLIFY = FALSE)
    }
    colData(referenceRGset) <- colData(referenceRGset)[commoncolumn]
    colData(GenomicMethylSet) <- colData(GenomicMethylSet)[commoncolumn]
    referencePd <- colData(referenceRGset)

    referenceRGset<-preprocessRaw(referenceRGset)
    referenceGMset<-mapToGenome(referenceRGset)
    reference_sex<-getSex(referenceGMset);rm(referenceRGset) #####predict Sex
    ##合并sex信息#
    reference_sex<- as.data.frame(reference_sex);
    reference_sex %<>% dplyr::select(predictedSex) %>% dplyr::rename(Sex=predictedSex)%>% dplyr::mutate(samplename=rownames(reference_sex))
    GM_sex<-as.data.frame(colData(GenomicMethylSet));
    GM_sex %<>% dplyr::select(Sex) %>% dplyr::mutate(samplename=rownames(GM_sex))
    combined_Sex<-rbind(GM_sex,reference_sex)
    
    combinedGMset <- combineArrays(GenomicMethylSet, referenceGMset, outType = referencePlatform)
    colData(combinedGMset) <- newpd
    colnames(combinedGMset) <- newpd$sampleNames
    GenomicMethylSet_colnames<-colnames(GenomicMethylSet)
    GenomicMethylSet_colData<-colData(GenomicMethylSet)
    rm(referenceGMset,GenomicMethylSet)
    combined_Sex<-combined_Sex[match(sampleNames(combinedGMset),rownames(combined_Sex)),]#对齐
    if (verbose) {
        message(strwrap("[estimateCellCounts2_GMset] Processing user and reference\n                        data together.\n", 
            width = 80, prefix = " ", initial = ""))
    }
    combinedGMset<-fixMethOutliers(combinedGMset);gc()
    combinedGMset<-preprocessQuantile(combinedGMset,fixOutliers = F,sex = combined_Sex$Sex);gc()
    
    referenceMset <- combinedGMset[, combinedGMset$studyIndex == 
        "reference"]
    colData(referenceMset) <- as(referencePd, "DataFrame")
    mSet <- combinedGMset[, combinedGMset$studyIndex == "user"]
    colData(mSet) <- as(GenomicMethylSet_colData, "DataFrame")
    rm(combinedGMset)
    if (probeSelect != "IDOL") {
        if (verbose) {
            message(strwrap("[estimateCellCounts2_GMset] Picking probes for\n                            composition estimation.\n", 
                width = 80, prefix = " ", initial = ""))
        }
        compData <- pickCompProbes(referenceMset, cellTypes = cellTypes, 
            compositeCellType = compositeCellType, probeSelect = probeSelect)
        coefs <- compData$coefEsts
        if (verbose) {
            message(strwrap("[estimateCellCounts2_GMset] Estimating  proportion\n                            composition (prop), if you provide cellcounts\n                            those will be provided as counts in the\n                            composition estimation.\n", 
                width = 80, prefix = " ", initial = ""))
        }
        prop <- projectCellType_CP(getBeta(mSet)[rownames(coefs), 
            ], coefs, lessThanOne = lessThanOne)
        prop <- round(prop, 4)
        rownames(prop) <- GenomicMethylSet_colnames;
        counts <- round(prop * cellcounts, 0)
        if (returnAll) {
            list(prop = prop, counts = counts, compTable = compData$compTable, 
                normalizedData = mSet)
        }
        else {
            list(prop = prop, counts = counts)
        }
    }else {
        if (verbose) {
            message(strwrap("[estimateCellCounts2_GMset] Using IDOL L-DMR probes for\n                            composition estimation.\n", 
                width = 80, prefix = " ", initial = ""))
        }
        p <- getBeta(referenceMset)
        pd <- as.data.frame(colData(referenceMset))
        rm(referenceMset)
        if (!is.null(cellTypes)) {
            if (!base::all(cellTypes %in% pd$CellType)) {
                stop(strwrap("elements of argument 'cellTypes' are not part of\n                            'referenceMset$CellType'", 
                  width = 80, prefix = " ", initial = ""))
            }
            keep <- which(pd$CellType %in% cellTypes)
            pd <- pd[keep, ]
            p <- p[, keep]
        }
        pd$CellType <- factor(pd$CellType, levels = cellTypes)
        ffComp <- genefilter::rowFtests(p, pd$CellType)
        tIndexes <- split(seq(along = pd$CellType), pd$CellType)
        prof <- vapply(tIndexes, function(i) rowMeans(p[, i]), 
            FUN.VALUE = numeric(dim(p)[1]))
        r <- rowRanges(p)
        compTable <- cbind(ffComp, prof, r, abs(r[, 1] - r[, 
            2]))
        names(compTable)[1] <- "Fstat"
        names(compTable)[base::c(-2, -1, 0) + ncol(compTable)] <- base::c("low", "high", "range")
        tstatList <- lapply(tIndexes, function(i) {
            x <- rep(0, ncol(p))
            x[i] <- 1
            return(genefilter::rowttests(p, factor(x)))
        })
        trainingProbes <- CustomCpGs
        trainingProbes <- trainingProbes[trainingProbes %in% 
            rownames(p)]
        p <- p[trainingProbes, ]
        pMeans <- colMeans(p)
        names(pMeans) <- pd$CellType
        form <- as.formula(sprintf("y ~ %s - 1", paste(levels(pd$CellType), 
            collapse = "+")))
        phenoDF <- as.data.frame(model.matrix(~pd$CellType - 
            1))
        colnames(phenoDF) <- sub("^pd\\$CellType", "", colnames(phenoDF))
        if (ncol(phenoDF) == 2) {
            X <- as.matrix(phenoDF)
            coefEsts <- t(solve(t(X) %*% X) %*% t(X) %*% t(p))
            coefs <- coefEsts
        }else {
            tmp <- validationCellType(Y = p, pheno = phenoDF, 
                modelFix = form)
            coefEsts <- tmp$coefEsts
            coefs <- coefEsts
        }
        compData <- list(coefEsts = coefEsts, compTable = compTable, 
            sampleMeans = pMeans)
        if (verbose) {
            message(strwrap("[estimateCellCounts2_GMset] Estimating  proportion\n                            composition (prop), if you provide cellcounts\n                            those will be provided as counts in the\n                            composition estimation.\n", 
                width = 80, prefix = " ", initial = ""))
        }
        prop <- projectCellType_CP(getBeta(mSet)[rownames(coefs), 
            ], coefs, lessThanOne = lessThanOne)
        prop <- round(prop, 4)
        rownames(prop) <- GenomicMethylSet_colnames;
        counts <- round(prop * cellcounts, 0)
        if (returnAll) {
            list(prop = prop, counts = counts, compTable = compTable, 
                normalizedData = mSet)
        }
        else {
            list(prop = prop, counts = counts)
        }
    }
}