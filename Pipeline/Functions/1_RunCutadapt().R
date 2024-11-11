RunCutadapt<-function(multithread, FWD, REV, dataname=NULL, UseCutadapt=FALSE) {
    
    if (UseCutadapt == FALSE) return(NULL)
    else if (UseCutadapt == TRUE){
        # For single-end reads, only look for R1 files
        fnFs <- sort(list.files(file.path(dataPath, "FASTQs", dataname), 
                              pattern="_R1_001.fastq", 
                              full.names = TRUE,
                              recursive = TRUE))
        
        if(length(fnFs) == 0) {
            stop("No fastq files found. Please check your dataPath and file naming pattern.")
        }
        
        FWD.orients <- allOrients(FWD)
        REV.orients <- allOrients(REV)

        # Extract sample names
        SampleNames <- sapply(strsplit(basename(fnFs), "_R1_001.fastq"), `[`, 1)
        
        # Create paths for filtered outputs (only forward reads)
        filtFs <- createOutputFilePaths(suffix="_F_filt.fastq.gz", 
                                      outputDirectoryPrefix="_prefilteredsequences")

        # Pre-filter only forward reads
        out <- dada2::filterAndTrim(fnFs, filtFs, maxN=0, rm.phix=TRUE, 
                                  compress=TRUE, multithread=multithread, verbose=TRUE)

        # Count primers in forward reads only
        PrimerCountSamplesList <- list()
        for (i in seq_along(filtFs)) {
            PrimerCountSamplesList[[i]] <- rbind(
                FWD.ForwardReads = sapply(FWD.orients, primerHits, fn = filtFs[[i]]),
                REV.ForwardReads = sapply(REV.orients, primerHits, fn = filtFs[[i]]))
        }

        # Setup cutadapt output paths
        path.cut <- file.path(path, "IntermediateOutputs", paste0(dataname,"_CutadaptedSeqs"))
        if(!dir.exists(path.cut)) dir.create(path.cut)
        fnFs.cut <- file.path(path.cut, basename(fnFs))

        FWD.RC <- dada2:::rc(FWD)
        REV.RC <- dada2:::rc(REV)

        # For single-end reads, look for both orientations of both primers
        SE.flags <- paste("-g", FWD,  # Forward primer at start
                         "-a", REV,    # Reverse primer at end
                         "-g", REV.RC, # Reverse-complement of reverse primer at start
                         "-a", FWD.RC) # Reverse-complement of forward primer at end

        # Run Cutadapt
        for(i in seq_along(fnFs)) {
            system2("cutadapt", args = c(
                SE.flags,
                "-n", 2,
                "-o", fnFs.cut[i],
                filtFs[i]))
        }

        # Check remaining primers in cutadapt output
        PrimerCountCutadaptedSamplesList <- list()
        for (i in seq_along(filtFs)) {
            PrimerCountCutadaptedSamplesList[[i]] <- rbind(
                FWD.ForwardReads = sapply(FWD.orients, primerHits, fn = fnFs.cut[[i]]),
                REV.ForwardReads = sapply(REV.orients, primerHits, fn = fnFs.cut[[i]]))
        }

        return(list("PrimerCountSamplesList"=PrimerCountSamplesList, 
                   "PrimerCountCutadaptedSamplesList"=PrimerCountCutadaptedSamplesList))
    }
}
