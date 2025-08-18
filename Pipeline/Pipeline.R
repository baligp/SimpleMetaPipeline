####################################
# CUTADAPT > DADA2 > IDTAXA pipeline (Simplified version)

# Set CRAN mirror
r <- getOption("repos")
r["CRAN"] <- "http://cran.us.r-project.org"
options(repos = r)

# Set seed
set.seed(0.1)

# Install lulu (still loaded but not used)
devtools::install_github("tobiasgf/lulu")

# Load required packages
library(tidyverse)
library(gridExtra)
library(seqinr)
library(ShortRead)
library(Biostrings)
library(dada2)
library(lulu)
library(DECIPHER)

# Load pipeline functions
function_files <- list.files(file.path(path, "BioinformaticPipeline", "Pipeline", "Functions"))
sapply(file.path(path, "BioinformaticPipeline", "Pipeline", "Functions", function_files), source)

# Determine starting point
StartingStep <- IdentifyLastInputPresent(c("CutadaptOutput", "DadaOutput", "IdtaxaOutput", "BlastOutput"))

# Set default data path
if (is.null(dataPath)) { dataPath <- path }

# 1. Cutadapt
if (StartingStep <= 1) {
    CutadaptOutput <- RunCutadapt(
        FWD = FWD,
        REV = REV,
        multithread = TRUE,
        dataname = dataname,
        UseCutadapt = UseCutadapt
    )
    CacheOutput(CutadaptOutput)
}

# 2. DADA2
if (StartingStep <= 2) {
    DadaOutput <- RunDADA2(
        truncLen = truncLen,
        trimLeft = trimLeft,
        maxN = maxN,
        maxEE = maxEE,
        truncQ = truncQ,
        DesiredSequenceLengthRange = DesiredSequenceLengthRange,
        dataname = dataname,
        multithread = multithread,
        pool = pool,
        MixedOrientation = MixedOrientation,
        NumberOfRuns = NumberOfRuns,
        ReadType = ReadType
    )
    CacheOutput(DadaOutput)
}

# 3. IDTAXA
if (StartingStep <= 3) {
    ConfirmInputsPresent("DadaOutput")

    IdtaxaOutput <- RunIdtaxa(
        IDTAXA = IDTAXA,
        trainingSet = trainingSet,
        TableToMergeTo = DadaOutput$SeqDataTable,
        SeqsToAssign = SeqsToAssign,
        threshold = threshold,
        parallel = parallel,
        multithread = multithread
    )
    CacheOutput(IdtaxaOutput)
}

# 4. BLAST (optional)
if (StartingStep <= 4) {
    ConfirmInputsPresent("CutadaptOutput")
    ConfirmInputsPresent("DadaOutput")
    ConfirmInputsPresent("IdtaxaOutput")

    if (Blast) {
        BlastOutput <- RunBLAST(
            dbname = dbname,
            clustering = "ESV",
            TableToMergeTo = IdtaxaOutput$SeqDataTable,
            assignmentThresholds = assignmentThresholds,
            Blastdesiredranks = Blastdesiredranks,
            multithread = multithread
        )
    } else {
        BlastOutput <- NULL
    }
    CacheOutput(BlastOutput)
}

# 5. Export Results
ConfirmInputsPresent("CutadaptOutput")
ConfirmInputsPresent("DadaOutput")
ConfirmInputsPresent("IdtaxaOutput")
ConfirmInputsPresent("BlastOutput")

# Save DADA2 diagnostics
glist <- lapply(DadaOutput$SecondaryOutputs$DadaPlots, ggplotGrob)
ggsave(file.path(path, "Results", paste0(dataname, "_DadaPlots.pdf")), marrangeGrob(glist, nrow = 1, ncol = 1))

write.csv(CutadaptOutput$PrimerCountSamplesList,
          file = file.path(path, "Results", paste0(dataname, "_PrimerCounts.csv")))
write.csv(CutadaptOutput$PrimerCountCutadaptedSamplesList,
          file = file.path(path, "Results", paste0(dataname, "_PrimerCountsAfterCutadapt.csv")))
write.csv(DadaOutput$SecondaryOutputs$DadaTables,
          file = file.path(path, "Results", paste0(dataname, "_DadaTable.csv")))
write.csv(DadaOutput$SecondaryOutputs$SeqLengthDist,
          file = file.path(path, "Results", paste0(dataname, "_DadaSeqLengthDistribution.csv")))

# Save final taxonomic assignments
if (!is.null(BlastOutput)) {
    saveRDS(BlastOutput$SeqDataTable, file = file.path(path, "Results", paste0(dataname, "_SeqDataTable.RDS")))
    write.csv(BlastOutput$SeqDataTable, file = file.path(path, "Results", paste0(dataname, "_SeqDataTable.csv")))
    write.csv(BlastOutput$FullBlastOutput, file = file.path(path, "Results", paste0(dataname, "_FullBlastOutput.csv")))
    pdf(file = file.path(path, "Results", paste0(dataname, "_TaxaAssignment.pdf")))
    plot(IdtaxaOutput$PlotData)
    dev.off()
} else {
    saveRDS(IdtaxaOutput$SeqDataTable, file = file.path(path, "Results", paste0(dataname, "_SeqDataTable.RDS")))
    write.csv(IdtaxaOutput$SeqDataTable, file = file.path(path, "Results", paste0(dataname, "_SeqDataTable.csv")))
    pdf(file = file.path(path, "Results", paste0(dataname, "_TaxaAssignment.pdf")))
    plot(IdtaxaOutput$PlotData)
    dev.off()
}
