#Script to remove adapters, primers, and padding and linkers from sequencing data.

# Initialisations

## Necessary libraries and software
cutadapt <- "D:/Users/e0453664/AppData/Local/miniconda3/Scripts/cutadapt.exe" # Set to path to cutadapt.
system2(cutadapt, args = "--version")

library(ShortRead); packageVersion("ShortRead")
library(seqTools); packageVersion("seqTools")
library(microseq); packageVersion("microseq")


## Loading the sequences
getwd()
setwd("./") # Set to working directory of where to save all working files.
path <- "./" # CHANGE ME to the directory containing the fastq files after unzipping.
list.files(path)

#Assuming, forward and reverse fastq filenames have format: SAMPLENAME_1.fastq and SAMPLENAME_2.fastq, we extract the sample names.
fnFs <- sort(list.files(path, pattern="_1.fastq.gz", full.names = TRUE))
fnRs <- sort(list.files(path, pattern="_2.fastq.gz", full.names = TRUE))
# Extract sample names, assuming filenames have format: SAMPLENAME_XXX.fastq
sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)

## Initialising functions
# Function to get all orientations of a sequence.
allOrients <- function(sequence) {
  # Create all orientations of the input sequence.
  require(Biostrings)
  dna <- DNAString(sequence)  # The Biostrings works w/ DNAString objects rather than character vectors.
  orients <- c(Forward = dna, Complement = complement(dna), Reverse = reverse(dna), 
               RevComp = reverseComplement(dna))
  return(sapply(orients, toString))  # Convert back to character vector.
}

# Function to count presence of a sequence in reads
sequenceHits <- function(sequence, fn) {
  nhits <- vcountPattern(sequence, sread(readFastq(fn)), fixed = FALSE)
  return(sum(nhits > 0))
}

## Pre-filtering to remove ambiguous bases (Ns)
fnFs.filtN <- file.path(path, "N_filtered", basename(fnFs)) # Put N-filtered files in filtN/ subdirectory.
fnRs.filtN <- file.path(path, "N_filtered", basename(fnRs))
filterAndTrim(fnFs, fnFs.filtN, fnRs, fnRs.filtN, maxN = 0, multithread = TRUE) # Set multithread=FALSE on windows.

# Adapter removal
# Declare adapter sequence
FWD_adapter <- "TCGTCGGCAGCGTCAGATGTGTATAAGAGACAG"
REV_adapter <- "GTCTCGTGGGCTCGGAGATGTGTATAAGAGACAG"

# Getting all orientations
FWD_adapter.orients <- allOrients(FWD_adapter)
REV_adapter.orients <- allOrients(REV_adapter)
FWD_adapter.orients

# Counting adapters containing reads before removal
rbind(FWD_adapter.ForwardReads = sapply(FWD_adapter.orients, sequenceHits, fn = fnFs.filtN[[1]]), 
      FWD_adapter.ReverseReads = sapply(FWD_adapter.orients, sequenceHits, fn = fnRs.filtN[[1]]), 
      REV_adapter.ForwardReads = sapply(REV_adapter.orients, sequenceHits, fn = fnFs.filtN[[1]]), 
      REV_adapter.ReverseReads = sapply(REV_adapter.orients, sequenceHits, fn = fnRs.filtN[[1]]))

# Removing adapters
path.cut_adapter <- file.path(path, "adapter_Removed")
if(!dir.exists(path.cut_adapter)) dir.create(path.cut_adapter)
fnFs.cut_adapter <- file.path(path.cut_adapter, basename(fnFs))
fnRs.cut_adapter <- file.path(path.cut_adapter, basename(fnRs))

FWD_adapter.RC <- dada2:::rc(FWD_adapter)
REV_adapter.RC <- dada2:::rc(REV_adapter)
# Trim FWD and the reverse-complement of REV off of R1 (forward reads)
R1.flags <- paste("-g", FWD_adapter, "-a", REV_adapter.RC) 
# Trim REV and the reverse-complement of FWD off of R2 (reverse reads)
R2.flags <- paste("-G", REV_adapter, "-A", FWD_adapter.RC) 

# Run cutadapt
outputStatsadapter <- capture.output(
  for(i in seq_along(fnFs)) {
    system2(cutadapt, args = c(R1.flags, R2.flags, "-n", 2, # -n 2 required to remove FWD and REV from reads
                               "-m", 1,
                               "-o", fnFs.cut_adapter[i], "-p", fnRs.cut_adapter[i], # output files
                               fnFs.filtN[i], fnRs.filtN[i])) # input files
  }
)
cat(outputStatsadapter, file="./Outputs/cutadapt_adapter_trimming_stats.txt", sep="\n", append = FALSE)

# Counting adapter containing reads after removal
rbind(FWD_adapter.ForwardReads = sapply(FWD_adapter.orients, sequenceHits, fn = fnFs.cut_adapter[[1]]), 
      FWD_adapter.ReverseReads = sapply(FWD_adapter.orients, sequenceHits, fn = fnRs.cut_adapter[[1]]), 
      REV_adapter.ForwardReads = sapply(REV_adapter.orients, sequenceHits, fn = fnFs.cut_adapter[[1]]), 
      REV_adapter.ReverseReads = sapply(REV_adapter.orients, sequenceHits, fn = fnRs.cut_adapter[[1]]))

# Primer removal
# Declare primer sequence
FWD_primer <- "GTGYCAGCMGCCGCGGTAA"  
REV_primer <- "CCGYCAATTYMTTTRAGTTT"

# Getting all orientations
FWD_primer.orients <- allOrients(FWD_primer)
REV_primer.orients <- allOrients(REV_primer)
FWD_primer.orients

# Counting primer containing reads before removal
rbind(FWD_primer.ForwardReads = sapply(FWD_primer.orients, sequenceHits, fn = fnFs.cut_adapter[[1]]), 
      FWD_primer.ReverseReads = sapply(FWD_primer.orients, sequenceHits, fn = fnRs.cut_adapter[[1]]), 
      REV_primer.ForwardReads = sapply(REV_primer.orients, sequenceHits, fn = fnFs.cut_adapter[[1]]), 
      REV_primer.ReverseReads = sapply(REV_primer.orients, sequenceHits, fn = fnRs.cut_adapter[[1]]))

# Removing primers
path.cut_primer <- file.path(path, "Primer_adapter_Removed")
if(!dir.exists(path.cut_primer)) dir.create(path.cut_primer)
fnFs.cut_primer <- file.path(path.cut_primer, basename(fnFs))
fnRs.cut_primer <- file.path(path.cut_primer, basename(fnRs))

FWD_primer.RC <- dada2:::rc(FWD_primer)
REV_primer.RC <- dada2:::rc(REV_primer)
# Trim FWD and the reverse-complement of REV off of R1 (forward reads)
R1.flags <- paste("-g", FWD_primer, "-a", REV_primer.RC) 
# Trim REV and the reverse-complement of FWD off of R2 (reverse reads)
R2.flags <- paste("-G", REV_primer, "-A", FWD_primer.RC) 

# Run cutadapt
outputStatsPrimer <- capture.output(
  for(i in seq_along(fnFs)) {
    system2(cutadapt, args = c(R1.flags, R2.flags, "-n", 2, # -n 2 required to remove FWD and REV from reads
                               "-m", 1,
                               "-o", fnFs.cut_primer[i], "-p", fnRs.cut_primer[i], # output files
                               fnFs.cut_adapter[i], fnRs.cut_adapter[i])) # input files
  }
)
cat(outputStatsPrimer, file="./Outputs/cutadapt_primer_trimming_stats.txt", sep="\n", append = FALSE)

# Counting primer containing reads after removal
rbind(FWD_primer.ForwardReads = sapply(FWD_primer.orients, sequenceHits, fn = fnFs.cut_primer[[1]]), 
      FWD_primer.ReverseReads = sapply(FWD_primer.orients, sequenceHits, fn = fnRs.cut_primer[[1]]), 
      REV_primer.ForwardReads = sapply(REV_primer.orients, sequenceHits, fn = fnFs.cut_primer[[1]]), 
      REV_primer.ReverseReads = sapply(REV_primer.orients, sequenceHits, fn = fnRs.cut_primer[[1]]))

# Padding and linker removal
# Declare padding + linker sequence
FWD_pl <- "AATGATACGGCGACCACCGAGATCTACAC"
REV_pl <- "CAAGCAGAAGACGGCATACGAGAT"

# Getting all orientations
FWD_pl.orients <- allOrients(FWD_pl)
REV_pl.orients <- allOrients(REV_pl)
FWD_pl.orients

# Counting padding + linker containing reads before removal
rbind(FWD_pl.ForwardReads = sapply(FWD_pl.orients, sequenceHits, fn = fnFs.cut_primer[[1]]), 
      FWD_pl.ReverseReads = sapply(FWD_pl.orients, sequenceHits, fn = fnRs.cut_primer[[1]]), 
      REV_pl.ForwardReads = sapply(REV_pl.orients, sequenceHits, fn = fnFs.cut_primer[[1]]), 
      REV_pl.ReverseReads = sapply(REV_pl.orients, sequenceHits, fn = fnRs.cut_primer[[1]]))

# Removing padding and linker
path.cut_pl <- file.path(path, "Final_Processed")
if(!dir.exists(path.cut_pl)) dir.create(path.cut_pl)
fnFs.cut_pl <- file.path(path.cut_pl, basename(fnFs))
fnRs.cut_pl <- file.path(path.cut_pl, basename(fnRs))

FWD_pl.RC <- dada2:::rc(FWD_pl)
REV_pl.RC <- dada2:::rc(REV_pl)
# Trim FWD and the reverse-complement of REV off of R1 (forward reads)
R1.flags <- paste("-g", FWD_pl, "-a", REV_pl.RC) 
# Trim REV and the reverse-complement of FWD off of R2 (reverse reads)
R2.flags <- paste("-G", REV_pl, "-A", FWD_pl.RC) 
# Run Cutadapt
system2(cutadapt, args = c(R1.flags, R2.flags, "-n", 2, # -n 2 required to remove FWD and REV from reads
                           "-m", 1,
                           "-o", fnFs.cut_pl[5], "-p", fnRs.cut_pl[5], # output files
                           fnFs.cut_primer[5], fnRs.cut_primer[5])) # input files
outputStatsPL <- capture.output(
  for(i in seq_along(fnFs)) {
    system2(cutadapt, args = c(R1.flags, R2.flags, "-n", 2, # -n 2 required to remove FWD and REV from reads
                               "-m", 1,
                               "-o", fnFs.cut_pl[i], "-p", fnRs.cut_pl[i], # output files
                               fnFs.cut_primer[i], fnRs.cut_primer[i])) # input files
  }
)
cat(outputStatsPL, file="./Outputs/cutadapt_pl_trimming_stats.txt", sep="\n", append = FALSE)

# Counting padding and linker containing reads after removal
rbind(FWD_pl.ForwardReads = sapply(FWD_pl.orients, sequenceHits, fn = fnFs.cut_pl[[1]]), 
      FWD_pl.ReverseReads = sapply(FWD_pl.orients, sequenceHits, fn = fnRs.cut_pl[[1]]), 
      REV_pl.ForwardReads = sapply(REV_pl.orients, sequenceHits, fn = fnFs.cut_pl[[1]]), 
      REV_pl.ReverseReads = sapply(REV_pl.orients, sequenceHits, fn = fnRs.cut_pl[[1]]))

# Summary calculations
## Getting all filenames
# Forward and reverse fastq filenames have the format SAMPLENAME_X.fastq.gz
cutFs_raw <- sort(list.files(path, pattern = "_1.fastq.gz", full.names = TRUE))
cutRs_raw <- sort(list.files(path, pattern = "_2.fastq.gz", full.names = TRUE))
cut_raw <- sort(list.files(path, pattern = "fastq.gz", full.names = TRUE))

path.filtN <- file.path(path, "N_filtered")
cutFs_filtN <- sort(list.files(path.filtN, pattern = "_1.fastq.gz", full.names = TRUE))
cutRs_filtN <- sort(list.files(path.filtN, pattern = "_2.fastq.gz", full.names = TRUE))
cut_filtN <- sort(list.files(path.filtN, pattern = "fastq.gz", full.names = TRUE))

cutFs_adapter <- sort(list.files(path.cut_adapter, pattern = "_1.fastq.gz", full.names = TRUE))
cutRs_adapter <- sort(list.files(path.cut_adapter, pattern = "_2.fastq.gz", full.names = TRUE))
cut_adapter <- sort(list.files(path.cut_adapter, pattern = "fastq.gz", full.names = TRUE))

cutFs_primer <- sort(list.files(path.cut_primer, pattern = "_1.fastq.gz", full.names = TRUE))
cutRs_primer <- sort(list.files(path.cut_primer, pattern = "_2.fastq.gz", full.names = TRUE))
cut_primer <- sort(list.files(path.cut_primer, pattern = "fastq.gz", full.names = TRUE))

cutFs_pl <- sort(list.files(path.cut_pl, pattern = "_1.fastq.gz", full.names = TRUE))
cutRs_pl <- sort(list.files(path.cut_pl, pattern = "_2.fastq.gz", full.names = TRUE))
cut_pl <- sort(list.files(path.cut_pl, pattern = "fastq.gz", full.names = TRUE))

# Extract sample names, assuming filenames have the format SAMPLENAME_X
get.sample.name <- function(fname) strsplit(basename(fname), "_")[[1]][1]
sample.names <- unname(sapply(cutFs_pl, get.sample.name))
head(sample.names)

samples <- sort(list.files(path, pattern = "fastq.gz", full.names = FALSE))

## Comparing number of reads for each sample during each step
raw_seqs <- numeric(length(cut_raw))
prefiltered_seqs <- numeric(length(cut_raw))
adapterremoved_seqs <- numeric(length(cut_raw))
primerremoved_seqs <- numeric(length(cut_raw))
padlinkerremoved_seqs <- numeric(length(cut_raw))

for (j in 1:length(samples)){
  raw_seqs[j] <- nrow(readFastq(cut_raw[j]))
  prefiltered_seqs[j] <- nrow(readFastq(cut_filtN[j]))
  adapterremoved_seqs[j] <- nrow(readFastq(cut_adapter[j]))
  primerremoved_seqs[j] <- nrow(readFastq(cut_primer[j]))
  padlinkerremoved_seqs[j] <- nrow(readFastq(cut_pl[j]))
}
fractionOfRawDataInFinal_seqs <- padlinkerremoved_seqs / raw_seqs

sample_summaries <- data.frame(samples, raw_seqs, prefiltered_seqs, adapterremoved_seqs, primerremoved_seqs, padlinkerremoved_seqs, fractionOfRawDataInFinal_seqs, stringsAsFactors = FALSE)
write.csv(sample_summaries, file = "./Outputs/Number of Sequence Summary (Adapter removal).csv")
summary(sample_summaries)
hist(sample_summaries$raw_seqs)
hist(sample_summaries$padlinkerremoved_seqs)

## Comparing file sizes for each step
raw_size <- numeric(length(cut_raw))
prefiltered_size <- numeric(length(cut_filtN))
adaptorremoved_size <- numeric(length(cut_adaptor))
primerremoved_size <- numeric(length(cut_primer))
padlinkerremoved_size <- numeric(length(cut_pl))

for (j in 1:length(samples)){
  raw_size[j] <- file.size(cut_raw[j])
  prefiltered_size[j] <- file.size(cut_filtN[j])
  adaptorremoved_size[j] <- file.size(cut_adaptor[j])
  primerremoved_size[j] <- file.size(cut_primer[j])
  padlinkerremoved_size[j] <- file.size(cut_pl[j])
}
fractionOfRawDataInFinal_size <- padlinkerremoved_size / raw_size

file_info <- data.frame(samples, raw_size, prefiltered_size, adaptorremoved_size, primerremoved_size, padlinkerremoved_size, fractionOfRawDataInFinal_size, stringsAsFactors = FALSE)
write.csv(file_info, file = "./Outputs/Size of Files Summary (Adapter removal).csv")
summary(file_info)
hist(file_info$raw_size)
hist(sample_summaries$padlinkerremoved_size)
