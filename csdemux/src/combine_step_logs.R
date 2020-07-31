#!/usr/bin/env Rscript

log <- file(snakemake@log[[1]],
            open = "wt")
sink(log, type = "message")
sink(log, type = "output", append = TRUE)

library(data.table)
library(bit64)

# step_files <- list.files("test/stats",
#                          pattern = "trim.txt",
#                          full.names = TRUE)
# 
step_files <- snakemake@input[["step_files"]]


names(step_files) <- sapply(step_files, function(x)
  unlist(strsplit(basename(x), split = ".", fixed = TRUE))[[1]])

FreadIndiv <- function(x) {
  indiv_data <- fread(x, header = FALSE)
  indiv_data[, .(
    type = gsub("[^[:alpha:]]+", "", V1),
    reads = as.numeric(unlist(strsplit(V2, " "))[[1]]),
    bases = as.numeric(unlist(strsplit(V3, " "))[[1]])
  ), by = 1:nrow(indiv_data)][
    , .(type, reads, bases)]
}

step_data <- rbindlist(lapply(step_files, FreadIndiv), idcol = "indiv")
fwrite(step_data, snakemake@output[["step_data"]])

sessionInfo()
