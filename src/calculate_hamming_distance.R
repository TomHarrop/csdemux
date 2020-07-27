#!/usr/bin/env Rscript

log <- file(snakemake@log[[1]],
            open = "wt")
sink(log, type = "message")
sink(log, type = "output", append = TRUE)

library(data.table)
library(bit64)
library(future.apply)

# from https://johanndejong.wordpress.com/2015/10/02/faster-hamming-distance-in-r-2/
hamming <- function(X, Y = NULL) {
    if (is.null(Y)) {
        uniqs <- unique(as.vector(X))
        U <- X == uniqs[1]
        H <- t(U) %*% U
        for ( uniq in uniqs[-1] ) {
            U <- X == uniq
            H <- H + t(U) %*% U
        }
    } else {
        uniqs <- union(X, Y)
        H <- t(X == uniqs[1]) %*% (Y == uniqs[1])
        for ( uniq in uniqs[-1] ) {
            H <- H + t(X == uniq) %*% (Y == uniq)
        }
    }
    nrow(X) - H
}

BcHam <- function(x, y){
    hamming(as.matrix(unlist(strsplit(x, split = ""))),
            as.matrix(unlist(strsplit(y, split = ""))))[1,1]
}

GetBcByDist <- function(bc, exp_bc) {
    data.table(
        found_barcode = bc,
        expected_barcode = exp_bc,
        dist = unlist(lapply(exp_bc, BcHam, y = bc)))}



# stats_file <- "output/010_demux/stats.txt"
# barcodes_file <- "data/combined_sampleinfo.csv"
# threads <- 8

threads <- snakemake@threads[[1]]
barcodes_file <- snakemake@input[["barcodes"]]
stats_file <- snakemake@input[["stats"]]

plan(multiprocess(workers = threads))    
options(future.globals.maxSize = +Inf)

# found barcodes
demux_stats <- fread(stats_file, 
                     skip = 4)[, .(barcode = `#Name`,
                                   Reads,
                                   Bases)]
setorder(demux_stats, -Reads)

# expected barcodes
expected_barcodes <- fread(barcodes_file)[, .(sample, barcode)]

# which found barcodes have 1 error vs. the expected barcodes
all_exp_bc <- expected_barcodes[, unique(barcode)]
barcodes_to_search <- demux_stats[!barcode %in% all_exp_bc, unique(barcode)]
foundbclist <- future_lapply(X = barcodes_to_search,
                             FUN = GetBcByDist,
                             exp_bc = all_exp_bc,
                             future.packages = "data.table")
foundbc <- rbindlist(foundbclist)
saveRDS(foundbc, snakemake@output[["foundbc"]])

# log
sessionInfo()
