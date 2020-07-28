#!/usr/bin/env Rscript

log <- file(snakemake@log[[1]],
            open = "wt")
sink(log, type = "message")
sink(log, type = "output", append = TRUE)

library(data.table)
library(bit64)
library(ggplot2)
library(scales)

BcLookup <- function(x){
    expected_barcodes[expected_barcode == x, unique(sample)[[1]]]
}

GetMmDt <- function(i){
    my_mm <- merge(barcode_hammings[dist == i],
                   abundant_reads,
                   by.x = "found_barcode",
                   by.y = "barcode")
    my_mm[, sample := BcLookup(expected_barcode),
          by = expected_barcode]
    mydup <- duplicated(my_mm, by = "found_barcode")
    if(any(mydup)){
        dupbc <- my_mm[mydup, unique(found_barcode)]
        print(my_mm[found_barcode == dupbc[[1]]])
        stop(
            paste0("Duplicated barcodes at i==", i,
                   "\nOnly the first duplicate is printed")
        )
    }
    return(my_mm)}

# barcode_hamming_file <- "output/010_demux/hamming_distances.Rds"
# stats_file <- "output/010_demux/stats.txt"
# barcodes_file <- "data/samples.csv"

barcode_hamming_file <- snakemake@input[["foundbc"]] 
barcodes_file <- snakemake@input[["barcodes"]]
stats_file <- snakemake@input[["stats"]]

# found barcodes
demux_stats <- fread(stats_file, 
                     skip = 4)[, .(barcode = `#Name`,
                                   Reads,
                                   Bases)]
setorder(demux_stats, -Reads)
total_reads <- demux_stats[, sum(Reads)]

# only look at reads that were abundant enough to care about
abundant_reads <- demux_stats[Reads / total_reads > 1e-4]

# expected barcodes
expected_barcodes <- fread(barcodes_file)[, .(sample,
                                              expected_barcode = barcode)]

# calculated hamming distances
barcode_hammings <- readRDS(barcode_hamming_file)

# the reads that matched an expected barcode perfectly
no_errors <- merge(expected_barcodes,
                   abundant_reads,
                   by.x = "expected_barcode",
                   by.y = "barcode",
                   all.x = TRUE,
                   all.y = FALSE)
no_errors[, dist := 0]
no_errors[, found_barcode := expected_barcode]

# reads that could be matched with dist < i
dist_errs <- rbindlist(lapply(1:2, GetMmDt))

# mung for plotting
countdt <- rbind(no_errors, dist_errs)
plotdt <- countdt[, .(sample_reads = sum(Reads, na.rm = TRUE)), by = .(sample, dist)]
plotdt[, readfrac := sample_reads / total_reads]
setorder(plotdt, -sample_reads)
sample_order <- plotdt[dist == 0, unique(sample)]

# unassigned reads
unassigned <- demux_stats[, sum(Reads)] - plotdt[, sum(sample_reads)]
plotdt2 <- rbind(plotdt,
      data.table(sample = "Unassigned",
                 dist = NA,
                 sample_reads = unassigned,
                 readfrac = NA))

plotdt2[, sample := factor(sample, levels = c(sample_order, "Unassigned"))]

gp <- ggplot(plotdt2, aes(y = sample_reads/1e6, x = sample, fill = as.factor(dist))) +
    theme_minimal(base_size = 6) + 
    coord_flip() + 
    xlab(NULL) + ylab("Reads (M)") +
    scale_y_continuous(expand = c(0, 0)) +
    scale_fill_viridis_d(na.value = 4,
                         guide = guide_legend(title = "Barcode distance")) +
    geom_col()

ggsave(snakemake@output[["plot"]],
       gp,
       width = 148,
       height = 210,
       units = "mm",
       device = cairo_pdf)

fwrite(countdt, snakemake@output[["report"]])

sessionInfo()
