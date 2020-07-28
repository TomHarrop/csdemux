#!/usr/bin/env Rscript

log <- file(snakemake@log[[1]],
            open = "wt")
sink(log, type = "message")
sink(log, type = "output", append = TRUE)

library(data.table)
library(bit64)
library(scales)
library(ggplot2)

# stats_file <- "output/010_demux/stats.txt"
stats_file <- snakemake@input[["stats"]]

demux_stats <- fread(stats_file, 
                     skip = 4)

setorder(demux_stats, -Bases)
demux_stats[, bco := factor(`#Name`, levels = unique(`#Name`))]

# cumulative fraction of barcodes
demux_stats[, cf := cumsum(Reads) / sum(Reads)]
demux_stats[, cs := cumsum(Reads)]

# find the inflexion point
demux_stats[, first_diff := c(NA, diff(Reads))]
demux_stats[, second_diff := c(NA, diff(first_diff))]
first_inflexion <- demux_stats[which.max(second_diff),
            as.integer(bco)]
mylab <- paste(first_inflexion, "libraries detected")

# draw plot
mycols <- viridis::viridis(3)
gp <- ggplot(demux_stats[1:1000], aes(x = as.integer(bco), y = cf, group = 1)) +
    theme_minimal(base_size = 14) + 
    scale_x_continuous(trans = log_trans(base = 4),
                       breaks = trans_breaks(function(x) log(x, 4),
                                             function(x) 4^x)) +
    xlab("Barcode rank") +
    ylab("Cumulative fraction of reads") +
    geom_vline(xintercept = first_inflexion,
               colour = mycols[[2]],
               linetype = 2) +
    geom_text(aes(x = first_inflexion,
                  y = 1,
                  label = mylab),
              colour = mycols[[2]],
              hjust = 1,
              nudge_x = -0.01) +
    geom_path(colour = mycols[[1]])

ggsave(snakemake@output[["plot"]],
       gp,
       width = 210,
       height = 148,
       units = "mm",
       device = cairo_pdf)

sessionInfo()
