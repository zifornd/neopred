#!/usr/bin/env Rscript

library(ggplot2)
library(reshape2)
library(dplyr)
library(RColorBrewer)
library(ComplexHeatmap)
library(optparse)

option_list = list(
make_option(c("-i", "--hla"), type="character", default=NULL,
            help="infiltration data from timer or cibersort", metavar="character"),
make_option(c("-m", "--meta"), type="character", default=NULL,
            help="meta information", metavar="character"),
make_option(c("-e", "--expression"), type="character", default=NULL,
            help="tpm expression data", metavar="character"),
make_option(c("-d", "--design"), type="character", default=NULL,
            help="specifying meta groups", metavar="character"),
make_option(c("-o", "--outdir"), type="character", default=NULL,
            help="output directory", metavar="character"),
make_option(c("-p", "--patID"), type="character", default=NULL,
            help="patient ID", metavar="character"),
make_option(c("-s", "--source"), type="character", default=NULL,
            help="hlaoncoplot R script", metavar="character")

);

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

hla <- read.table(opt$hla, header = TRUE, sep = "\t")
meta <- read.table(opt$meta,sep = ",", header = TRUE, row.names = 1, check.names = FALSE, stringsAsFactors = F)
expr <- read.table(opt$expression, header = TRUE, row.names = 1, sep = ",", check.names = FALSE)

overlapped_sam <- intersect(hla$subject, colnames(expr))

hla <- hla[hla$subject %in% overlapped_sam,]
meta <- meta[rownames(meta) %in% overlapped_sam,]
expr <- expr[overlapped_sam]

source(opt$source)
hla_plot <- hla_oncoplot(hla, expr, meta, opt$design, opt$patID)
pdf(file = paste(opt$outdir,"hla_frequency_plot.pdf",sep = ""), width=8, height=6)
print(hla_plot)
dev.off()
