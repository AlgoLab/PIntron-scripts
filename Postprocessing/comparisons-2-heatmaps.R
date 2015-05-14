#!/usr/bin/env Rscript

## if needed:
#install.packages(c("tidyr", "dplyr", "reshape2", "pheatmap"))

library(tidyr)
library(dplyr)
library(reshape2)
library(pheatmap)



patients <- c("P152", "P174", "P209")
timepoints <- c("T0", "T12", "T24", "T48", "T96")

load.single.report <- function(report.filename) {
    message("Loading report from ", report.filename)
    df <- read.csv(report.filename, header=FALSE, stringsAsFactors=FALSE)
    colnames(df) <- c("Gene ID", "Gene Name", "Transcript Name", "Patient", "Timepoint", "Value")

    df$Value[df$Value == "True"] <- "1"
    df$Value[df$Value == "False"] <- "0"
    df$Value <- as.integer(df$Value)

    base.df <- expand.grid(`Gene ID`=unique(df$`Gene ID`),
                           `Gene Name`=unique(df$`Gene Name`),
                           `Transcript Name`=unique(df$`Transcript Name`),
                           Patient=patients, Timepoint=timepoints,
                           KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE
                           )

    missing.points <- dplyr::setdiff(base.df, df %>% select(-Value))
    if (nrow(missing.points)>0) {
        missing.points$Value <- 0
        df <- dplyr::bind_rows(df, missing.points)
    }

    df$Dataset <- paste(df$Patient, df$Timepoint)
    df
}

produce.heatmap <- function(df, out.filename) {
    message("Printing heatmap to file ", out.filename)
    df.wide <- (    df
                %>% dplyr::select(-`Gene ID`, -`Gene Name`, -`Dataset`)
                %>% dcast(`Transcript Name` ~ `Timepoint` + `Patient`, value.var="Value")
                )
    rownames(df.wide) <- df.wide$`Transcript Name`
    df.wide$`Transcript Name` <- NULL
    df.wide.mat <- as.matrix(df.wide)
    png(filename=out.filename,
        width=800,
        height=max(400, 40*length(unique(df$`Transcript Name`))))
    pheatmap(df.wide.mat,
             color=c("#AA3939", "#2D882D"),
             cluster_cols=FALSE,
             cluster_rows=length(unique(df$`Transcript Name`))>1,
             gaps_col=(1:(length(timepoints)-1))*length(patients),
             breaks=c(0, 0.5, 1),
             legend_breaks=c(0.25, 0.75),
             legend_labels=c("Missing", "Present"),
             treeheight_row=0,
             display_numbers=T,
             number_format="%d",
             fontsize=14,
             number_color="white", fontsize_number=11,
             main=paste("Gene: ", unique(df$`Gene Name`),
                 " [", unique(df$`Gene ID`), "]", sep=""),
             )
    dev.off()
    invisible(NULL)
}

args <- commandArgs(trailingOnly = TRUE)
for (report.filename in args) {
    df <- load.single.report(report.filename)
    out.filename <- paste(sub("^(.*)\\.[^.]*$", "\\1", report.filename), ".png", sep="")
    produce.heatmap(df, out.filename)
}
