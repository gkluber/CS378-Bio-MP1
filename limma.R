# Title     : COVID-19 LIMMA
# Objective : Run LIMMA on Cite-Seq COVID-19 data
# Created by: gklub
# Created on: 10/8/2020

#if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")

#BiocManager::install("edgeR")

library(edgeR)
library(readr)
library(data.table)

patients <- read.csv(file = 'SRARunTable.csv')
#head(patients)

genecounts <- NULL
for (run in patients[,1])
{
  #writeLines(paste(run, "\n"))
  filename <- paste0("count/", run, ".tsv")
  print(filename)
  #count <- read.table(file = filename, sep = '\t', header = TRUE)
  count <- read_tsv(filename, col_names = FALSE)
  print(head(count))
  #args <- vector(mode="list", length = nrow(count))
  #targs <- list()
  #for (row in seq_len(nrow(count)))
  #{
  #  targs[[count[row, 1]]] <- count[row, 2]
  #}
  #tpose <- do.call(data.frame, targs)
  tpose <- transpose(count)
  #print(tpose[1:20])
  colnames(tpose) <- tpose[1,]
  tpose <- tpose[2,]
  rownames(tpose) <- as.list(run)
  print(tpose[1:20])
  if (is.null(genecounts))
  {
    genecounts <<- tpose
  }
  else
  {
    genecounts <<- rbind(genecounts, tpose)
  }
}

genes <- transpose(genecounts)
rownames(genes) = colnames(genecounts)
colnames(genes) = rownames(genecounts)
head(genes)
