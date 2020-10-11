# Title     : COVID-19 LIMMA
# Objective : Run LIMMA on Cite-Seq COVID-19 data
# Created by: gklub
# Created on: 10/8/2020

#if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")

#BiocManager::install("edgeR")

library(edgeR)
library(readr)
library(tidyr)

patients <- read.table(file = 'SraRunTable.csv', sep=',', header = TRUE)

#print(head(patients))

genes <- NULL
for (run in patients[,1])
{
  writeLines(paste(run, "\n"))
  filename <- paste0("count/", run, ".tsv")
  print(filename)
  count <- read.table(file = filename, sep = '\t', row.names = 1, header = TRUE)
  #count <- read.tsv(filename, col_names = FALSE)
  #print(head(count))
  #args <- vector(mode="list", length = nrow(count))
  #targs <- list()
  #for (row in seq_len(nrow(count)))
  #{
  #  targs[[count[row, 1]]] <- count[row, 2]
  #}
  #tpose <- do.call(data.frame, targs)
  #print(paste("dimensions = (", paste(dim(tpose), collapse = ","), ")"))
  #print(tpose[1:20])

  colnames(count) <- as.list(run)
  
  #print(tpose[1:20])
  
  
  if (is.null(genes))
  {
    genes <<- count
  }
  else
  {
    genes <<- cbind(genes, count)
  }
}
print(head(genes))

# Preprocessing for LIMMA
d0 <- DGEList(data.matrix(genes))
d0 <- calcNormFactors(d0)

# Filter out low-expressed genes
cutoff <- 1
drop <- which(apply(cpm(d0), 1, max) < cutoff)
d <- d0[-drop,]

# Define groups
group <- factor(patients$TREATMENT)
png(file="mds.png")
plotMDS(d, col = as.numeric(group))
dev.off()

# Voom
mm <- model.matrix(~0 + group)
png(file="voom.png")
y <- voom(d, mm, plot = T)
dev.off()

# Actual LIMMA now with empirical Bayes
fit <- lmFit(y, mm)
print(head(coef(fit)))
contr <- makeContrasts(groupMock.72.hpi - groupSARS.CoV.2.72.hpi, levels = colnames(coef(fit)))
tmp <- contrasts.fit(fit, contr)
tmp <- eBayes(tmp)

# Get most differentially expressed genes
top.table <- topTable(tmp, sort.by = "P", n = Inf)
write.csv(top.table, "total_diff.csv")

# Get number of DE genes and output them
#numDE <- length(which(top.table$adj.P.Val < 0.01))
de <- top.table[which(top.table$adj.P.Val < 0.001), ]
de <- top.table[which(abs(top.table$logFC) >= 2), ]
write.csv(de, "filtered_diff.csv")

# Run MDS using the DE genes
png(file="mds_diff.png")
plotMDS(d0[rownames(de),], col = as.numeric(group))
dev.off()