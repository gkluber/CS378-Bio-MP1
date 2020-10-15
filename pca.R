#if (!requireNamespace('BiocManager', quietly = TRUE))
#    install.packages('BiocManager')

#install.packages("spatial")
#install.packages("stringi")

#BiocManager::install('PCAtools')

library(PCAtools)

# load in differentially expressed genes
diff = read.table("filtered_diff.csv", sep=",", header=T, row.names=1)
diff = row.names(diff)
print(diff[1:10])

# load in datasets
patients <- read.table(file = 'SraRunTable.csv', sep=',', row.names=1, header = TRUE)
genes <- read.table(file="combined_counts.tsv", sep='\t', row.names=1, header = TRUE)
genes = genes[diff,]
print(head(genes))

p <- pca(genes, metadata = patients, removeVar=.1)

print("plotting scree plot")
s <- screeplot(p, axisLabSize = 18, 
                    titleLabSize = 22,
                    title="SCREE plot", 
                    subtitle="Principal Components Composed of Differentially Expressed Genes")
png(file="scree.png")
plot(s)
dev.off()

print("plotting biplot")
png(file="biplot.png")
b <- biplot(p,
    title = "Biplot of Treated and Control",
    colby = "TREATMENT",
    colkey = c(Mock.72.hpi="red", SARS.CoV.2.72.hpi="green"))
plot(b)
dev.off()

print("plotting pairs plot")
png(file="pp.png",
    width = 800,
    height = 800)
b <- pairsplot(p,
    title = "Pairs Plot of Principal Components",
    titleLabSize = 20,
    colby = "TREATMENT",
    colkey = c(Mock.72.hpi="red", SARS.CoV.2.72.hpi="green"),
    axisLabSize=7)
plot(b)
dev.off()