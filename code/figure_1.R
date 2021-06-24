
# Transcriptomic analysis
metadata <- read.delim('~/Desktop/repos/Klebsiella_2021/data/metadata.tsv', sep='\t', header=TRUE)
clinical_srr <- as.vector(subset(metadata, type == 'clinical')[,1])
laboratory_srr <- as.vector(subset(metadata, type == 'laboratory')[,1])
rm(metadata)

current <- paste0('~/Desktop/active_projects/klebsiella/data/transcript_mapping/', clinical_srr[1], '.mapped.tsv')
clinical_transcription <- read.delim(current, sep='\t', header=FALSE, row.names=1)
colnames(clinical_transcription) <- paste0(clinical_srr[1],'_', 1:ncol(clinical_transcription))
for (x in 2:length(clinical_srr)) {
  current <- paste0('~/Desktop/active_projects/klebsiella/data/transcript_mapping/', clinical_srr[x], '.mapped.tsv')
  current <- read.delim(current, sep='\t', header=FALSE, row.names=1)
  colnames(current) <- paste0(clinical_srr[x],'_', 1:ncol(current))
  clinical_transcription <- merge(clinical_transcription, current, by='row.names')
  rownames(clinical_transcription) <- clinical_transcription$Row.names
  clinical_transcription$Row.names <- NULL}
rm(clinical_srr, x, current)
clinical_transcription <- as.data.frame(t(clinical_transcription))
clinical_transcription$`*` <- NULL

current <- paste0('~/Desktop/active_projects/klebsiella/data/transcript_mapping/', laboratory_srr[1], '.mapped.tsv')
laboratory_transcription <- read.delim(current, sep='\t', header=FALSE, row.names=1)
colnames(laboratory_transcription) <- paste0(laboratory_srr[1],'_', 1:ncol(laboratory_transcription))
for (x in 2:length(laboratory_srr)) {
  current <- paste0('~/Desktop/active_projects/klebsiella/data/transcript_mapping/', laboratory_srr[x], '.mapped.tsv')
  current <- read.delim(current, sep='\t', header=FALSE, row.names=1)
  colnames(current) <- paste0(laboratory_srr[x],'_', 1:ncol(current))
  laboratory_transcription <- merge(laboratory_transcription, current, by='row.names')
  rownames(laboratory_transcription) <- laboratory_transcription$Row.names
  laboratory_transcription$Row.names <- NULL}
rm(laboratory_srr, x, current)
laboratory_transcription <- as.data.frame(t(laboratory_transcription))
laboratory_transcription$`*` <- NULL

# Remove low abundance samples
clinical_transcription <- as.data.frame(t(clinical_transcription))
laboratory_transcription <- as.data.frame(t(laboratory_transcription))
min_samp1 <- round(median(as.vector(colSums(clinical_transcription))) * 0.8)
min_samp2 <- round(median(as.vector(colSums(laboratory_transcription))) * 0.8)
min_samp <- min(c(min_samp1, min_samp2))
clinical_transcription <- clinical_transcription[,which(as.vector(colSums(clinical_transcription)) >= min_samp)]
laboratory_transcription <- laboratory_transcription[,which(as.vector(colSums(laboratory_transcription)) >= min_samp)]
rm(min_samp, min_samp1, min_samp2)

# Subsample data
library(vegan)
sub_samp1 <- round(median(as.vector(colSums(clinical_transcription))) * 0.8)
sub_samp2 <- round(median(as.vector(colSums(laboratory_transcription))) * 0.8)
sub_samp <- min(c(sub_samp1, sub_samp2))
for (x in 1:ncol(clinical_transcription)) {clinical_transcription[,x] <- as.vector(rrarefy(clinical_transcription[,x], sample=sub_samp))}
for (x in 1:ncol(laboratory_transcription)) {laboratory_transcription[,x] <- as.vector(rrarefy(laboratory_transcription[,x], sample=sub_samp))}
clinical_transcription <- as.data.frame(t(clinical_transcription))
laboratory_transcription <- as.data.frame(t(laboratory_transcription))
rm(sub_samp, sub_samp1, sub_samp2, x)

# Add gene names and aggregate
#gene_dict <- read.delim('~/Desktop/active_projects/klebsiella/genomes/Kpneumoniae_mgh78578.genes.tsv', sep='\t', header=TRUE)
#clinical_transcription <- merge(gene_dict, clinical_transcription, by.x='gene', by.y='row.names')
#clinical_transcription$gene <- NULL
#clinical_transcription <- aggregate(. ~ description, data=clinical_transcription, FUN=sum)
#rownames(clinical_transcription) <- clinical_transcription$description
#clinical_transcription$description <- NULL
#clinical_transcription <- as.data.frame(t(clinical_transcription))
#laboratory_transcription <- merge(gene_dict, laboratory_transcription, by.x='gene', by.y='row.names')
#laboratory_transcription$gene <- NULL
#laboratory_transcription <- aggregate(. ~ description, data=laboratory_transcription, FUN=sum)
#laboratory_transcription$description <- NULL
#rownames(laboratory_transcription) <- laboratory_transcription$description
#laboratory_transcription <- as.data.frame(t(laboratory_transcription))
#rm(gene_dict)

# Generate summary stat table
summ_table <- col_t_welch(laboratory_transcription, clinical_transcription,
                          alternative='two.sided', mu=0, conf.level=0.95)
summ_table$p.corr <- p.adjust(summ_table$pvalue, method='BH')
summ_table$FoldChange <- summ_table$mean.x/ summ_table$mean.y
summ_table$log2.Fold.Change <- as.numeric(log2(summ_table$FoldChange))

#-----------------------------------------------------------------------------------------#

# Generate figure
library(EnhancedVolcano)
library(magrittr)
library(matrixTests)


pdf(file='~/Desktop/repos/Klebsiella_2021/results/Figure_1.pdf', width=8, height=6)

EnhancedVolcano(summ_table,
                lab = rownames(summ_table),
                x = 'log2.Fold.Change',
                y = 'pvalue',
                xlim = c(-6,6),
                ylim = c(0,10),
                xlab = bquote(~Log[2]~'fold change'),
                pCutoff = 0.05,
                FCcutoff = 1.0,
                pointSize = 4.0,
                labSize = 6.0,
                colAlpha = 1,
                legendPosition = 'right',
                legendLabSize = 12,
                legendIconSize = 4.0,
                drawConnectors = TRUE,
                widthConnectors = 0.75)

dev.off()


