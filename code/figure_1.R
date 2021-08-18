
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
write.table(clinical_transcription, file='~/Desktop/active_projects/klebsiella/data/transcript_mapping/clinical_subsampled_transript.tsv', 
            quote=FALSE, sep='\t', row.names=TRUE, col.names=TRUE)
write.table(laboratory_transcription, file='~/Desktop/active_projects/klebsiella/data/transcript_mapping/laboratory_subsampled_transript.tsv', 
            quote=FALSE, sep='\t', row.names=TRUE, col.names=TRUE)

# Calculate median values and save
clinical_median <- as.vector(apply(clinical_transcription, 2, median))
clinical_median <- as.data.frame(cbind(colnames(clinical_transcription), clinical_median))
write.table(clinical_median, file='~/Desktop/active_projects/klebsiella/data/transcript_mapping/clinical_median_transript.tsv', 
            quote=FALSE, sep='\t', row.names=FALSE, col.names=FALSE)
laboratory_median <- as.vector(apply(laboratory_transcription, 2, median))
laboratory_median <- as.data.frame(cbind(colnames(laboratory_transcription), laboratory_median))
write.table(laboratory_median, file='~/Desktop/active_projects/klebsiella/data/transcript_mapping/laboratory_median_transript.tsv', 
            quote=FALSE, sep='\t', row.names=FALSE, col.names=FALSE)
rm(clinical_median, laboratory_median)

#---------------------------------------------# 

# Read in preprocessed data here to save time
clinical_transcription <- read.delim('~/Desktop/active_projects/klebsiella/data/transcript_mapping/clinical_subsampled_transript.tsv', 
                                     sep='\t', header=TRUE, row.names=1)
laboratory_transcription <- read.delim('~/Desktop/active_projects/klebsiella/data/transcript_mapping/laboratory_subsampled_transript.tsv', 
                                     sep='\t', header=TRUE, row.names=1)

# Generate summary stat table
library(matrixTests)
summ_table <- col_t_welch(laboratory_transcription, clinical_transcription,
                          alternative='two.sided', mu=0, conf.level=0.95)
summ_table$p.corr <- p.adjust(summ_table$pvalue, method='BH')
summ_table$FoldChange <- summ_table$mean.x/ summ_table$mean.y
summ_table$log2.Fold.Change <- as.numeric(log2(summ_table$FoldChange))
#write.table(summ_table, file='~/Desktop/repos/Klebsiella_2021/results/Table_S1.tsv', 
#            quote=FALSE, sep='\t', row.names=TRUE, col.names=TRUE)

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

#-----------------------------------------------------------------------------------------#

# Generate figure
library(EnhancedVolcano)
library(magrittr)

quick_barchart <- function(gene, name='', panel='') {
  clinical <- as.vector(clinical_transcription[,gene])
  laboratory <- as.vector(laboratory_transcription[,gene])
  pval <- wilcox.test(clinical, laboratory, exact=FALSE)$p.value
  
  clinical <- log2(clinical + 1)
  clinical_median <- median(clinical)
  clinical_q25 <- as.numeric(quantile(clinical, 0.25))
  clinical_q75 <- as.numeric(quantile(clinical, 0.75))
  laboratory <- log2(laboratory + 1)
  laboratory_median <- median(laboratory)
  laboratory_q25 <- as.numeric(quantile(laboratory, 0.25))
  laboratory_q75 <- as.numeric(quantile(laboratory, 0.75))
  
  ymax <- round(max(c(clinical, laboratory)) * 1.2)
  yjump <- ymax / 5
  lab_dist <- -(ymax * 0.12)
  sig_pos <- (ymax + lab_dist)
  name <- paste0(name, '\n(', gene, ')')
  panel_lab <- ymax * 1.1
  
  par(mar=c(2.5,3,2,1), xpd=FALSE, las=1, mgp=c(1.9,0.8,0), lwd=1.5, xaxt='n')
  barplot(c(clinical_median, laboratory_median), xlim=c(0,2.6), ylim=c(0,ymax), ylab='Transcript (Log2)', 
          col=c('#B13AED','#76EEC6'), cex.axis=0.7, yaxt='n', cex.lab=0.7, main=name, cex.main=0.6)
  axis(side=2, at=seq(0,ymax,yjump), cex.axis=0.6, lwd=1.5)
  box()
  if (clinical_median != 0) {
    segments(x0=0.7, y0=clinical_q25, y1=clinical_q75)
    segments(x0=0.9, y0=clinical_q25, x1=0.5)
    segments(x0=0.9, y0=clinical_q75, x1=0.5)
  }
  if (laboratory_median != 0) {
    segments(x0=1.9, y0=laboratory_q25, y1=laboratory_q75)
    segments(x0=2.1, y0=laboratory_q25, x1=1.7)
    segments(x0=2.1, y0=laboratory_q75, x1=1.7)
  }
  par(xpd=TRUE)
  text(x=c(0.7,1.9), y=lab_dist, labels=c('Clinical\nisolates','Laboratory\nstrains'), cex=0.8)
  text(x=-0.85, y=panel_lab, font=2, labels=panel, cex=1.1)
  par(xpd=FALSE)
  
  
  if (pval < 0.001) {
    segments(x0=0.7, x1=1.9, y0=sig_pos)
    text(x=1.3, y=sig_pos*1.05, '***', cex=1.3, font=2)
    pval <- 1
  }
  if (pval <= 0.01) {
    segments(x0=0.7, x1=1.9, y0=sig_pos)
    text(x=1.3, y=sig_pos*1.05, '**', cex=1.3, font=2)
    pval <- 1
  }
  if (pval <= 0.05) {
    segments(x0=0.7, x1=1.9, y0=sig_pos)
    text(x=1.3, y=sig_pos*1.05, '*', cex=1.3, font=2)
    pval <- 1
  }
  }

pdf(file='~/Desktop/repos/Klebsiella_2021/results/Figure_1A.pdf', width=8, height=5)
EnhancedVolcano(summ_table,
                title = '',
                subtitle = '',
                lab = rownames(summ_table),
                axisLabSize = 10,
                x = 'log2.Fold.Change',
                y = 'pvalue',
                xlim = c(-5,5),
                ylim = c(0,10),
                xlab = bquote(~Log[2]~'fold change'),
                pCutoff = 0.05,
                FCcutoff = 2,
                pointSize = 3.0,
                selectLab = c('KPN_03340','KPN_01151','KPN_00435','KPN_04433','KPN_00795','KPN_pKPN4p07058','KPN_03976','KPN_03279','KPN_00656','KPN_03160'),
                labSize = 2.5,
                labCol = 'black',
                labFace = 'bold',
                col = c('black', 'chartreuse4', 'blue4', 'firebrick'),
                legendLabSize = 10,
                legendIconSize = 4.0,
                colAlpha = 4/5,
                drawConnectors = FALSE,
                colConnectors = 'black')
dev.off()


# Barcharts for genes of interest
pdf(file='~/Desktop/repos/Klebsiella_2021/results/Figure_1B-E.pdf', width=6, height=2)
layout(matrix(c(1,2,3,4), nrow=1, ncol=4, byrow=TRUE))
quick_barchart('KPN_00795', 'hutC: Transcription Factor', 'B')
quick_barchart('KPN_03340', 'gcvH: Glycine cleavage system', 'C')
text(x=2, y=12, font=2, labels='A', cex=1.1)
quick_barchart('KPN_pKPN4p07058', 'Aminoglycoside N-acetyltransferase', 'D')
quick_barchart('KPN_04433', 'Putative stress-response protein', 'E')
dev.off()

pdf(file='~/Desktop/repos/Klebsiella_2021/results/Figure_S1.pdf', width=5, height=1.5)
layout(matrix(c(1,2,3), nrow=1, ncol=3, byrow=TRUE))
quick_barchart('KPN_03279', 'mrkA: Type-3 fimbriae', 'A')
quick_barchart('KPN_00435', 'rpmE2: 50S ribosomal protein L31', 'B')
quick_barchart('KPN_03976', 'rpmG: 50S ribosomal protein L33', 'C')
dev.off()


# zur regulon
# KPN_00435	rpmE2: 50S ribosomal protein L31
# KPN_03976	rpmG: 50S ribosomal protein L33
