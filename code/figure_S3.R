

# Read in preprocessed data here to save time
clinical_transcription <- read.delim('~/Desktop/active_projects/klebsiella/data/transcript_mapping/clinical_subsampled_transript.tsv', 
                                     sep='\t', header=TRUE, row.names=1)
laboratory_transcription <- read.delim('~/Desktop/active_projects/klebsiella/data/transcript_mapping/laboratory_subsampled_transript.tsv', 
                                       sep='\t', header=TRUE, row.names=1)
clinical_valta <- clinical_transcription$KPN_04269
laboratory_valta <- laboratory_transcription$KPN_04269
rm(clinical_transcription, laboratory_transcription)

pval <- round(wilcox.test(clinical_valta, laboratory_valta)$p.value, 4)
clinical_valta <- log2(clinical_valta)
laboratory_valta <- log2(laboratory_valta)

pdf(file='~/Desktop/repos/Klebsiella_2021/results/figures/Figure_S3.pdf', width=4, height=5)
par(mar=c(3,4,2,1), las=1, mgp=c(2.5,0.7,0), xpd=FALSE, lwd=2)
plot(0, type='n', ylim=c(0,15), xlim=c(0,2), main='Valine Tranaminase (KPN_04269)',
     ylab='Transcript Abundance (Log2)', xlab='', xaxt='n', yaxt='n', cex.lab=1.4)
axis(2, at=seq(0,15,3), lwd=2, cex.axis=0.9)
boxplot(clinical_valta, cex=0, lwd=3, at=0.5, col='#B13AED', ylab='', staplewex=0.6, 
        boxwex=1, lty=1, medlwd=3, xaxt='n', yaxt='n', add=TRUE)
boxplot(laboratory_valta, cex=0, lwd=3, at=1.5, col='#76EEC6', ylab='', staplewex=0.6, 
        boxwex=1, lty=1, medlwd=3, xaxt='n', yaxt='n', add=TRUE)
mtext(c('Clinical', 'Laboratory'), side=1, at=c(0.5,1.5), padj=1, cex=1.2)
segments(x0=0.5, y0=14, x1=1.5)
text(x=1, y=14.6, '***', font=2, cex=1.5)
dev.off()
