
# Transpose transcriptome replicates for new riptide function
clinical <- as.data.frame(t(read.delim('/home/mjenior/Desktop/active_projects/klebsiella/data/transcript_mapping/clinical_subsampled_transript.tsv', sep='\t', header=TRUE, row.names=1)))
write.table(clinical, file='/home/mjenior/Desktop/active_projects/klebsiella/data/transcript_mapping/clinical_replicates.tsv', 
            quote=FALSE, sep='\t', row.names=TRUE, col.names=FALSE)
laboratory <- as.data.frame(t(read.delim('/home/mjenior/Desktop/active_projects/klebsiella/data/transcript_mapping/laboratory_subsampled_transript.tsv', sep='\t', header=TRUE, row.names=1)))
write.table(laboratory, file='/home/mjenior/Desktop/active_projects/klebsiella/data/transcript_mapping/laboratory_replicates.tsv', 
            quote=FALSE, sep='\t', row.names=TRUE, col.names=FALSE)
rm(clinical, laboratory)

#-----------------------------------------------------------------------------------------------------------#

# Read in data
clinical <- '~/Desktop/repos/Klebsiella_2021/data/clinical_maxfit_reps/flux_samples.tsv'
clinical <- read.delim(clinical, sep='\t', header=TRUE, row.names=1)
clinical <- as.data.frame(apply(clinical, 2, as.numeric))
rownames(clinical) <- paste0('clinical_', c(1:nrow(clinical)))
laboratory <- '~/Desktop/repos/Klebsiella_2021/data/laboratory_maxfit_reps/flux_samples.tsv'
laboratory <- read.delim(laboratory, sep='\t', header=TRUE, row.names=1)
laboratory <- as.data.frame(apply(laboratory, 2, as.numeric))
rownames(laboratory) <- paste0('laboratory_', c(1:nrow(laboratory)))


median(clinical[,'ILETA'])
median(laboratory[,'ILETA'])

median(clinical[,'ALATA_L'])
median(laboratory[,'ALATA_L'])

median(clinical[,'VPAMTr'])
median(laboratory[,'VPAMTr'])

median(clinical[,'VALTA'])
median(laboratory[,'VALTA'])

ALATA_L --> VPAMTr --> VALTA

# Create metadata table
clinical_metadata <- cbind(rownames(clinical), rep('clinical', nrow(clinical)))
colnames(clinical_metadata) <- c('sample','condition')
laboratory_metadata <- cbind(rownames(laboratory), rep('laboratory', nrow(laboratory)))
colnames(laboratory_metadata) <- c('sample','condition')
metadata <- as.data.frame(rbind(clinical_metadata, laboratory_metadata))
rm(clinical_metadata, laboratory_metadata)

# Save biomass fluxes
clinical_biomass <- clinical$BIOMASS_
laboratory_biomass <- laboratory$BIOMASS_
biomass_pval <- round(wilcox.test(clinical_biomass, laboratory_biomass, exact=FALSE)$p.value, 4)

# Separate into conserved reactions
core_rxns <- intersect(colnames(clinical), colnames(laboratory))
clinical_core <- clinical[, core_rxns]
clinical_core$BIOMASS_ <- NULL
laboratory_core <- laboratory[, core_rxns]
laboratory_core$BIOMASS_ <- NULL
rm(core_rxns)

# Assmble flux sample table
flux_samples <- as.data.frame(rbind(clinical_core, laboratory_core))
flux_samples <- flux_samples + abs(min(flux_samples))

# Calculate dissimilarity
library(vegan)
flux_dist <- vegdist(flux_samples, method='bray')

# Test difference
flux_samples$condition <- as.factor(metadata$condition)
pval <- adonis(flux_dist ~ condition, data=flux_samples, perm=999, method='bray')
pval <- pval$aov.tab[[6]][1]
pval <- as.character(round(pval, 4))

# Ordination analysis - center point
flux_nmds <- as.data.frame(metaMDS(flux_dist, k=2, trymax=50)$points)
flux_x <- (abs(max(flux_nmds$MDS1)) - abs(min(flux_nmds$MDS1))) / 2
flux_y <- (abs(max(flux_nmds$MDS2)) - abs(min(flux_nmds$MDS2))) / 2
flux_nmds$MDS1 <- flux_nmds$MDS1 - flux_x
flux_nmds$MDS2 <- flux_nmds$MDS2 - flux_y
flux_x <- max(abs(max(flux_nmds$MDS1)), abs(min(flux_nmds$MDS1))) + 0.01
flux_y <- max(abs(max(flux_nmds$MDS2)), abs(min(flux_nmds$MDS2))) + 0.01

# Subset points
flux_nmds <- merge(x=metadata, y=flux_nmds, by.x='sample', by.y='row.names')
clinical_points <- subset(flux_nmds, condition == 'clinical')
laboratory_points <- subset(flux_nmds, condition == 'laboratory')

# Calculate dentroids
clinical_centroids <- aggregate(cbind(clinical_points$MDS1, clinical_points$MDS2)~clinical_points$condition, data=clinical_points, mean)
laboratory_centroids <- aggregate(cbind(laboratory_points$MDS1, laboratory_points$MDS2)~laboratory_points$condition, data=laboratory_points, mean)

# Mean wwithin-group dissimilarity
flux_groups <- as.factor(c(rep('clinical', nrow(clinical_core)), rep('laboratory', nrow(laboratory_core))))
meandist(flux_dist, grouping=flux_groups)
#             clinical  laboratory
# clinical   0.004 0.006
# laboratory 0.006 0.001

#----------------------------------------------------------------------#

# Generate figures
pdf(file='~/Desktop/repos/Klebsiella_2021/results/Figure_3.pdf', width=5, height=3)
layout(matrix(c(1,2,2,
                1,2,2), nrow=2, ncol=3, byrow=TRUE))

library(vioplot)
#pdf(file='~/Desktop/repos/Klebsiella_2021/results/biomass.pdf', width=2.5, height=5)
par(mar=c(5,3,1,1), xpd=FALSE, las=1, mgp=c(1.8,0.7,0), lwd=1.7)
vioplot(clinical_biomass, laboratory_biomass, col=c('#B13AED','#76EEC6'), 
        ylim=c(0, 3), ylab='Sampled Biomass Flux', lwd=1.7, drawRect=FALSE, yaxs='i', yaxt='n')
axis(side=2, at=seq(0,3,0.5), cex.axis=0.8, lwd=1.7)
legend('bottomright', legend=as.expression(bquote(paste(italic('p'),'-value << 0.001'))), 
       bty='n', pt.cex=0, cex=0.8)
segments(x0=1, y0=2.7, x1=2)
text(x=1.5, y=2.8, '***', font=2, cex=1.5)
par(xpd=TRUE)
text(x=c(0.9,1.9), y=-0.45, labels=c('Clinical','Laboratory'), cex=1.2, srt=45)
text(x=-0.1, y=3, 'A', font=2, cex=1.2)
par(xpd=FALSE)
box()
#dev.off()

library(scales)
#pdf(file='~/Desktop/repos/Klebsiella_2021/results/flux_samples_nmds.pdf', width=5.5, height=5)
par(mar=c(3.5,3.5,1,1), las=1, mgp=c(2.5,0.6,0), lwd=1.7)
plot(x=flux_nmds$MDS1, y=flux_nmds$MDS2, xlim=c(-0.006,0.006), ylim=c(-0.004, 0.004),
     xlab='NMDS Axis 1', ylab='NMDS Axis 2', pch=19, cex.lab=1.1, cex=0, cex.axis=0.8)
segments(x0=clinical_points$MDS1, y0=clinical_points$MDS2, 
         x1=clinical_centroids[1,2], y1=clinical_centroids[1,3], col='ivory3')
points(x=clinical_points$MDS1, y=clinical_points$MDS2, bg=alpha('#B13AED',0.8), pch=21, cex=1.7)
segments(x0=laboratory_points$MDS1, y0=laboratory_points$MDS2, 
         x1=laboratory_centroids[1,2], y1=laboratory_centroids[1,3], col='ivory3')
points(x=laboratory_points$MDS1, y=laboratory_points$MDS2, bg=alpha('#76EEC6',0.8), pch=21, cex=1.7)
legend('topright', legend='Core Metabolism', bty='n', pt.cex=0, cex=1.3)
legend('topleft', legend=c('Clinical isolates','Laboratory strains'), 
       pt.bg=c('#B13AED','#76EEC6'), pch=21, pt.cex=1.6, cex=1.1, box.lwd=1.7)
legend('bottomright', legend=as.expression(bquote(paste(italic('p'),'-value = 0.001 ***'))), 
       bty='n', pt.cex=0, cex=0.8)
box()
par(xpd=TRUE)
text(x=-0.0082, y=0.0043, 'B', font=2, cex=1.2)
par(xpd=FALSE)
#dev.off()

dev.off()

