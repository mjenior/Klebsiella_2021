
metadata <- read.delim('~/Desktop/repos/Klebsiella_2021/data/metadata.tsv', sep='\t', header=TRUE)
metadata$study <- NULL
srr_files <- c('SRR1720484','SRR1720485','SRR1720486','SRR1720487','SRR1720488','SRR1720489',
               'SRR1720490','SRR1720491','SRR1720492','SRR1720493','SRR1720494','SRR1720495','SRR1720496',
               'SRR1720497','SRR1720498','SRR1720499','SRR1720500','SRR1720501','SRR1720502','SRR1720503',
               'SRR1720504','SRR1720505','SRR1720506','SRR1720507','SRR1720508','SRR1720509','SRR1720510',
               'SRR1720511','SRR1720512','SRR1720513','SRR1720514','SRR1720515','SRR1720516','SRR1720517',
               'SRR1720518','SRR1720519','SRR2147282','SRR2147283','SRR2147284','SRR2147285','SRR2147286',
               'SRR2147287','SRR2147288','SRR2147289','SRR2147290','SRR2147291','SRR2147292','SRR2147293',
               'SRR2147294','SRR2147295','SRR2147296','SRR8260120','SRR8260121','SRR8603247','SRR8603248')
size <- 20
flux_samples <- read.delim('~/Desktop/active_projects/klebsiella/data/SRR1720483/flux_samples.tsv', sep='\t', header=TRUE)
sub_sample <- sample(1:nrow(flux_samples), size, replace=FALSE)
flux_samples <- flux_samples[sub_sample,]
flux_samples$sample <- rep('SRR1720483', nrow(flux_samples))
rownames(flux_samples) <- paste0('SRR1720483_', 1:nrow(flux_samples))
flux_samples <- as.data.frame(t(flux_samples))
for (x in srr_files) {
  fluxes <- paste0('~/Desktop/active_projects/klebsiella/data/', x, '/flux_samples.tsv')
  fluxes <- read.delim(fluxes, sep='\t', header=TRUE)
  sub_sample <- sample(1:nrow(fluxes), size, replace=FALSE)
  fluxes <- fluxes[sub_sample,]
  fluxes$sample <- rep(x, nrow(fluxes))
  rownames(fluxes) <- paste0(x, '_', 1:nrow(fluxes))
  fluxes <- as.data.frame(t(fluxes))
  flux_samples <- merge(flux_samples, fluxes, by='row.names')
  rownames(flux_samples) <- flux_samples$Row.names
  flux_samples$Row.names <- NULL
}
flux_samples <- as.data.frame(t(flux_samples))
rm(fluxes)

# Prep for machine learning
flux_samples$BIOMASS_ <- NULL
flux_samples$type <- metadata$type[match(flux_samples$sample, metadata$id)]
metadata <- flux_samples[, c('sample', 'type')]
rownames(metadata) <- rownames(flux_samples)
flux_samples$sample <- NULL
lab_samples <- subset(flux_samples, type=='laboratory')
lab_samples$type <- NULL
lab_samples[] <- lapply(lab_samples, function(x) {
  if(is.factor(x)) as.numeric(as.character(x)) else x})
clinic_samples <- subset(flux_samples, type=='clinical')
clinic_samples$type <- NULL
clinic_samples[] <- lapply(clinic_samples, function(x) {
  if(is.factor(x)) as.numeric(as.character(x)) else x})

# Downsample
sub_sample <- sample(1:nrow(lab_samples), 200, replace=FALSE)
lab_samples <- lab_samples[sub_sample,]
sub_sample <- sample(1:nrow(clinic_samples), 200, replace=FALSE)
clinic_samples <- clinic_samples[sub_sample,]

# Merge data for supervised learning
lab_samples$condition <- 1
clinic_samples$condition <- 0
flux_samples <- as.data.frame(rbind(lab_samples, clinic_samples))
flux_samples$condition <- as.factor(flux_samples$condition)

# Run AUCRF and obtain feature lists
library(AUCRF)
set.seed(906801)
all_aucrf <- AUCRF(condition ~ ., data=flux_samples, pdel=0, k0=10)
#print(all_aucrf)

# Assemble feature table
top_rxns_importance <- all_aucrf$ranking[1:all_aucrf$Kopt]
rf_rxns <- as.data.frame(cbind(labels(top_rxns_importance), as.vector(top_rxns_importance)))
colnames(rf_rxns) <- c('id','mda')
rf_rxns$mda <- as.numeric(as.character(rf_rxns$mda))
rf_rxns <- rf_rxns[order(rf_rxns$mda),]
rf_rxns <- subset(rf_rxns, mda >= 1)
#write.table(rf_rxns, file='~/Desktop/repos/Klebsiella_2021/data/lab_clinic_mda.tsv', 
#            quote=FALSE, sep='\t', row.names=FALSE, col.names=TRUE)
#rf_rxns <- read.delim('~/Desktop/repos/Klebsiella_2021/data/lab_clinic_mda.tsv', sep='\t', header=TRUE)

big_wilcox <- function(group1, group2) {
  pvals <- c()
  reps <- round(length(group1) * 10)
  subsamp <- round(length(group1) / 10)
  for (x in c(1:reps)) {
    test_1 <- sample(group1, size=subsamp)
    test_2 <- sample(group2, size=subsamp)
    pvals[x] <- wilcox.test(test_1, test_2, exact=FALSE)$p.value}
  pval <- median(pvals)
  return(pval)
}

lab_valta <- abs(lab_samples[,'VALTA'])
clinic_valta <- abs(clinic_samples[,'VALTA'])
valta_pval <- big_wilcox(lab_valta, clinic_valta)

lab_argtex <- abs(lab_samples[,'ARGtex'])
clinic_argtex <- abs(clinic_samples[,'ARGtex'])
argtex_pval <- big_wilcox(lab_argtex, clinic_argtex)

# Generate figure
pdf(file='~/Desktop/repos/Klebsiella_2021/results/Figure_S1.pdf', width=8, height=4)
layout(matrix(c(1,2), nrow=1, ncol=2))

par(mar=c(3.5,3.5,1,1), mgp=c(2.2, 0.6, 0), new=FALSE, xpd=FALSE, lwd=2, las=1, lwd.tick=2)
plot(all_aucrf, xlim=c(0,250), ylim=c(0.99,1.0), cex.axis=0.7, pch=1)
box()
par(xpd=TRUE)
text(x=-60, y=1.0005, 'A', cex=1.2, font=2)
par(xpd=FALSE)

xmax <- max(rf_rxns$mda) * 1.2
par(mar=c(3,5,1,1), mgp=c(1.8, 0.6, 0), new=FALSE, xpd=FALSE, lwd=2, las=1, cex.axis=0.8)
dotchart(rf_rxns$mda,  labels=rf_rxns$id, xlab='Mean Decrease Accuracy (%)', xlim=c(0,xmax),  
         pch=16, lwd=1.7, xaxs='i', pt.cex=0.1, cex=0.8)
points(x=rf_rxns$mda, y=c(1:nrow(rf_rxns)), pch=21, cex=1.4, bg='red')
legend('bottomright', legend='MDA >= 1.0', pt.cex=0, bg='white')
par(xpd=TRUE)
lab_pos <- -(xmax*0.45)
text(x=lab_pos, y=nrow(rf_rxns)+1, 'B', cex=1.2, font=2)
par(xpd=FALSE)

dev.off()


# Subset to informative features
flux_samples <- flux_samples[, rf_rxns$id]

# Ordination analysis
library(vegan)
flux_dist <- vegdist(flux_samples, method='bray')
test <- merge(x=metadata, y=flux_samples, by='row.names')
rownames(test) <- test$Row.names
test$Row.names <- NULL
test$sample <- NULL
nmds_pval <- adonis(flux_dist ~ type, data=test, perm=999, method='bray')
nmds_pval <- nmds_pval$aov.tab[[6]][1]
nmds_pval <- as.character(round(nmds_pval, 3))

flux_nmds <- as.data.frame(metaMDS(flux_dist, k=2, trymax=25)$points)
flux_x <- (abs(max(flux_nmds$MDS1)) - abs(min(flux_nmds$MDS1))) / 2
flux_y <- (abs(max(flux_nmds$MDS2)) - abs(min(flux_nmds$MDS2))) / 2
flux_nmds$MDS1 <- flux_nmds$MDS1 - flux_x
flux_nmds$MDS2 <- flux_nmds$MDS2 - flux_y
flux_x <- max(abs(max(flux_nmds$MDS1)), abs(min(flux_nmds$MDS1))) + 0.01
flux_y <- max(abs(max(flux_nmds$MDS2)), abs(min(flux_nmds$MDS2))) + 0.01

flux_nmds <- merge(metadata, flux_nmds, by='row.names')
rownames(flux_nmds) <- flux_nmds$Row.names
flux_nmds$Row.names <- NULL
flux_nmds$sample <- NULL
clinical_nmds_points <- subset(flux_nmds, type == 'clinical')
laboratory_nmds_points <- subset(flux_nmds, type == 'laboratory')

# Generate figure
library(scales)
library(vioplot)
pdf(file='~/Desktop/repos/Klebsiella_2021/results/Figure_2.pdf', width=6, height=3)
layout(matrix(c(1,1,2,3), nrow=1, ncol=4))

par(mar=c(3.5,3.5,0.5,0.5), las=1, mgp=c(2.2,0.7,0), lwd=2, xaxt='l')
plot(x=flux_nmds$MDS1, y=flux_nmds$MDS2, xlim=c(-flux_x,flux_x), ylim=c(-flux_y, flux_y),
     xlab='NMDS Axis 1', ylab='NMDS Axis 2', pch=19, cex.lab=1.1, cex=0, cex.axis=0.9)
points(x=laboratory_nmds_points$MDS1, y=laboratory_nmds_points$MDS2, bg=alpha('darkcyan',0.8), pch=21, cex=1.7)
points(x=clinical_nmds_points$MDS1, y=clinical_nmds_points$MDS2, bg=alpha('white',0.8), pch=21, cex=1.7)
legend('topleft', legend=c('Clinical isolate','Laboratory strain'), bg='white',
       pt.bg=c('white', 'darkcyan'), pch=21, pt.cex=1.6, cex=1, box.lwd=2)
legend('bottomleft', legend=as.expression(bquote(paste(italic('p'),'-value = 0.001 ***'))), cex=0.9, bty='n')
par(xpd=TRUE)
text(x=-0.5, y=0.25, 'A', cex=1.2, font=2)
par(xpd=FALSE)

par(mar=c(2.8,3,1.7,0.5), xpd=FALSE, las=1, mgp=c(1.8,0.7,0), lwd=2, xaxt='n')
vioplot(clinic_valta, lab_valta, col=c('white', 'darkcyan'), main='Valine transaminase (VALTA)', cex.main=0.8,
        ylim=c(0, 20), ylab='Context-specific Reaction Flux', lwd=1.7, drawRect=FALSE, yaxs='i', cex.axis=0.9)
segments(x0=1, y0=17, x1=2)
text(x=1.5, y=17.7, '***', font=2, cex=2)
par(xpd=TRUE)
text(x=c(1,2), y=-1.4, labels=c('Clinical\nisolate','Laboratory\nstrain'), cex=0.9)
text(x=-0.2, y=20.8, 'B', cex=1.2, font=2)
par(xpd=FALSE)

vioplot(clinic_argtex, lab_argtex, col=c('white', 'darkcyan'), main='Argenine transport (ARGtex)', cex.main=0.8,
        ylim=c(0, 5), ylab='Context-specific Reaction Flux', lwd=1.7, drawRect=FALSE, yaxs='i', cex.axis=0.9)
segments(x0=1, y0=4, x1=2)
text(x=1.5, y=4.2, '*', font=2, cex=2)
par(xpd=TRUE)
text(x=c(1,2), y=-0.35, labels=c('Clinical\nisolate','Laboratory\nstrain'), cex=0.9)
text(x=-0.1, y=5.2, 'C', cex=1.2, font=2)
par(xpd=FALSE)

dev.off()

