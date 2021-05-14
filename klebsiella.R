
metadata <- read.delim('~/Desktop/active_projects/klebsiella/metadata.tsv', sep='\t', header=TRUE)
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
flux_samples$type <- metadata$type[match(flux_samples$sample, metadata$id)]
metadata <- flux_samples[, c('sample', 'type')]
rownames(metadata) <- rownames(flux_samples)
flux_samples$sample <- NULL
condition <- flux_samples$type
invitro_samples <- subset(flux_samples, type=='in_vitro')
invitro_samples$type <- NULL
invitro_samples <- as.data.frame(apply(invitro_samples, 2, as.numeric))
invivo_samples <- subset(flux_samples, type=='in_vivo')
invivo_samples$type <- NULL
invivo_samples <- as.data.frame(apply(invivo_samples, 2, as.numeric))
flux_samples$type <- NULL
flux_samples <- as.data.frame(apply(flux_samples, 2, as.numeric))

# Merge data for supervised learning
invitro_samples$condition <- 1
invivo_samples$condition <- 0
all_samples <- rbind(invitro_samples, invivo_samples)
all_samples$condition <- as.factor(all_samples$condition)

# Run AUCRF and obtain feature lists
library(AUCRF)
set.seed(906801)
all_aucrf <- AUCRF(condition ~ ., data=all_samples, pdel=0, k0=10)
print(all_aucrf)
rm(all_samples)

# Assemble feature table
top_rxns_importance <- all_aucrf$ranking[1:all_aucrf$Kopt]
rf_rxns <- as.data.frame(cbind(labels(top_rxns_importance), as.vector(top_rxns_importance)))
colnames(rf_rxns) <- c('id','mda')
rm(all_aucrf, top_rxns_importance)









library(vioplot)
plot_flux <- function(reaction, rxn_name='test', ymax=20) {
  invitro_flux <- invitro_samples[,reaction]
  invivo_flux <- invivo_samples[,reaction]
  if (rxn_name == 'test') {rxn_name = reaction}
  
  pvals <- c()
  for (x in c(1:2500)) {
    test_1 <- sample(invitro_flux, size=25)
    test_2 <- sample(invivo_flux, size=25)
    pvals[x] <- round(wilcox.test(test_1, test_2, exact=FALSE)$p.value, 3)}
  pval <- median(pvals)
  print(pval)
  
  par(mar=c(2,3,1.5,0.5), xpd=FALSE, las=1, mgp=c(1.8,0.7,0), lwd=2)
  vioplot(invitro_flux, invivo_flux, col=c('blue3', 'darkorange2'), main=rxn_name, cex.main=1,
          ylim=c(-ymax, ymax), ylab='Predicted Flux', lwd=1.7, drawRect=FALSE, yaxt='n', yaxs='i')
  axis(side=2, at=seq(-ymax, ymax, ymax/4), cex.axis=0.7, lwd=2)
  abline(h=0, lwd=2, lty=2, col='gray25')
  vioplot(invitro_flux, invivo_flux, col=c('blue3', 'darkorange2'), main=rxn_name, cex.main=1, add=TRUE,
          ylim=c(-ymax, ymax), ylab='Predicted Flux', lwd=1.7, drawRect=FALSE, yaxt='n', yaxs='i')
  box(lwd=2)
}


png(filename='~/Desktop/active_projects/klebsiella/results/aucrf_mda.png', units='in', width=3.5, height=5, res=300)
rf_results <- read.delim('~/Desktop/active_projects/klebsiella/data/rf.mda.invitro_invivo.tsv', sep='\t', header=TRUE)
rf_results <- rf_results[order(rf_results$mda),] 
rf_results$name <- gsub('_', ' ', rf_results$name)
par(mar=c(3, 0.5, 0.5, 0.5), mgp=c(1.4, 0.5, 0), xpd=FALSE, lwd=1.7)
dotchart(rf_results$mda, bg='chartreuse3', xlim=c(0,150),  
         pch=21, lwd=1.7, pt.cex=1.5, cex=0.8)
text(x=-1, y=seq(1.4,11.4,1), labels=rf_results$name, cex=0.8, pos=4)
mtext('Mean Decrease Accuracy', side=1, padj=2.5, cex=0.9)
dev.off()




png(filename='~/Desktop/active_projects/klebsiella/results/VALTA.png', 
    units='in', width=3, height=4, res=300)
plot_flux('VALTA', rxn_name='Valine transaminase')
segments(x0=1, y0=10, x1=2)
text(x=1.5, y=12, '***', font=2, cex=1.3) # 0
par(xpd=TRUE)
text(x=c(1,2), y=-22.5, labels=c('K. pneumoniae\nin vitro','K. pneumoniae\nin vivo'), cex=0.8, font=3)
par(xpd=FALSE)
dev.off()

png(filename='~/Desktop/active_projects/klebsiella/results/SUCDi.png', 
    units='in', width=3, height=4, res=300)
plot_flux('SUCDi', rxn_name='Succinate dehydrogenase')
segments(x0=1, y0=15, x1=2)
text(x=1.5, y=17, '***', font=2, cex=1.3) # 0
par(xpd=TRUE)
text(x=c(1,2), y=-22.5, labels=c('K. pneumoniae\nin vitro','K. pneumoniae\nin vivo'), cex=0.8, font=3)
par(xpd=FALSE)
dev.off()

png(filename='~/Desktop/active_projects/klebsiella/results/FUM.png', 
    units='in', width=3, height=4, res=300)
plot_flux('FUM', rxn_name='Fumarase')
segments(x0=1, y0=15, x1=2)
text(x=1.5, y=17, '***', font=2, cex=1.3) # 0
par(xpd=TRUE)
text(x=c(1,2), y=-22.5, labels=c('K. pneumoniae\nin vitro','K. pneumoniae\nin vivo'), cex=0.8, font=3)
par(xpd=FALSE)
dev.off()

png(filename='~/Desktop/active_projects/klebsiella/results/PDH.png', 
    units='in', width=3, height=4, res=300)
plot_flux('PDH', rxn_name='Pyruvate dehydrogenase')
segments(x0=1, y0=15, x1=2)
text(x=1.5, y=17, '***', font=2, cex=1.3) # 0
par(xpd=TRUE)
text(x=c(1,2), y=-22.5, labels=c('K. pneumoniae\nin vitro','K. pneumoniae\nin vivo'), cex=0.8, font=3)
par(xpd=FALSE)
dev.off()

png(filename='~/Desktop/active_projects/klebsiella/results/ADK1.png', 
    units='in', width=3, height=4, res=300)
plot_flux('ADK1', rxn_name='Adenylate kinase', ymax=5)
segments(x0=1, y0=2.5, x1=2)
text(x=1.5, y=3, '**', font=2, cex=1.3) # 0.002
par(xpd=TRUE)
text(x=c(1,2), y=-5.6, labels=c('K. pneumoniae\nin vitro','K. pneumoniae\nin vivo'), cex=0.8, font=3)
par(xpd=FALSE)
dev.off()

png(filename='~/Desktop/active_projects/klebsiella/results/PPK.png', 
    units='in', width=3, height=4, res=300)
plot_flux('PPK', rxn_name='Polyphosphate kinase', ymax=5)
segments(x0=1, y0=2.5, x1=2)
text(x=1.5, y=3, '*', font=2, cex=1.3) # 0.012
par(xpd=TRUE)
text(x=c(1,2), y=-5.6, labels=c('K. pneumoniae\nin vitro','K. pneumoniae\nin vivo'), cex=0.8, font=3)
par(xpd=FALSE)
dev.off()

png(filename='~/Desktop/active_projects/klebsiella/results/EX_arg__L_e.png', 
    units='in', width=3, height=4, res=300)
plot_flux('EX_arg__L_e', rxn_name='L-Arginine exchange', ymax=5)
segments(x0=1, y0=2.5, x1=2) 
text(x=1.5, y=3, '*', font=2, cex=1.3) # 0.011
par(xpd=TRUE)
text(x=c(1,2), y=-5.6, labels=c('K. pneumoniae\nin vitro','K. pneumoniae\nin vivo'), cex=0.8, font=3)
par(xpd=FALSE)
dev.off()

png(filename='~/Desktop/active_projects/klebsiella/results/MDH.png', 
    units='in', width=3, height=4, res=300)
plot_flux('MDH', rxn_name='Malate dehydrogenase', ymax=10)
segments(x0=1, y0=8, x1=2)
text(x=1.5, y=9, '*', font=2, cex=1.3) # 0.015
par(xpd=TRUE)
text(x=c(1,2), y=-11.25, labels=c('K. pneumoniae\nin vitro','K. pneumoniae\nin vivo'), cex=0.8, font=3)
par(xpd=FALSE)
dev.off()













metadata <- read.delim('~/Desktop/active_projects/klebsiella/metadata.tsv', sep='\t', header=TRUE)
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
flux_samples <- read.delim('~/Desktop/active_projects/klebsiella/data/SRR1720483/flux_samples.tsv', sep='\t', header=TRUE, row.names=1)
sub_sample <- sample(1:nrow(flux_samples), size, replace=FALSE)
flux_samples <- flux_samples[sub_sample,]
flux_samples$sample <- rep('SRR1720483', nrow(flux_samples))
rownames(flux_samples) <- paste0('SRR1720483_', 1:nrow(flux_samples))
flux_samples <- as.data.frame(t(flux_samples))
for (x in srr_files) {
  fluxes <- paste0('~/Desktop/active_projects/klebsiella/data/', x, '/flux_samples.tsv')
  fluxes <- read.delim(fluxes, sep='\t', header=TRUE, row.names=1)
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
rm(fluxes, size, sub_sample)

samples <- flux_samples$sample
flux_samples$sample <- NULL
flux_samples[] <- lapply(flux_samples, function(x) {
  if(is.factor(x)) as.numeric(as.character(x)) else x})
flux_samples <- flux_samples + abs(min(flux_samples))
flux_dist <- vegdist(flux_samples, method='bray')
flux_samples$sample <- samples
rm(samples)

test <- merge(x=metadata, y=flux_samples, by.x='id', by.y='sample')
test$id <- NULL
test$sample <- NULL
rownames(test) <- rownames(flux_samples)
nmds_pval <- adonis(flux_dist ~ type, data=test, perm=99, method='bray')
nmds_pval <- nmds_pval$aov.tab[[6]][1]
nmds_pval <- as.character(round(nmds_pval, 3))
test <- test[,c('type','ADK1')]

flux_nmds <- as.data.frame(metaMDS(flux_dist, k=2, trymax=25)$points)
flux_x <- (abs(max(flux_nmds$MDS1)) - abs(min(flux_nmds$MDS1))) / 2
flux_y <- (abs(max(flux_nmds$MDS2)) - abs(min(flux_nmds$MDS2))) / 2
flux_nmds$MDS1 <- flux_nmds$MDS1 - flux_x
flux_nmds$MDS2 <- flux_nmds$MDS2 - flux_y
flux_x <- max(abs(max(flux_nmds$MDS1)), abs(min(flux_nmds$MDS1))) + 0.01
flux_y <- max(abs(max(flux_nmds$MDS2)), abs(min(flux_nmds$MDS2))) + 0.01



flux_nmds <- merge(test, flux_nmds, by='row.names')
rownames(flux_nmds) <- flux_nmds$Row.names
flux_nmds$Row.names <- NULL
flux_nmds$ADK1 <- NULL


invivo_nmds_points <- subset(flux_nmds, type == 'in_vivo')
invitro_nmds_points <- subset(flux_nmds, type == 'in_vitro')


pdf(file='~/Desktop/active_projects/klebsiella/flux_samples_nmds.pdf', width=4.5, height=4)
par(mar=c(3.5,3.5,0.5,0.5), las=1, mgp=c(2.2,0.7,0), lwd=2)
plot(x=flux_nmds$MDS1, y=flux_nmds$MDS2, xlim=c(-0.015,0.015), ylim=c(-0.008, 0.008),
     xlab='NMDS Axis 1', ylab='NMDS Axis 2', pch=19, cex.lab=1.1, cex=0, cex.axis=0.9)
points(x=invitro_nmds_points$MDS1, y=invitro_nmds_points$MDS2, bg=alpha('darkcyan',0.8), pch=21, cex=1.7)
points(x=invivo_nmds_points$MDS1, y=invivo_nmds_points$MDS2, bg=alpha('white',0.8), pch=21, cex=1.7)
legend('topright', legend=c('in vivo','in vitro'), text.font=3,
       pt.bg=c('white', 'darkcyan'), pch=21, pt.cex=1.6, cex=1.1, box.lwd=2)
text(x=0.006, y=-0.0075, as.expression(bquote(paste(italic('p'),'-value = 0.01 **'))), cex=0.9, pos=4)
box()
dev.off()











