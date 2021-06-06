
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

# Remove biomass objective
flux_samples$BIOMASS_ <- NULL




samples <- flux_samples$sample
flux_samples$sample <- NULL
flux_samples[] <- lapply(flux_samples, function(x) {
  if(is.factor(x)) as.numeric(as.character(x)) else x})
flux_samples <- flux_samples + abs(min(flux_samples))


library(ape)
mod <- rda(flux_samples, scale=TRUE)
biplot(mod, scaling=3, type=c('text', 'points'))


library(vegan)
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

library(scales)
pdf(file='~/Desktop/repos/Klebsiella_2021/results/flux_samples_nmds.pdf', width=4.5, height=4)
par(mar=c(3.5,3.5,0.5,0.5), las=1, mgp=c(2.2,0.7,0), lwd=2)
plot(x=flux_nmds$MDS1, y=flux_nmds$MDS2, xlim=c(-0.015,0.015), ylim=c(-0.008, 0.008),
     xlab='NMDS Axis 1', ylab='NMDS Axis 2', pch=19, cex.lab=1.1, cex=0, cex.axis=0.9)
points(x=invitro_nmds_points$MDS1, y=invitro_nmds_points$MDS2, bg=alpha('darkcyan',0.8), pch=21, cex=1.7)
points(x=invivo_nmds_points$MDS1, y=invivo_nmds_points$MDS2, bg=alpha('white',0.8), pch=21, cex=1.7)
legend('topright', legend=c('in vivo','in vitro'), text.font=3, bg='white',
       pt.bg=c('white', 'darkcyan'), pch=21, pt.cex=1.6, cex=1.1, box.lwd=2)
text(x=0.006, y=-0.0075, as.expression(bquote(paste(italic('p'),'-value = 0.01 **'))), cex=0.9, pos=4)
box()
dev.off()
