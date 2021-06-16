
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
rm(fluxes, size, sub_sample)

# Format and remove biomass objective
flux_samples <- as.data.frame(t(flux_samples))
flux_samples$BIOMASS_ <- NULL
flux_samples <- merge(flux_samples, metadata, by.x='sample', by.y='id')
samples <- flux_samples$sample
flux_samples$sample <- NULL
flux_groups <- as.factor(flux_samples$type)
flux_samples$type <- NULL
flux_samples[] <- lapply(flux_samples, function(x) {
  if(is.factor(x)) as.numeric(as.character(x)) else x})

# Limit features with machine learning
library(randomForest)
rf_obj <- randomForest(flux_groups ~ ., data=flux_samples, importance=TRUE, err.rate=TRUE, ntree=1500, mtry=15)
rf_obj <- importance(rf_obj, type=1, scale=TRUE)
rf_obj <- as.data.frame(subset(rf_obj, rf_obj > (abs(min(rf_obj)))))
flux_samples <- flux_samples[, rownames(rf_mda)]
rm(rf_obj, flux_groups)




# Calculate dissimilarities and PCoA axes
library(vegan)
library(ape)
flux_samples <- flux_samples + abs(min(flux_samples))
flux_dist <- vegdist(flux_samples, method='bray')
flux_pcoa_model <- pcoa(flux_dist, correction='none', rn=NULL)




flux_pcoa_points <- as.data.frame(flux_pcoa_model$vectors[,c(1,2)])
flux_pcoa_values <- flux_pcoa_model$values
flux_pcoa_values <- flux_pcoa_values$Relative_eig[c(1,2)] * 100.0
colnames(flux_pcoa_points) <- round(flux_pcoa_values, digits=2)
rm(flux_pcoa_model, flux_pcoa_values)


# Test difference
pval <- adonis(flux_dist ~ samples, flux_samples, perm=999)$aov.tab[[6]][1]



# Center points
flux_x <- (abs(max(flux_pcoa_points[,1])) - abs(min(flux_pcoa_points[,1]))) / 2
flux_y <- (abs(max(flux_pcoa_points[,2])) - abs(min(flux_pcoa_points[,2]))) / 2
flux_pcoa_points[,1] <- flux_pcoa_points[,1] - flux_x
flux_pcoa_points[,2] <- flux_pcoa_points[,2] - flux_y
flux_xlim <- max(abs(max(flux_pcoa_points[,1])), abs(min(flux_pcoa_points[,1]))) + 0.01
flux_ylim <- max(abs(max(flux_pcoa_points[,2])), abs(min(flux_pcoa_points[,2]))) + 0.01

# Combine with metadata
flux_pcoa_points$samples <- samples
flux_pcoa_points <- merge(flux_pcoa_points, metadata, by.x='samples', by.y='id')

# Subset points
clinical_pcoa_points <- subset(flux_pcoa_points, type == 'clinical')
laboratory_pcoa_points <- subset(flux_pcoa_points, type == 'laboratory')

# Prep for plotting
flux_pcoa_points[,2] <- flux_pcoa_points[,2] * -1
x_axis_lab <- gsub('X', '', as.character(colnames(flux_pcoa_points)[2]))
x_axis_lab <- paste ('PC1 (', x_axis_lab, '%)', sep='')
y_axis_lab <- gsub('X', '', as.character(colnames(flux_pcoa_points)[3]))
y_axis_lab <- paste ('PC2 (', y_axis_lab, '%)', sep='')

# Generate figure
clinical_col <- 'red2'
laboratory_col <- 'blue3'
library(scales)

pdf(file='~/Desktop/repos/Klebsiella_2021/results/flux_samples_pcoa.pdf', width=4.5, height=4)
par(mar=c(4.1,4.1,1,1), las=1, mgp=c(2.9,0.7,0), lwd=2, lick.lwd=2)
plot(x=flux_pcoa_points[,2], y=flux_pcoa_points[,3], xlim=c(-0.01,0.01), ylim=c(-0.01,0.01),
     xlab=x_axis_lab, ylab=y_axis_lab, pch=19, cex.lab=1.4, cex=0)
points(x=laboratory_pcoa_points[,2], y=laboratory_pcoa_points[,3], bg=alpha(laboratory_col,0.75), pch=21, cex=2.4)
points(x=clinical_pcoa_points[,2], y=clinical_pcoa_points[,3], bg=alpha(clinical_col,0.75), pch=21, cex=2.4)
legend('topleft', legend=c('Clinical isolates','Laboratory isolates'), 
       pt.bg=c(clinical_col,laboratory_col), pch=21, pt.cex=2, pt.lwd=1.5, cex=1.1, bty='n')
legend('bottomright', legend=as.expression(bquote(paste(italic('p'),'-value = 0.001***'))), 
       pch=1, cex=1.2, pt.cex=0, bty='n')
box()

# Biplot vectors




dev.off()








flux_samples$sample <- NULL
pca_mod <- rda(flux_samples, scale=TRUE)



biplot(pca_mod, scaling=3, type=c('text','points'))


pcoa_vectors <- scores(pca_mod, choices=1:2, display='sp') 

pcoa_magnitude <- as.vector(apply(pcoa_vectors, 1, function (x) sqrt((abs(x[1])^2) + (abs(x[2])^2))))
pcoa_vectors <- as.data.frame(pcoa_vectors)
pcoa_vectors$magnitude <- pcoa_magnitude
pcoa_vectors <- pcoa_vectors[order(-pcoa_vectors$magnitude),]
pcoa_vectors <- subset(pcoa_vectors, magnitude >= 1.484914)




# scores() choices= indicates which axes are to be selected, make sure to specify the scaling if its different than 2 
arrows(0, 0, pcoa_vectors[,1], pcoa_vectors[,2], length=0, lty=1, col='red')


summary(pca_mod)
screeplot(pca_mod) 



library(stats)
fit <- prcomp(flux_samples)
biplot(fit, scaling=3, type=c('text','points'))



library(FactoMineR)
res.pca <- PCA(flux_samples, graph=FALSE)










