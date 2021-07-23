

clinical_fluxes <- read.delim('~/Desktop/repos/Klebsiella_2021/data/clinical_maxfit/flux_samples.tsv', 
                              sep='\t', header=TRUE, row.names=1)
sub_sample <- sample(1:nrow(clinical_fluxes), 250, replace=FALSE)
clinical_fluxes <- clinical_fluxes[sub_sample,]
rownames(clinical_fluxes) <- paste0('clinical_', 1:nrow(clinical_fluxes))

laboratory_fluxes <- read.delim('~/Desktop/repos/Klebsiella_2021/data/laboratory_maxfit/flux_samples.tsv', 
                                sep='\t', header=TRUE, row.names=1)
sub_sample <- sample(1:nrow(laboratory_fluxes), 250, replace=FALSE)
laboratory_fluxes <- laboratory_fluxes[sub_sample,]
rownames(laboratory_fluxes) <- paste0('laboratory_', 1:nrow(laboratory_fluxes))
rm(sub_sample)




big_wilcox <- function(group1, group2) {
  print(paste('Median 1:', median(group1)))
  print(paste('Median 2:', median(group2)))
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

# Test specific reactions
lab_valtex <- laboratory_fluxes[,'VALtex']
clinic_valtex <- clinical_fluxes[,'VALtex']
valtex_pval <- big_wilcox(lab_valta, clinic_valta)


lab_HIStex <- laboratory_fluxes[,'HIStex']
clinic_HIStex <- clinical_fluxes[,'HIStex']
HIStex_pval <- big_wilcox(lab_HIStex, clinic_HIStex)



# Find unique reactions
lab_only_rxns <- setdiff(colnames(laboratory_fluxes), colnames(clinical_fluxes))
clinic_only_rxns <- setdiff(colnames(clinical_fluxes), colnames(laboratory_fluxes))


# Subset to shared reactions
shared_rxns <- intersect(colnames(laboratory_fluxes), colnames(clinical_fluxes))
laboratory_fluxes <- laboratory_fluxes[,shared_rxns]
clinical_fluxes <- clinical_fluxes[,shared_rxns]
rm(shared_rxns)


# Merge data for supervised learning
laboratory_fluxes$condition <- 1
clinical_fluxes$condition <- 0
flux_samples <- as.data.frame(rbind(laboratory_fluxes, clinical_fluxes))
flux_samples$condition <- as.factor(flux_samples$condition)
flux_samples$BIOMASS_ <- NULL

# Run AUCRF and obtain feature lists
library(AUCRF)
set.seed(906801)
Kopt <- 15
all_aucrf <- AUCRF(condition ~ ., data=flux_samples, pdel=0, k0=Kopt)
plot(all_aucrf)

# Assemble feature table
top_rxns_importance <- all_aucrf$ranking[1:100]
rf_rxns <- as.data.frame(cbind(labels(top_rxns_importance), as.vector(top_rxns_importance)))
colnames(rf_rxns) <- c('id','mda')
rf_rxns$mda <- as.numeric(as.character(rf_rxns$mda))
rf_rxns <- subset(rf_rxns, mda >= 5.0)
rf_rxns <- rf_rxns[order(-rf_rxns$mda),]



