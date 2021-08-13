### Reaction Essentiality
# Loading rxn essentiality files

SRR = c('SRR8603247', 'SRR8603248', 'SRR1720519', 'SRR1720518', 'SRR1720517', 'SRR1720516', 'SRR1720515', 'SRR1720514', 'SRR1720513', 'SRR1720512', 'SRR1720511', 'SRR1720510', 'SRR1720509', 'SRR1720508', 'SRR1720507', 'SRR1720506', 'SRR1720505', 'SRR1720504', 'SRR1720503', 'SRR1720502', 'SRR1720501', 'SRR1720500', 'SRR1720499', 'SRR1720498', 'SRR1720497', 'SRR1720496', 'SRR1720495', 'SRR1720494', 'SRR1720493', 'SRR1720492', 'SRR1720491', 'SRR1720490', 'SRR1720489', 'SRR1720488', 'SRR1720487', 'SRR1720486', 'SRR1720485', 'SRR1720484', 'SRR1720483', 'SRR8260120', 'SRR8260121', 'SRR2147296', 'SRR2147295', 'SRR2147294', 'SRR2147293', 'SRR2147292', 'SRR2147291', 'SRR2147290', 'SRR2147289', 'SRR2147288', 'SRR2147287', 'SRR2147286', 'SRR2147285', 'SRR2147284', 'SRR2147283', 'SRR2147282')

myList = c()

for (x in SRR) {
  x1 = paste('~/Papin Lab/Research Data/K. pneumoniae/Reaction Essentiality Output/', x, '.txt', sep="")
  nam = paste(x)
  myList = c(myList, assign(nam, as.vector(read.delim(toString(x1), sep=',', header=FALSE)$V1)))
}

# make dataframe

typeof(SRR1720483)

overall = unique(myList)

typeof(overall)

essential = matrix(0, nrow=length(overall), ncol=length(SRR))

colnames(essential) = (SRR)
rownames(essential) = (overall)

numb_SRR = list(SRR8603247, SRR8603248, SRR1720519, SRR1720518, SRR1720517, SRR1720516, SRR1720515, SRR1720514, SRR1720513, SRR1720512, SRR1720511, SRR1720510, SRR1720509, SRR1720508, SRR1720507, SRR1720506, SRR1720505, SRR1720504, SRR1720503, SRR1720502, SRR1720501, SRR1720500, SRR1720499, SRR1720498, SRR1720497, SRR1720496, SRR1720495, SRR1720494, SRR1720493, SRR1720492, SRR1720491, SRR1720490, SRR1720489, SRR1720488, SRR1720487, SRR1720486, SRR1720485, SRR1720484, SRR1720483, SRR8260120, SRR8260121, SRR2147296, SRR2147295, SRR2147294, SRR2147293, SRR2147292, SRR2147291, SRR2147290, SRR2147289, SRR2147288, SRR2147287, SRR2147286, SRR2147285, SRR2147284, SRR2147283, SRR2147282)

typeof(numb_SRR)

n = 1
match = 0 

for (t in 1:length(numb_SRR)) {
  temp_SRR = matrix(unlist((numb_SRR[t])))
  print(t)
  for (i in 1:length(overall)) {
    if (overall[i] %in% temp_SRR) {
      essential[i,n] = 1
      match = match + 1 
    }
    else {
      essential[i,n] = 0 
    }
  }
  n = n + 1 
}

write.table(essential, file = '~/Papin Lab/Research Data/K. pneumoniae/Figures/essential_rxns.txt', sep="\t")

# separating core and accessory essential genes

core_essential = matrix(, ncol=ncol(essential))
accessory_essential = matrix(, ncol=ncol(essential))
row_names_core = matrix()
row_names_accessory = matrix()

for (i in 1:nrow(essential)) {
  if (sum(essential[i ,])==ncol(essential)) {
    core_essential = rbind(core_essential, essential[i ,])
    row_names_core = rbind(row_names_core, overall[i])
  }
  else {
    accessory_essential = rbind(accessory_essential, essential[i ,])
    row_names_accessory = rbind(row_names_accessory, overall[i])
  }
}

rownames(core_essential) = row_names_core
row_names_core = row_names_core[-1 ,]
rownames(accessory_essential) = row_names_accessory
row_names_accessory = row_names_accessory[-1 ,]

core_essential = core_essential[-1 ,]
accessory_essential = accessory_essential[-1 ,]

write.table(accessory_essential, file = '~/Papin Lab/Research Data/K. pneumoniae/Figures/essential_accessory_rxns.txt', sep="\t")

# separating accessory essential genes with a 55% threshold

accessory_yes = matrix(, ncol=ncol(accessory_essential))
accessory_no = matrix(, ncol=ncol(accessory_essential))
row_names_accessory_yes = matrix()
row_names_accessory_no = matrix()

in_vitro = accessory_essential[, 40:56]
in_vivo = accessory_essential[, 1:39]

for (i in 0:nrow(accessory_essential)) {
  if (34<sum(in_vivo[i ,]) & sum(in_vitro[i ,])<2.5) {
    accessory_yes = rbind(accessory_yes, accessory_essential[i ,])
    row_names_accessory_yes = rbind(row_names_accessory_yes, row_names_accessory[i])
    print(sum(in_vitro[i ,]))
  }
  if (5>sum(in_vivo[i ,]) & sum(in_vitro[i ,])>14) {
    accessory_yes = rbind(accessory_yes, accessory_essential[i ,])
    row_names_accessory_yes = rbind(row_names_accessory_yes, row_names_accessory[i])
    print(sum(in_vitro[i ,]))
  }
}

row_names_accessory_yes = row_names_accessory_yes[-1 ,]
row_names_accessory_no = row_names_accessory_no[-1 ,]

row_names_accessory_yes = c('L-valine transport via diffusion', 'L-valine reversible transport via proton symport')

accessory_yes = accessory_yes[-1 ,]
accessory_yes = accessory_yes[-1 ,]
#accessory_yes = accessory_yes[-13 ,]
#accessory_yes = accessory_yes[-10 ,]
#accessory_yes = accessory_yes[-3 ,]
accessory_no = accessory_no[-1 ,]

rownames(accessory_yes) = row_names_accessory_yes
rownames(accessory_no) = row_names_accessory_no

# the creation of the heatmap

library(pheatmap)

my_hclust_gene <- hclust(dist(essential), method = "complete")

# load package
library(dendextend)

as.dendrogram(my_hclust_gene) %>%
  plot(horiz = TRUE)

my_gene_col <- cutree(tree = as.dendrogram(my_hclust_gene), k = 2)

my_gene_col <- data.frame(cluster = ifelse(test = my_gene_col == 1, yes = "cluster 1", no = "cluster 2"))

head(my_gene_col)

col_clust = data.frame("T" = rep(c("Clinical", "Laboratory"), c(39,17)))
row.names(col_clust) = colnames(essential)

my_color = list(T = c(Clinical = "darkorchid2", Laboratory = "aquamarine2"))

png('~/Papin Lab/Research Data/K. pneumoniae/Figures/Essential_rxns.png', width=5000, height=10000, unit="px", res=150)

pheatmap(essential, cellwidth = 10, cellheight = 10, annotation_colors = my_color, show_colnames=FALSE, annotation_col = col_clust, cluster_cols = FALSE, legend=FALSE, fontsize_number = 0.9* fontsize)

dev.off()
