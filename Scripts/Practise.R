library('ggplot2')
library('pheatmap')
library('vegan')
library('compositions')
library('ggrepel')
library('phyloseq')
library('tidyverse')

### Pour changer le dendogram et l'ordre des lignes (ou colonnes), utiliser fix(pheatmap) et changer la ligne 70 :
# trees_row cluster_mat(mat, distance = clustering_distance_rows, 
# method = clustering_method)
# par
# trees_row = clust


fix(pheatmap)

codes_bioinfo <- read.table(file = "reference/codes_bioinfo_plantes.txt",header = T,sep="\t")
shared <- read.table("stats/shared_new_01_clean.tsv",header =T, check.names=F)
shared <- shared[,-3]
shared <- shared[,c(2:length(colnames(shared)))]

rownames(shared) = shared[,1]
shared <- shared[,c(2:length(colnames(shared)))]

shared[shared == 0] <- 0.1

shared_clr <- clr(shared)

distance <- vegdist(shared_clr,method = "euclidian")
clust <- hclust(distance)

distance <- as.matrix(distance)


codes_bioinfo

pr <- prcomp(distance)
comp <- as.data.frame(pr$rotation)
comp$Code_Miseq <- rownames(comp)

comp <- merge(x = comp, y = codes_bioinfo, by = "Code_Miseq", all.y = TRUE)


ggplot(comp,aes(PC1,PC2,color = Plante)) + geom_point(aes(shape=Site),size=3)+
  scale_colour_manual(values=c("#80b1d3","#b3de69","#9F7CFC","#D5CE00","#ff4f4f","#fdb462")) +
   xlab("PC1 (26.8%)") + ylab('PC2 (22.7%)')

my_colour = list(Site = c(Kaw="#f0f0f0",Nouragues="#bdbdbd",Petit_Saut="#636363"),
                 Plante = c(Aechmea_aquilega = "#80b1d3",CF_Aechmea_mertensii = "#b3de69",Catopsis_berteroniana = "#9F7CFC",Heliconia = "#DB7093", Vriesea_pleistoca ="#ff4f4f",Vriesea_splendens="#fdb462"))

#%pca
fviz_pca_ind(pr)






map$tree_row$order
codes_bioinfo$Code.bioinfo

## having rownames in good order
map$tree_row$labels[map$tree_row$order]

ordering_rows <- map$tree_row$labels[map$tree_row$order]

# reordering labels of bioinfo codes
col_order <- rownames(shared[clust$order,])

codes_bioinfo[order(factor(codes_bioinfo$Code_Miseq ,levels = order )),]

codes_bioinfo = codes_bioinfo[order(factor(codes_bioinfo$Code_Miseq,levels =  col_order )),]


rownames(annotation_row)


annotation_row = data.frame(
  Plante = codes_bioinfo$Esp..ce,
  Site = codes_bioinfo$Site)

my_colour = list(Site = c(Kaw="#f0f0f0",Nouragues="#bdbdbd",Petit_Saut="#636363"),
                 Plante = c(Aechmea_aquilega = "#80b1d3",CF_Aechmea_mertensii = "#b3de69",Catopsis_berteroniana = "#9F7CFC",Heliconia = "#D5CE00", Vriesea_pleistoca ="#ff4f4f",Vriesea_splendens="#fdb462"))


rownames(annotation_row) = col_order

# Fichier avec le numéro des OTUs à afficher dans l'ordre.
row_order = read.table('HEATMAPS/ORDER_OTUS')
row_order <- row_order$V1
shared <- shared[,row_order]


pheatmap(shared,cutree_rows = 13, border_color = NA, annotation_row = annotation_row,cluster_cols = F, clustering_method = "complete",annotation_colors = my_colour)













# TEEST

df_shared_clr <- as.data.frame(shared_clr)


ggplot(df_shared_clr,aes(df_shared_clr[,1])) + 











# pca %


bray_curtis <- vegdist(shared,diag = F)







shared_2 <- read.table("stats/shared_new_01_clean.tsv",header =T, check.names=F)

shared_2 <- shared_2[,-3]
shared_2 <- shared_2[,c(2:length(colnames(shared_2)))]

rownames(shared_2) = shared_2[,1]
shared_2 <- shared_2[,c(2:length(colnames(shared_2)))]

shared_2$Code.Miseq <- rownames(shared_2)

cooodes <- read.table('reference/codes_bioinfo_plantes.txt',header=T,sep="\t")

merging <- merge(x=shared_2,y=cooodes,by="Code.Miseq")

rownames(merging) <- merging$group

random <- rep(c(1,2),len=92)
merging$random <- random

matrix <- merging[, 2:90]

esp <- merging$Esp..ce

site <- merging$Site

data = merging[,c('Esp..ce','Site')]

adonis2(distance ~ esp+site ,data= data, permutations = 9999,by="margin")



anosim(distance, esp, distance = "bray", permutations = 9999)














bray_curtis

# coordinates

coordinates <- read.table('reference/coordinates',header = T,sep="\t",row.names = 1)

coordinates_dist <- stats::dist(coordinates)
plot(hclust(coordinates_dist))

coordinates_dist <- as.matrix(coordinates_dist) 
for (i in rownames(coordinates_dist)) {
  print(i)
}
shared_2$Co
