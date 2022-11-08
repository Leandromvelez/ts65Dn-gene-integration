# Set working directory (setwd) where preprocessed files are saved

# load required packages (lines for installing packages from Bioconductor included)
library(ggVennDiagram)
library(qgraph)
library(factoextra)
library(FactoMineR)
# if (!require("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
# BiocManager::install("qvalue")
library(qvalue)
library(reshape2)
# BiocManager::install("GO.db")
library('GO.db')
# BiocManager::install("impute")
library('impute')
# BiocManager::install("preprocessCore")
library('preprocessCore')
library(WGCNA)
library(ggplot2)
library(MetBrewer)
library(forcats)

# Load GTEx preprocessed data from working directory
load('GTEx NA included env.RData')
working_dataset=GTEx_subfiltered
row.names(working_dataset) = working_dataset$gene_tissue
working_dataset$gene_tissue=NULL
working_dataset = as.data.frame(t(working_dataset))
working_dataset[1:5,1:5]

# Load differential expression (DE) analysis data from mouse model (4 tissues)
liver = read.csv('mouse DEGS/liver degs.csv')
liver$tissue = paste0('liver')
gwat = read.csv('mouse DEGS/GWAT degs.csv')
gwat$tissue = paste0('gwat')
iwat = read.csv('mouse DEGS/iwat degs.csv')
iwat$tissue = paste0('iwat')
muscle = read.csv('mouse DEGS/sk musc degs.csv')
muscle$tissue = paste0('muscle')

# make a list of significant genes (p < 0.01) from 4 tissues 
gene_list <- list(Liver= liver$Symbol[liver$Liver.P.Value<0.01],
                  Iwat = iwat$Symbol[iwat$iWAT.P.Value<0.01],
                  Gwat = gwat$Symbol[gwat$gWAT.P.Value<0.01],
                  Muscle = muscle$Symbol[muscle$SkMusc.P.Value<0.01])
pdf(file = 'Venn Diagram Significant Genes across tissues P less 0.01.pdf')
# plot a venn diagram to show relationship between genes (cutoff p<0.01) in the 4 tissues
ggVennDiagram(gene_list, label_alpha = 0) + scale_fill_distiller(palette = "RdYlBu") + ggtitle('Tissue-sharing of DEGs P <0.01')
dev.off()

# make a list of significant genes (p < 0.05) from 4 tissues 
gene_list <- list(Liver= liver$Symbol[liver$Liver.P.Value<0.05],
                  Iwat = iwat$Symbol[iwat$iWAT.P.Value<0.05],
                  Gwat = gwat$Symbol[gwat$gWAT.P.Value<0.05],
                  Muscle = muscle$Symbol[muscle$SkMusc.P.Value<0.05])
pdf(file = 'Venn Diagram Significant Genes across tissues P less 0.05.pdf')
# plot a venn diagram to show relationship between genes (cutoff p<0.05) in the 4 tissues
ggVennDiagram(gene_list, label_alpha = 0) + scale_fill_distiller(palette = "RdYlBu") + ggtitle('Tissue-sharing of DEGs P <0.05')
dev.off()

# make a list of significant genes (p < 0.001) from 4 tissues 
gene_list <- list(Liver= liver$Symbol[liver$Liver.P.Value<0.001],
                  Iwat = iwat$Symbol[iwat$iWAT.P.Value<0.001],
                  Gwat = gwat$Symbol[gwat$gWAT.P.Value<0.001],
                  Muscle = muscle$Symbol[muscle$SkMusc.P.Value<0.001])
pdf(file = 'Venn Diagram Significant Genes across tissues P less 0.001.pdf')
# plot a venn diagram to show relationship between genes (cutoff p<0.001) in the 4 tissues
ggVennDiagram(gene_list, label_alpha = 0) + scale_fill_distiller(palette = "RdYlBu") + ggtitle('Tissue-sharing of DEGs P <0.001')
dev.off()

# extract DE genes with p<0.001 from 4 tissues, and combine them into a dataframe
df1 = as.data.frame(liver$Symbol[liver$Liver.P.Value<0.001])
colnames(df1) = 'gene_symbol'
df1$tissue = paste0('Liver')

df2 = as.data.frame(iwat$Symbol[iwat$iWAT.P.Value<0.001])
colnames(df2) = 'gene_symbol'
df2$tissue = paste0('iwat')

df3 = as.data.frame(gwat$Symbol[gwat$gWAT.P.Value<0.001])
colnames(df3) = 'gene_symbol'
df3$tissue = paste0('gwat')

df4 = as.data.frame(muscle$Symbol[muscle$SkMusc.P.Value<0.001])
colnames(df4) = 'gene_symbol'
df4$tissue = paste0('muscle')

full_sig_degs = as.data.frame(rbind(df1, df2, df3, df4))
full_sig_degs = na.omit(full_sig_degs)
table(full_sig_degs$tissue)

# read mouse gene info containing human orthologues
orth_table = read.delim('Mouse Gene info with Human Orthologues.txt')
# keep human orthologues of DE genes with p<0.001 from 4 tissues (liver, iwat, gwat, muscle)
full_sig_degs$human_orth = orth_table$human_orth[match(full_sig_degs$gene_symbol, orth_table$Symbol)]
# Change names of tissues (i.e. iwat to Adipose - Subcutaneous) to make integration with GTEx human data
full_sig_degs$human_tissue = ifelse(full_sig_degs$tissue=='Liver', 'Liver', '')
full_sig_degs$human_tissue = ifelse(full_sig_degs$tissue=='iwat', 'Adipose - Subcutaneous', paste0(full_sig_degs$human_tissue))
full_sig_degs$human_tissue = ifelse(full_sig_degs$tissue=='gwat', 'Adipose - Visceral (Omentum)', paste0(full_sig_degs$human_tissue))
full_sig_degs$human_tissue = ifelse(full_sig_degs$tissue=='muscle', 'Muscle - Skeletal', paste0(full_sig_degs$human_tissue))
full_sig_degs = na.omit(full_sig_degs)
full_sig_degs$gene_tissueH = paste0(full_sig_degs$human_orth, '_', full_sig_degs$human_tissue)
gene_set = full_sig_degs$gene_tissueH
# extract human orthologues of DE genes from GTEx working dataset into a dataframe
isogenes = working_dataset[,colnames(working_dataset) %in% gene_set]
ii = na.omit(isogenes)

## get bicor and p values from multiple biweight midcorrelations (all human orthologues of DE genes, and 4 tissues) in GTEx working dataset
cc1 = bicorAndPvalue(ii, ii, use = 'p')
cc3 = cc1$bicor
cc3[is.na(cc3)] = 0
cc4 = ifelse(cc1$p<0.01, '*', '')

colnames(ii)[1:10]

## data preparation for heatmap correlation structure
anno = data.frame(row.names(cc3), Group=gsub(".*_","", row.names(cc3)))
row.names(anno) = row.names(cc3)
library(RColorBrewer)
anno$row.names.cc3.=NULL
pdf(file = 'global cor structure of DEGS.pdf')
breaksList = seq(-1, 1, by = .1)
# heatmap correlation structure of human orthologues from DE genes (4 tissues) in GTEx expression dataset.
pheatmap::pheatmap(cc3, annotation_row = anno, display_numbers = F, labels_row = F, fontsize_number = 0, labels_col = F, color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdGy")))(length(breaksList)) )
dev.off()

#extract bicor and p values into a dataframe
cors_df = reshape2::melt(as.matrix(cc1$bicor))  
colnames(cors_df) = c('gene_tissue_1', 'gene_tissue_2', 'bicor')
nn1 = reshape2::melt(as.matrix(cc1$p))  
cors_df$pvalue = nn1$value
#add q value (FDR)
qq1 = qvalue::qvalue(cors_df$pvalue)
cors_df$qvalue = qq1$qvalues

cors_df = cors_df[order(cors_df$qvalue, decreasing = F),]
cors_df$gene_symbol_1 = gsub("\\_.*","",cors_df$gene_tissue_1)
cors_df$gene_symbol_2 = gsub("\\_.*","",cors_df$gene_tissue_2)
cors_df$tissue_1 = gsub(".*_","",cors_df$gene_tissue_1)
cors_df$tissue_2 = gsub(".*_","",cors_df$gene_tissue_2)
head(cors_df)
# make a dataframe excluding self correlated genes (bicor = 1)
tt1= cors_df[!cors_df$gene_tissue_1==cors_df$gene_tissue_2,]

# group by origin tissue (tissue_1) gene and get average bicor
tt2 = tt1 %>% dplyr::select(gene_tissue_1, bicor) %>% dplyr::group_by(gene_tissue_1) %>% dplyr::summarise(mean_abs=mean(abs(bicor), na.rm=T))
# order by decreasing mean bicor value
top_genes = tt2[order(tt2$mean_abs, decreasing = T),]
# extract top 25 genes
topgset = as.vector(top_genes$gene_tissue_1[1:25])
# extract the top 25 genes from the bicor and p value dataframe
nn1 = tt1[tt1$gene_tissue_1 %in% topgset,]

# save top 25 genes
write.csv(nn1, file = 'gene enrichments with cumulative DEGs.csv', row.names = F)

############# Genetic correlation Network

### get genes associated with major pathways in adipose and liver (Fig. 10C)
# read uniprot human database with gene names and go terms associated
pathways_annots = read.delim('uniprot-human-genes and goterms mapping.tab')
# extract genes associated with the top pathway from centrality estimates: adipose ('endoplasmic reticulum') and liver ('interferon-alpha')
pathw_set_adip = pathways_annots[grepl('endoplasmic reticulum', pathways_annots$Gene.ontology..biological.process.),]
pathw_set_adip$gene_tissue = paste0(pathw_set_adip$Gene.names...primary.., '_', 'Adipose - Subcutaneous')
pathw_set_liv = pathways_annots[grepl('interferon-alpha', pathways_annots$Gene.ontology..biological.process.),]
pathw_set_liv$gene_tissue = paste0(pathw_set_liv$Gene.names...primary.., '_', 'Liver')
# combine top 200 genes from liver and adipose into a dataframe
full_paths = c(pathw_set_adip$gene_tissue[1:200], pathw_set_liv$gene_tissue[1:200])
# use the top genes to filter 
tt3 = tt1[tt1$gene_tissue_2 %in% full_paths,]
tt3 = na.omit(tt3)

# read human secreted proteins table
sec_prots = read.delim('human secreted proteins.tab')
# keep only secreted proteins from origin tissue (gene_symbol_1)
tt3 = tt3[tt3$gene_symbol_1 %in% sec_prots$Gene.names...primary..,]
# group by gene_tissue1 (origin gene tissue) and get average bicor value
tt2 = tt3 %>% dplyr::select(gene_tissue_1, bicor) %>% dplyr::group_by(gene_tissue_1) %>% dplyr::summarise(mean_abs=mean(abs(bicor), na.rm=T))
# order by decreasing mean bicor value
top_genes = tt2[order(tt2$mean_abs, decreasing = T),]
# extract top 20 genes
topgset = as.vector(top_genes$gene_tissue_1[1:20])
# extract the top 20 genes from the bicor and p value dataframe
nn1 = tt1[tt1$gene_tissue_1 %in% topgset,]
# add column with -log(qvalue)
nn1$logq = -log(nn1$qvalue)
## add column with mean bicor values
nn1$meanbics = tt2$mean_abs[match(nn1$gene_tissue_1, tt2$gene_tissue_1)]
nn1 = na.omit(nn1)
length(unique(nn1$gene_tissue_1))
pdf(file = 'top-ranked candidates mean connectivity.pdf')
## top ranked gene connectivity
ggplot(nn1, aes(x=fct_reorder(gene_tissue_1, meanbics, .desc = TRUE), y=logq, fill=tissue_1)) + theme_classic() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))  + geom_violin(width=0.6) + geom_boxplot(width=0.2, position = position_dodge(width=0.6), alpha=0.3, color='grey') + xlab('') + ylab('gene-gene connectivity -log(qvalue)')
dev.off()

# extract top ranked from nn1 (gene connectivity) plus genes associated with the top pathway from centrality estimates: adipose ('endoplasmic reticulum') and liver ('interferon-alpha')
set1 = c(paste(nn1$gene_tissue_1[1]), full_paths)
## extract genes from GTEx working dataset
all_tog = working_dataset[,colnames(working_dataset) %in% set1]
## get bicor and p values from multiple biweight midcorrelations in GTEx working dataset
cc1 = bicor(all_tog, all_tog, use = 'p')
map1 = as.data.frame(cc1)
map1 = reshape2::melt(as.matrix(map1))
map1$tissue1 = gsub(".*_","",map1$Var1)
map1$tissue2 = gsub(".*_","",map1$Var2)
head(map1)
# highlight ATRN from adipose for network plot
map1$tissue_col = ifelse(grepl('ATRN', map1$Var1), 'seagreen1', 'darkorange2')
# highlight liver and adipose tissue for network plot
map1$tissue_col = ifelse(grepl('Liver', map1$tissue1), 'darkorchid2', paste0(map1$tissue_col))

map1$value = as.numeric(map1$value)
map1$value[map1$value > 0.999999] <- 0
map2 = reshape2::dcast(map1, Var1 ~ Var2, value.var = 'value')
row.names(map2) = map2$Var1
map2$Var1 = NULL
table(map1$tissue1)
colkey1 = colnames(map2)
names(colkey1) = map1$tissue_col[match(colkey1, map1$Var1)]
## undirected network, highlighting centrality of ATRN gene
pdf(file = paste('Undirected network qvalLess 1e-3 - ', "ATRN", '.pdf'))
qgraph(map2, minimum = 0.25, cut = 0.6, vsize = 2, color=names(colkey1), legend = F, borders = TRUE, layout='spring', posCol = "dodgerblue3", negCol = "firebrick3", label.cex=F, directed=F, labels = F) + ggtitle('')
dev.off()
## now repeat the undirected network, adding gene names
pdf(file = paste('Undirected network qvalLess 1e-3 - ', "ATRN", ' with labels.pdf'))
qgraph(map2, minimum = 0.2, cut = 0.6, vsize = 2, color=names(colkey1), legend = F, borders = TRUE, layout='spring', posCol = "dodgerblue3", negCol = "firebrick3", label.cex=3, directed=F, labels = gsub("\\_.*","",colnames(map2))) + ggtitle('')
dev.off()
