#############################################################
#Name - Architha Ellappalayam
#Project - Proteomics Analysis on Breast Cancer Data
#Script - Applying SNFTool on Krug et.al. protegenomics data  
#Version - 0.1 
#############################################################

#Set the working directory
setwd("M:/AppData/FolderRedirection/Desktop/fourth_project/Breast_Cancer_Proteomics/")

#Downloading the required packages 
if (!requireNamespace("BiocManager", quietly=TRUE))
  install.packages("BiocManager")
BiocManager::install("readxl")
BiocManager::install("SNFtool")
BiocManager::install("gplots")
BiocManager::install("dendextend")
BiocManager::install("randomcoloR")
BiocManager::install("limma")
BiocManager::install("devtools")
BiocManager::install("pathfindR")
BiocManager::install("clusterProfiler")
BiocManager::install("ggplot2")
BiocManager::install("corrplot")
BiocManager::install('EnhancedVolcano')
BiocManager::install("ReactomePA")
BiocManager::install("org.Hs.eg.db")
BiocManager::install("pheatmap")

#Loading the required libraries 
library(BiocManager)
library(readxl)
library(SNFtool)
library(gplots)
library(matrixStats)
library(proxy)
library(dplyr)
library(dendextend)
library(randomcoloR)
library(limma)
library(devtools)
library(pathfindR)
library(org.Hs.eg.db)
library(ggplot2)
library(clusterProfiler)
library(corrplot)
library(ggrepel)
library(org.Hs.eg.db)
library(pheatmap)

################### Krug et.al.#######################

#Loading the proteogenomis data matrix from Krug et.al (2020) 
#Metadata from Krug et al. 
metadata_krug <- readxl::read_excel("Krug_et_al/Metadata_Supplemental_Table_1.xlsx", sheet = 2)

#Proteome data fom Krug et.al. 
proteome_krug <- readxl::read_excel("Krug_et_al/NIHMS1687926-supplement-Supplemental_Table_2 (1).xlsx", sheet = 3)
proteome_krug <- as.matrix(proteome_krug)
rownames(proteome_krug) <- proteome_krug[,2]

#Make a proteome matrix from the data above 
proteome_matrix_krug <- proteome_krug[,c(15:136)]
proteome_matrix_krug[is.na(proteome_matrix_krug)] <- 0
proteome_matrix_krug <- as.data.frame(proteome_matrix_krug)
dim(proteome_matrix_krug)
#10107 122

#Total of 122 samples present 

#Running the SNFTool on the Proteomic datase Krug et.al. 
#Setting the parameter for the SNF analysis 
K = 10;		# number of neighbors, usually (10~30)
alpha = 0.5;  	# hyperparameter, usually (0.3~0.8)
T = 20; 	# Number of Iterations, usually (10~20)

#Adding labels to the samples - PAm50 subtype 
truelabel_krug <- c(rep("Basal",29), rep("LumA", 57), rep("LumB", 17), rep("Her2", 14), rep("Normal-like", 5))

# #Make separate lists for each PAM50 data type m
samples_basal_krug <- metadata_krug$Sample.ID[metadata_krug$PAM50 == "Basal"]
samples_luma_krug <- metadata_krug$Sample.ID[metadata_krug$PAM50 == "LumA"]
samples_lumb_krug <- metadata_krug$Sample.ID[metadata_krug$PAM50 == "LumB"]
samples_normal_krug <- metadata_krug$Sample.ID[metadata_krug$PAM50 == "Normal-like"]
samples_her2_krug <- metadata_krug$Sample.ID[metadata_krug$PAM50 == "Her2"]

# Select the columns from data_matrix that match the samples in samples_to_select
subset_basal <- as.data.frame(proteome_matrix_krug[, intersect(colnames(proteome_matrix_krug), samples_basal_krug)])
subset_luma <- as.data.frame(proteome_matrix_krug[, intersect(colnames(proteome_matrix_krug), samples_luma_krug)])
subset_lumb <- as.data.frame(proteome_matrix_krug[, intersect(colnames(proteome_matrix_krug), samples_lumb_krug)])
subset_normal <- as.data.frame(proteome_matrix_krug[, intersect(colnames(proteome_matrix_krug), samples_normal_krug)])
subset_her2 <- as.data.frame(proteome_matrix_krug[, intersect(colnames(proteome_matrix_krug), samples_her2_krug)])

#Convert the dataframe value to numeric 
subset_basal <- mutate_all(subset_basal, function(x) as.numeric(as.character(x)))
subset_luma <- mutate_all(subset_luma, function(x) as.numeric(as.character(x)))
subset_lumb <- mutate_all(subset_lumb, function(x) as.numeric(as.character(x)))
subset_normal <- mutate_all(subset_normal, function(x) as.numeric(as.character(x)))
subset_her2 <- mutate_all(subset_her2, function(x) as.numeric(as.character(x)))

#Transposing the matrix 
subset_basal <- t(subset_basal)
subset_her2 <- t(subset_her2)
subset_luma <- t(subset_luma)
subset_lumb <- t(subset_lumb)
subset_normal <- t(subset_normal)

########### SNF tool for Basal and Luminal samples ####################

#Splitting the basal samples into two matrices 
# create a vector of row numbers to include in the first dataframe
indices_basal <- sample(ncol(subset_basal), floor(ncol(subset_basal)/2))

# split the dataframe based on the row indices
subset_basal_1 <- subset_basal[,indices_basal]
subset_basal_2 <- subset_basal[,-indices_basal]

#Splitting the basal samples into two matrices 
# create a vector of row numbers to include in the first dataframe
indices_luma <- sample(ncol(subset_luma), floor(ncol(subset_luma)/2))

# split the dataframe based on the row indices
subset_luma_1 <- subset_luma[,indices_luma]
subset_luma_2 <- subset_luma[,-indices_luma]

#Constructing the pair-wise distnace matrix 
basal_dist_matrix_1 = (dist2(as.matrix(subset_basal_1),as.matrix(subset_basal_1)))^(1/2)
basal_dist_matrix_2 = (dist2(as.matrix(subset_basal_2),as.matrix(subset_basal_2)))^(1/2)

luma_dist_matrix_1 = (dist2(as.matrix(subset_luma_1),as.matrix(subset_luma_1)))^(1/2)
luma_dist_matrix_2 = (dist2(as.matrix(subset_luma_2),as.matrix(subset_luma_2)))^(1/2)

#Constructing the similarity graphs
W_Basal_1 = affinityMatrix(basal_dist_matrix_1, K, alpha)
W_Basal_2 = affinityMatrix(basal_dist_matrix_2, K, alpha)

W_luma_1 = affinityMatrix(luma_dist_matrix_1, K, alpha)
W_luma_2 = affinityMatrix(luma_dist_matrix_2, K, alpha)


## then the overall matrix can be computed by similarity network fusion(SNF):
W_Basal = SNF(list(W_Basal_1,W_Basal_2), K, T)
W_Luma = SNF(list(W_luma_1,W_luma_2), K, T)


# create a dendrogram for hierarchical clustering
dend_basal <- hclust(as.dist(1-W_Basal))

# plot the heatmap with dendrogram
plot(dend_basal)
heatmap(W_Basal)

# create a dendrogram for hierarchical clustering
dend_luma <- hclust(as.dist(1-W_Luma))

# plot the heatmap with dendrogram
plot(dend_luma)
heatmap(W_Basal)

#Basal samples 
metadata_basal <- metadata_krug[metadata_krug$PAM50 == "Basal",]
er_bars <- ifelse(metadata_basal$ER.Updated.Clinical.Status == "positive", "firebrick3", "beige")
pr_bars <- ifelse(metadata_basal$PR.Clinical.Status== "positive", "red", "yellow")
nmf_bars <- ifelse(metadata_basal$NMF.Cluster == "Basal-I", "cyan", "green")
erbb2_bars <- ifelse(metadata_basal$ER.Updated.Clinical.Status == "negative", "blue", "yellow")
tnbc_bars <- ifelse(metadata_basal$TNBC.Updated.Clinical.Status == "positive", "orange", "violet")
bars <- cbind(er_bars, pr_bars, nmf_bars, erbb2_bars, tnbc_bars)
par(mar= c(13,2,6,2))
plot(dend_basal, xlab = NA, sub = NA, main = "Clustered dendrogram of Basal samples from Krug et.al.")
colored_bars(colors = bars, dend = dend_basal, rowLabels = c("ER", "PR", "NMF", "ERBB2", "TNBC"))
plot.new()
legend("topright",legend=c("ER Positive","ER negative", "PR positive", "PR negative", "NMF - Basal", "NMF - HER2", "ERBB2 negative", "ERBB2 unknown", "TNBC positive", "TNBC negative"),pch=15,col=c("firebrick3","beige", "red", "yellow", "cyan", "green", "blue", "yellow", "orange", "violet"), inset = c(-0.1,0))


#Luminal samples 
metadata_luma <- metadata_krug[metadata_krug$PAM50 == "LumA",]
er_bars_luma <- ifelse(metadata_luma$ER.Updated.Clinical.Status == "positive", "firebrick3", "beige")
pr_bars_luma <- ifelse(metadata_luma$PR.Clinical.Status== "positive", "red", "yellow")
nmf_bars_luma <- ifelse(metadata_luma$NMF.Cluster == "Basal-I", "cyan", "green")
erbb2_bars_luma <- ifelse(metadata_luma$ER.Updated.Clinical.Status == "negative", "blue", "yellow")
tnbc_bars_luma <- ifelse(metadata_luma$TNBC.Updated.Clinical.Status == "positive", "orange", "violet")
bars_luma <- cbind(er_bars_luma, pr_bars_luma, nmf_bars_luma, erbb2_bars_luma, tnbc_bars_luma)
plot(dend_luma)
colored_bars(colors = bars_luma, dend = dend_luma, rowLabels = c("ER", "PR", "NMF", "ERBB2", "TNBC"))
legend("top",legend=c("ER Positive","ER negative", "PR positive", "PR negative", "Basal", "HER2", "ERBB2 negative", "ERBB2 unknown", "TNBC positive", "TNBC negative"),pch=15,col=c("firebrick3","beige", "red", "yellow", "cyan", "green", "blue", "yellow", "orange", "violet"), inset = c(-0.1,0))
plot(0,0)

## You can display clusters in the data by the following function
## where C is the number of clusters.
C = 2 # number of clusters
group = spectralClustering(W_Basal,C); # the final subtypes information
group_luma = spectralClustering(W_Luma, C)

## Visualize the clusters present in the given similarity matrix
## as well as some sample information
## In this presentation no clustering method is ran the samples
## are ordered in function of their group label present in the group arguments
displayClustersWithHeatmap(W_Basal, group)
displayClustersWithHeatmap(W_Luma, group)



###################### Differential expression analysis - Krug et al ###################
#loading the excel sheet with two sample sets 
ds_samples_one_krug <- read.csv("genes_basal_one_krug.csv") #12 samples
ds_samples_two_krug <- read.csv("genes_basal_two_krug.csv") #17 samples 

ds_basal_one_krug <- subset_basal[rownames(subset_basal) %in% ds_samples_one_krug$SET_1,]
#12 10107
ds_basal_one_krug <- t(ds_basal_one_krug)
colnames(ds_basal_one_krug) <- rep("Basal_Cluster_One", 12)

ds_basal_two_krug <- subset_basal[rownames(subset_basal) %in% ds_samples_two_krug$SET_2,]
#17 10107
ds_basal_two_krug <- t(ds_basal_two_krug)
colnames(ds_basal_two_krug) <- rep("Basal_Cluster_Two", 17)

#Cmbining the matrices into one 
basal_log_norm_krug <- cbind(ds_basal_one_krug, ds_basal_two_krug)

#Design matrix
design_basal_set_krug <- as.factor(colnames(cbind(ds_basal_one_krug, ds_basal_two_krug)))
design_basal_set_krug


#Standard limma model fit
design.trainset.krug <- model.matrix(~0 + design_basal_set_krug)
colnames(design.trainset.krug) <- c("Basal_Cluster_One", "Basal_Cluster_Two")
fit.mydat <- lmFit(cbind(ds_basal_one_krug, ds_basal_two_krug), design.trainset.krug)

#Checking whether the two sample groups are the same or not 
contrasts.matrix.mydat.krug <- makeContrasts(Basal_Cluster_Two - Basal_Cluster_One, levels = design.trainset.krug)

#Fitting the data 
mydat_fits_krug <- contrasts.fit(fit.mydat, contrasts.matrix.mydat.krug)
mydat_ebfit_krug <- eBayes(mydat_fits_krug)

#Viewing the sorted gene list based on the logFoldChange 
genelist_mydat_krug <- topTable(mydat_ebfit_krug, number = 33000, sort.by = "logFC", adjust.method = "BH")

#PLotting the heatmap of the significant probes 
#total nmumber of unique significant probes - 63
sig_probes_krug <- genelist_mydat_krug[genelist_mydat_krug$logFC > 0.58 & genelist_mydat_krug$adj.P.Val < 0.05,]
sig_probes_one_krug <- genelist_mydat_krug[genelist_mydat_krug$logFC < -0.58 & genelist_mydat_krug$adj.P.Val < 0.05,]

##Collecting the expression data of the sginficiant probes
significant_probe_krug <- c(rownames(sig_probes_krug), rownames(sig_probes_one_krug))
significant_basal_expr_krug <- (basal_log_norm_krug[rownames(basal_log_norm_krug) %in% significant_probe_krug,])
significant_basal_expr_krug <-  avereps(significant_basal_expr_krug)
dim(significant_basal_expr_krug)
#1020 29

significant_basal_expr_up_krug <- significant_basal_expr_krug[rownames(significant_basal_expr_krug) %in% rownames(sig_probes_krug),]
significant_basal_expr_down_krug <- significant_basal_expr_krug[rownames(significant_basal_expr_krug) %in% rownames(sig_probes_one_krug),]

highGene <- rownames(significant_basal_expr_up_krug) #15 85
lowGene <- rownames(significant_basal_expr_down_krug) #467 85


#Making a volcano plot of the differentially expressed genes 
xCol = col2rgb(ifelse(rownames(genelist_mydat_krug) %in% highGene, "Red", ifelse(rownames(genelist_mydat_krug) %in% lowGene, "Blue",'Gray')))
plot(genelist_mydat_krug$logFC,main = "Volcano Plot of differntially expressed genes in Basal Cluster One",
     -log10(genelist_mydat_krug$adj.P.Val),cex.axis=0.8,
     col= rgb(xCol[1,],xCol[2,],xCol[3,],maxColorValue=255),
     pch=20,
     cex = 1,cex.lab = 1.3,
     ylab = '-log10(adj.P.Val )',
     xlab = 'logFC',cex.axis=0.7)
legend("topright", legend = c("Down-regulated genes (n = 181)", "Up-regulated genes (n = 839)"),fill = c("Blue", "Red"))


#Heatmap 
heatmap(significant_basal_expr_krug, col = bluered(700))
heatmap.2(significant_basal_expr_krug, trace = "none", col = bluered(1000), tracecol = NA)

#Gene Ontology and Pathway analysis using pathfindR package 
sig_probes_pathfindr <- sig_probes_krug[,-c(2,3,4,6)]
sig_probes_pathfindr <- cbind(rownames(sig_probes_pathfindr), sig_probes_pathfindr)

sig_probes_one_pathfindr <- sig_probes_one_krug[,-c(2,3,4,6)]
sig_probes_one_pathfindr <- cbind(rownames(sig_probes_one_pathfindr), sig_probes_one_pathfindr)


################ SNF Tool on Basal samples from Anurag et al #################
#Loading the proteogenomis data matrix from Anurag et.al (2020) 
#Metadata from Anurag et al. 
metadata_anurag <- readxl::read_excel("Anurag_et_al/metadata.xlsx", sheet = 2)

#Proteome data fom Krug et.al. 
proteome_anurag <- readxl::read_excel("Anurag_et_al/proteome_anurag.xlsx", sheet = 3)
proteome_anurag <- as.matrix(proteome_anurag)
rownames(proteome_anurag) <- proteome_anurag[,1]

#Make a proteome matrix from the data above 
proteome_matrix_anurag <- proteome_anurag[,c(36:106)]
proteome_matrix_anurag[is.na(proteome_matrix_anurag)] <- 0
proteome_matrix_anurag <- as.data.frame(proteome_matrix_anurag)
dim(proteome_matrix_anurag)
#11063 71

#Filter the metadata based on sample names from proteome matrix 
metadata_anurag <- metadata_anurag[metadata_anurag$sample_id %in% colnames(proteome_matrix_anurag),]
dim(metadata_anurag)
#71 49

#Select only the Baseline samples 
metadata_anurag <- metadata_anurag[metadata_anurag$Collection_Event == "Baseline",]

K = 10;		# number of neighbors, usually (10~30)
alpha = 0.5;  	# hyperparameter, usually (0.3~0.8)
T = 20; 	# Number of Iterations, usually (10~20)

#Selecting samples which are Basal 
samples_basal_anurag <- metadata_anurag$sample_id

# Select the columns from data_matrix that match the samples in samples_to_select
subset_basal_anurag <- as.data.frame(proteome_matrix_anurag[, intersect(colnames(proteome_matrix_anurag), samples_basal_anurag)])

#Convert the dataframe value to numeric 
subset_basal_anurag <- mutate_all(subset_basal_anurag, function(x) as.numeric(as.character(x)))
subset_basal_anurag[is.na(subset_basal_anurag)] <- 0

#Transpose
subset_basal_anurag <- t(subset_basal_anurag)


#Splitting the basal samples into two matrices 
# create a vector of row numbers to include in the first dataframe
indices_basal_anurag <- sample(ncol(subset_basal_anurag), floor(ncol(subset_basal_anurag)/2))

# split the dataframe based on the row indices
subset_basal_1_anurag <- subset_basal_anurag[,indices_basal_anurag]
subset_basal_2_anurag <- subset_basal_anurag[,-indices_basal_anurag]

#Constructing the pair-wise distnace matrix 
basal_dist_matrix_1_anurag = (dist2(as.matrix(subset_basal_1_anurag),as.matrix(subset_basal_1_anurag)))^(1/2)
basal_dist_matrix_2_anurag = (dist2(as.matrix(subset_basal_2_anurag),as.matrix(subset_basal_2_anurag)))^(1/2)

#Constructing the similarity graphs
W_Basal_1_anurag = affinityMatrix(basal_dist_matrix_1_anurag, K, alpha)
W_Basal_2_anurag = affinityMatrix(basal_dist_matrix_2_anurag, K, alpha)

## then the overall matrix can be computed by similarity network fusion(SNF):
W_Basal_anurag = SNF(list(W_Basal_1_anurag,W_Basal_2_anurag), K, T)

## You can display clusters in the data by the following function
## where C is the number of clusters.
C = 2 # number of clusters
group_anurag = spectralClustering(W_Basal_anurag,C);

## In this presentation no clustering method is ran the samples
## are ordered in function of their group label present in the group arguments
displayClustersWithHeatmap(W_Basal_anurag, group_anurag)

# create a dendrogram for hierarchical clustering
dend_basal_anurag <- hclust(as.dist(1-W_Basal_anurag))

# plot the heatmap with dendrogram
plot(dend_basal_anurag)


#Adding colored bars to dendrograms 
pcr_anurag <- ifelse(metadata_anurag$pCR_response == "Yes", "green", "red")

#RCB - Residual Cancer Burden
rcb_anurag <- factor(metadata_anurag$RCB)
n_rcb_anurag <- length(unique(rcb_anurag))
color_rcb <- c("yellow","cyan","blue", "orange","pink")
col_rcb_anurag <- color_rcb[rcb_anurag]

#PAM50 subtype 
pam50_anurag <- factor(metadata_anurag$PAM50)
n_pam50_anurag <- length(unique(pam50_anurag))
color_pam50 <- c("#d239dd", "#f760ba", "#ebbef7", "#88dbef", "#60eaba", "#c253e8")
col_pam50_anurag <- color_pam50[pam50_anurag]

#TNBC type
tnbc_anurag <- factor(metadata_anurag$TNBCType)
n_tnbc_anurag <- length(unique(tnbc_anurag))
color_tnbc <- c("#dbc63f" ,"#14e88c" ,"#c6497b" ,"#b24513", "#5be549", "#3bb2d3", "#c7d7fc", "#9914ff")
col_tnbc_anurag <- color_tnbc[tnbc_anurag]

#Race of the patient
race_anurag <- factor(metadata_anurag$Race)
n_race_anurag <- length(unique(race_anurag))
color_race <- c("#c69201", "#dd2c7f", "#fff8bf", "#866fc6", "#f9a55c", "#93171b")
col_race_anurag <- color_race[race_anurag]

#randomColor(count = 2, hue = "monochrome")
par(mar=c(9,6,2,0))
par(mar= c(13,2,6,2))
plot(dend_basal_anurag, sub = NA, xlab = NA)
colored_bars(colors = cbind(col_tnbc_anurag,col_race_anurag,col_pam50_anurag,col_rcb_anurag, pcr_anurag), dend = dend_basal_anurag, rowLabels = c("TNBC","Race"  ,"PAM50", "RCB","pCR"))


plot(1, 1, axes = "none")
legend(x = 0.9, y = 1.4, legend = c("pCR", "no-pCR"), fill = c("green", "red"), text.width = 0.13)
legend(x = 0.9, y = 1.26, legend = c("RCB-0", "RCB-I", "RCB-II", "RCB-III"), fill = c("yellow","cyan","blue", "orange","pink"), text.width = 0.13)
legend(x = 0.9, y = 1.02, legend = c("Basal","HER2", "LumA", "LumB","Normal" ), fill = c("#d239dd", "#f760ba", "#ebbef7", "#88dbef", "#60eaba", "#c253e8"), text.width = 0.13)
legend(x = 0.9, y = 0.75, legend = c("AfricanAmerican","Caucasian", "Other"), fill = c("#c69201", "#dd2c7f", "#fff8bf", "#866fc6", "#f9a55c", "#93171b"), text.width = 0.13)
legend(x = 0.9, y = 0.55, legend = c("BL1", "BL2", "IM", "LAR", "M", "MSL", "NA", "UNS"), fill = c("#dbc63f" ,"#14e88c" ,"#c6497b" ,"#b24513", "#5be549", "#3bb2d3", "#c7d7fc", "#9914ff"), text.width = 0.13)


plotAlluvial(W_Basal_anurag, 1:5, col = "blue")

#save.image("Anurag_Krug_data.RData")

#################### Differential expression analysis######################

#loading the excel sheet with two sample sets 
ds_samples_one_anurag <- read.csv("gene_clusters_anurag.csv") #28 samples
ds_samples_two_anurag <- read.csv("set_2_anurag.csv") #27 samples 

ds_basal_one_anurag <- subset_basal_anurag[rownames(subset_basal_anurag) %in% ds_samples_one_anurag$SET_1,]
#28 11063
ds_basal_one_anurag <- t(ds_basal_one_anurag)
colnames(ds_basal_one_anurag) <- rep("Basal_Cluster_One", 28)

ds_basal_two_anurag <- subset_basal_anurag[rownames(subset_basal_anurag) %in% ds_samples_two_anurag$SET_2,]
#27 11063
ds_basal_two_anurag <- t(ds_basal_two_anurag)
colnames(ds_basal_two_anurag) <- rep("Basal_Cluster_Two", 27)

#Cmbining the matrices into one 
basal_log_norm_anurag <- cbind(ds_basal_one_anurag, ds_basal_two_anurag)

#Design matrix
design_basal_set_anurag <- as.factor(colnames(cbind(ds_basal_one_anurag, ds_basal_two_anurag)))
design_basal_set_anurag


#Standard limma model fit
design.trainset.anurag <- model.matrix(~0 + design_basal_set_anurag)
colnames(design.trainset.anurag) <- c("Basal_Cluster_One", "Basal_Cluster_Two")
fit.mydat <- lmFit(cbind(ds_basal_one_anurag, ds_basal_two_anurag), design.trainset.anurag)

#Checking whether the two sample groups are the same or not 
contrasts.matrix.mydat.anurag <- makeContrasts(Basal_Cluster_Two - Basal_Cluster_One, levels = design.trainset.anurag)

#Fitting the data 
mydat_fits_anurag <- contrasts.fit(fit.mydat, contrasts.matrix.mydat.anurag)
mydat_ebfit_anurag <- eBayes(mydat_fits_anurag)

#Viewing the sorted gene list based on the logFoldChange 
genelist_mydat_anurag <- topTable(mydat_ebfit_anurag, number = 33000, sort.by = "logFC", adjust.method = "BH")

#PLotting the heatmap of the significant probes 
#total nmumber of unique significant probes - 63
sig_probes_anurag <- genelist_mydat_anurag[genelist_mydat_anurag$logFC > 0.58 & genelist_mydat_anurag$adj.P.Val < 0.05,]
sig_probes_one_anurag <- genelist_mydat_anurag[genelist_mydat_anurag$logFC < -0.58 & genelist_mydat_anurag$adj.P.Val < 0.05,]

##Collecting the expression data of the sginficiant probes
significant_probe_anurag <- c(rownames(sig_probes_anurag), rownames(sig_probes_one_anurag))
significant_basal_expr_anurag <- (basal_log_norm_anurag[rownames(basal_log_norm_anurag) %in% significant_probe_anurag,])
significant_basal_expr_anurag <-  avereps(significant_basal_expr_anurag)
dim(significant_basal_expr_anurag)
#308 55

significant_basal_expr_up_anurag <- significant_basal_expr_anurag[rownames(significant_basal_expr_anurag) %in% rownames(sig_probes_anurag),]
significant_basal_expr_down_anurag <- significant_basal_expr_anurag[rownames(significant_basal_expr_anurag) %in% rownames(sig_probes_one_anurag),]

highGene <- rownames(significant_basal_expr_up_anurag) #4 85
lowGene <- rownames(significant_basal_expr_down_anurag) #59 85

ggplot(data=genelist_mydat_anurag, aes(x=logFC, y= -log10(adj.P.Val), col = diffexpressed, label = delabel)) + geom_point() + theme_bw() + scale_color_manual(values=c("blue", "grey", "red")) + geom_text() + geom_text_repel(max.overlaps = 0) 
ggplot(data=genelist_mydat_anurag, aes(x=logFC, y= -log10(adj.P.Val), col = diffexpressed, label = delabel)) + geom_point() + theme_bw() + scale_color_manual(values=c("blue", "grey", "red")) + geom_text() + geom_text_repel(max.overlaps = 0) + theme(axis.title = element_text(size = 15))  + theme(legend.text = element_text(size = 15)) + theme(legend.title = element_text(size = 15))


genelist_mydat_anurag$genes <- rownames(genelist_mydat_anurag)
genelist_mydat_anurag$diffexpressed <- "NO"
# if log2Foldchange > 0.6 and pvalue < 0.05, set as "UP" 
genelist_mydat_anurag$diffexpressed[genelist_mydat_anurag$logFC > 0.58 & genelist_mydat_anurag$adj.P.Val < 0.05] <- "UP"
# if log2Foldchange < -0.6 and pvalue < 0.05, set as "DOWN"
genelist_mydat_anurag$diffexpressed[genelist_mydat_anurag$logFC < -0.58 & genelist_mydat_anurag$adj.P.Val < 0.05] <- "DOWN"

genelist_mydat_anurag$delabel <- NA
genelist_mydat_anurag$delabel[genelist_mydat_anurag$diffexpressed != "NO"] <- genelist_mydat_anurag$genes[genelist_mydat_anurag$diffexpressed != "NO"]



#Making a volcano plot of the differentially expressed genes 
xCol = col2rgb(ifelse(rownames(genelist_mydat_anurag) %in% highGene, "Red", ifelse(rownames(genelist_mydat_anurag) %in% lowGene, "Blue",'Gray')))
plot(genelist_mydat_anurag$logFC,main = "Volcano Plot of differntially expressed genes in Basal samples - Anurag et.al.",
     -log10(genelist_mydat_anurag$adj.P.Val),cex.axis=0.8,
     col= rgb(xCol[1,],xCol[2,],xCol[3,],maxColorValue=255),
     pch=20,
     cex = 1,cex.lab = 1.3,
     ylab = '-log10(adj.P.Val )',
     xlab = 'logFC',cex.axis=0.7)
legend("topright", legend = c("Down-regulated genes (n = 2082)", "Up-regulated genes (n = 2210)"),fill = c("Blue", "Red"))


#Heatmap 
heatmap(significant_basal_expr_anurag, col = redgreen(100))
heatmap.2(significant_basal_expr_anurag, trace = "none", col = redgreen(100))

#Gene Ontology and Pathway analysis using pathfindR package 
sig_probes_pathfindr <- sig_probes_anurag[,-c(2,3,4,6)]
sig_probes_pathfindr <- cbind(rownames(sig_probes_pathfindr), sig_probes_pathfindr)
  
sig_probes_one_pathfindr <- sig_probes_one_anurag[,-c(2,3,4,6)]
sig_probes_one_pathfindr <- cbind(rownames(sig_probes_one_pathfindr), sig_probes_one_pathfindr)
   
output_df_basal_up_anurag <- pathfindR::run_pathfindR(sig_probes_pathfindr, gene_sets = "GO-BP")
enrichment_chart(output_df, top_terms = 20)
output_cluster <- cluster_enriched_terms(output_df)
enrichment_chart(output_cluster, plot_by_cluster = TRUE, top_terms = 5)

term_gene_heatmap(output_df)
UpSet_plot(output_df)
term_gene_graph(output_df)

############### Making volcano plots prettier ################
ggplot(data=genelist_mydat_krug, aes(x=logFC, y= -log10(adj.P.Val), col = diffexpressed, label = delabel)) + geom_point() + theme_bw() + scale_color_manual(values=c("blue", "grey", "red")) + geom_text() + geom_text_repel(max.overlaps = 0) 
ggplot(data=genelist_mydat_krug, aes(x=logFC, y= -log10(adj.P.Val), col = diffexpressed, label = delabel)) + geom_point() + theme_bw() + scale_color_manual(values=c("blue", "grey", "red")) + geom_text() + geom_text_repel(max.overlaps = 0) + theme(axis.title = element_text(size = 15))  + theme(legend.text = element_text(size = 15)) + theme(legend.title = element_text(size = 15))


genelist_mydat_krug$genes <- rownames(genelist_mydat_krug)
genelist_mydat_krug$diffexpressed <- "NO"
# if log2Foldchange > 0.6 and pvalue < 0.05, set as "UP" 
genelist_mydat_krug$diffexpressed[genelist_mydat_krug$logFC > 0.58 & genelist_mydat_krug$adj.P.Val < 0.05] <- "UP"
# if log2Foldchange < -0.6 and pvalue < 0.05, set as "DOWN"
genelist_mydat_krug$diffexpressed[genelist_mydat_krug$logFC < -0.58 & genelist_mydat_krug$adj.P.Val < 0.05] <- "DOWN"

genelist_mydat_krug$delabel <- NA
genelist_mydat_krug$delabel[genelist_mydat_krug$diffexpressed != "NO"] <- genelist_mydat_krug$genes[genelist_mydat_krug$diffexpressed != "NO"]

#Heatmap legend 
library(circlize)
col_fun = colorRamp2(c(0, 0.5, 1), c("khaki1", "beige", "red2"))
lgd = Legend(col_fun = col_fun, title = "Correlation")
plot.new()

pushViewport(viewport(width = 0.9, height = 0.9))
grid.rect()  # border
draw(lgd, x = unit(1, "cm"), y = unit(1, "cm"), just = c("left", "bottom"))
draw(lgd, x = unit(0.5, "npc"), y = unit(0.5, "npc"))
draw(lgd, x = unit(1, "npc"), y = unit(1, "npc"), just = c("right", "top"))
popViewport()


Genes2 <- select(org.Hs.eg.db, lowGene, 'ENTREZID', 'SYMBOL')
x <- enrichPathway(gene=Genes2$ENTREZID,pvalueCutoff=0.05, readable=T)

barplot(x, showCategory=25)

#Saving the data matrices as an RData object 
#save.image("proteomics_data.RData")
################## Clustering based on genes - Proteomics ############

#Standard normalize these matrices 
#subset_basal_norm <- standardNormalization(subset_basal)
#subset_luma_norm <- standardNormalization(subset_luma)
#subset_lumb_norm <- standardNormalization(subset_lumb)
#subset_normal_norm <- standardNormalization(subset_normal)
#subset_her2_norm <- standardNormalization(subset_her2)

#Calculate the pair-wise distance 
basal_dist_matrix = (dist2(as.matrix(subset_basal),as.matrix(subset_basal)))^(1/2)
luma_dist_matrix = (dist2(as.matrix(subset_luma),as.matrix(subset_luma)))^(1/2)
lumb_dist_matrix = (dist2(as.matrix(subset_lumb),as.matrix(subset_lumb)))^(1/2)
normal_dist_matrix = (dist2(as.matrix(subset_normal),as.matrix(subset_normal)))^(1/2)
her2_dist_matrix = (dist2(as.matrix(subset_her2),as.matrix(subset_her2)))^(1/2)

#Constructing the similarity graphs 
W1 = affinityMatrix(basal_dist_matrix, K, alpha)
W2 = affinityMatrix(luma_dist_matrix, K, alpha)
W3 = affinityMatrix(lumb_dist_matrix, K, alpha)
W4 = affinityMatrix(normal_dist_matrix, K, alpha)
W5 = affinityMatrix(her2_dist_matrix, K, alpha)

## then the overall matrix can be computed by similarity network fusion(SNF):
W = SNF(list(W1,W2, W3, W4, W5), K, T)

#Truelabel - is the PAM50 subtypes 
truelabel = c(matrix(1,100,1),matrix(2,100,1));

## You can display clusters in the data by the following function
## where C is the number of clusters.
C = 2 # number of clusters
group = spectralClustering(W_Basal,C); # the final subtypes information

## Get a matrix containing the group information
## for the samples such as the SpectralClustering result and the True label
M_label=cbind(group,truelabel_krug)
colnames(M_label)=c("spectralClustering","TrueLabel")


## Use the getColorsForGroups function to assign a color to each group
## NB is more than 8 groups, you will have to input a vector
## of colors into the getColorsForGroups function
M_label_colors=t(apply(M_label,1,getColorsForGroups))
## or choose you own colors for each label, for example:
M_label_colors=cbind("spectralClustering"=getColorsForGroups(M_label[,"spectralClustering"],
                                                             colors=c("blue","green")),"TrueLabel"=getColorsForGroups(M_label[,"TrueLabel"],

                                                                    
                                                                                                                      