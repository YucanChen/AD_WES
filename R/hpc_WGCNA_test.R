library(openxlsx)
library(WGCNA)
options(stringsAsFactors = FALSE)
enableWGCNAThreads()

load("./ALLENexpdata_exploration.RData")

burdenresult_dementia<-read.csv("./metadata/burdenresult_dementia_annotatedforLR.csv")
burdenresult_AD<-read.csv("./metadata/burdenresult_AD_annotatedforLR.csv")
burdenresult_earlyAD<-read.csv("./metadata/burdenresult_earlyAD_annotatedforLR.csv")

datexp=as.matrix(ALLENexp.ADNorm)
datexp=datexp[which(rowSums(datexp) > 0),]

### WGCNA NetworkConstruction using AD vs Normal data ###

#===============================================================================
#
#  Code chunk 1: input data preparation
#
#===============================================================================

# prepare regional expression data

## The expression data from the FWM region is greatly different from those of the other three regions in network topology, and thus was not included in the later analyses.

setLabels = c("TCx","HIP","PCx")

TCx_colsamples<-ALLENexp.ADNorm.col[which(ALLENexp.ADNorm.col$structure_acronym==setLabels[1]),]
TCx_colsamples$rnaseq_profile_id=paste("X",TCx_colsamples$rnaseq_profile_id,sep = '')
datexp_TCx=datexp[,which(colnames(datexp) %in% TCx_colsamples$rnaseq_profile_id)]

HIP_colsamples<-ALLENexp.ADNorm.col[which(ALLENexp.ADNorm.col$structure_acronym==setLabels[2]),]
HIP_colsamples$rnaseq_profile_id=paste("X",HIP_colsamples$rnaseq_profile_id,sep = '')
datexp_HIP=datexp[,which(colnames(datexp) %in% HIP_colsamples$rnaseq_profile_id)]

PCx_colsamples<-ALLENexp.ADNorm.col[which(ALLENexp.ADNorm.col$structure_acronym==setLabels[3]),]
PCx_colsamples$rnaseq_profile_id=paste("X",PCx_colsamples$rnaseq_profile_id,sep = '')
datexp_PCx=datexp[,which(colnames(datexp) %in% PCx_colsamples$rnaseq_profile_id)]

datexp_TCx=datexp_TCx[which(rowSums(datexp_TCx)>1),]
datexp_HIP=datexp_HIP[which(rowSums(datexp_HIP)>1),]
datexp_PCx=datexp_PCx[which(rowSums(datexp_PCx)>1),]

## filter genes expressed in all three brain regions

datexp_allregion_genes<-c(rownames(datexp_TCx),rownames(datexp_HIP),rownames(datexp_PCx))
datexp_allregion_genes.table<-as.data.frame(table(datexp_allregion_genes))
summary(datexp_allregion_genes.table$Freq)
datexp_allregion_genes.table<-datexp_allregion_genes.table[which(datexp_allregion_genes.table$Freq==3),]
summary(datexp_allregion_genes.table$Freq)
datexp_allregion_genes<-datexp_allregion_genes.table$datexp_allregion_genes

datexp_TCx=datexp_TCx[which(rownames(datexp_TCx) %in% datexp_allregion_genes),]
datexp_HIP=datexp_HIP[which(rownames(datexp_HIP) %in% datexp_allregion_genes),]
datexp_PCx=datexp_PCx[which(rownames(datexp_PCx) %in% datexp_allregion_genes),]

datexp_TCx<-t(datexp_TCx)
datexp_HIP<-t(datexp_HIP)
datexp_PCx<-t(datexp_PCx)

datexp<-datexp[which(rownames(datexp) %in% datexp_allregion_genes),]
datexp<-t(datexp)

Disease_state=unlist(lapply(condition,function(x){if(x=="NoDementia"){0}else{1}}))
brain_region<-ALLENexp.ADNorm.col$structure_acronym
meta=cbind(id=rownames(datexp),Disease_state=Disease_state,Brain_region=brain_region)

nSets=3

## Form multi-set expression data

multiExpr = vector(mode = "list", length = nSets)

multiExpr[[1]] = list(data = as.data.frame(datexp_TCx))
names(multiExpr[[1]]$data) = colnames(datexp_TCx)
rownames(multiExpr[[1]]$data) = rownames(datexp_TCx)
multiExpr[[2]] = list(data = as.data.frame(datexp_HIP))
names(multiExpr[[2]]$data) = colnames(datexp_HIP)
rownames(multiExpr[[2]]$data) = rownames(datexp_HIP)
multiExpr[[3]] = list(data = as.data.frame(datexp_PCx))
names(multiExpr[[3]]$data) = colnames(datexp_PCx)
rownames(multiExpr[[3]]$data) = rownames(datexp_PCx)

# check the data and remove outlier

## check that the data has the correct format for many functions operating on multiple sets:

exprSize = checkSets(multiExpr)
exprSize

## cluster the samples on their Euclidean distance, separately in each set.

sampleTrees = list()
for (set in 1:nSets)
{
  sampleTrees[[set]] = hclust(dist(multiExpr[[set]]$data), method = "average")
}

pdf(file = "./results/Plots/ALLEN_ADvsNorm_SampleClustering.pdf", width = 12, height = 12);
par(mfrow=c(1,1))
par(mar = c(0, 4, 2, 0))
for (set in 1:nSets)
  plot(sampleTrees[[set]], main = paste("Sample clustering on all genes in", setLabels[set]),xlab="", sub="", cex = 0.7,cex.main=1, font.main=1 )
dev.off()

## Choose the "base" cut height for the female data set

baseHeight = 83
# Adjust the cut height for the male data set for the number of samples
cutHeights = c(83,83,83)
# Re-plot the dendrograms including the cut lines
pdf(file = "./results/Plots/ALLEN_ADvsNorm_SampleClustering_cutheightadded.pdf", width = 12, height = 12);
par(mfrow=c(2,1))
par(mar = c(0, 4, 2, 0))
for (set in 1:nSets)
{
  plot(sampleTrees[[set]], main = paste("Sample clustering on all genes in", setLabels[set]),
       xlab="", sub="", cex = 0.7);
  abline(h=cutHeights[set], col = "red");
}
dev.off()

## actual outlier removal

for (set in 1:nSets)
{
  # Find clusters cut by the line
  labels = cutreeStatic(sampleTrees[[set]], cutHeight = cutHeights[set])
  # Keep the largest one (labeled by the number 1)
  keep = (labels==1)
  multiExpr[[set]]$data = multiExpr[[set]]$data[keep, ]
}
collectGarbage()

## Check the size of the leftover data

exprSize = checkSets(multiExpr)
exprSize

## again cluster the samples on their Euclidean distance, separately in each set.

sampleTrees = list()
for (set in 1:nSets)
{
  sampleTrees[[set]] = hclust(dist(multiExpr[[set]]$data), method = "average")
}

pdf(file = "./results/Plots/ALLEN_ADvsNorm_SampleClustering_outlierdeleted.pdf", width = 12, height = 12);
par(mfrow=c(1,1))
par(mar = c(0, 4, 2, 0))
for (set in 1:nSets)
  plot(sampleTrees[[set]], main = paste("Sample clustering on all genes in", setLabels[set]),xlab="", sub="", cex = 0.7,cex.main=1, font.main=1 )
dev.off()

#===============================================================================
#
#  Code chunk 2: step-by-step network construction
#
#===============================================================================

# Choosing the soft-thresholding power: analysis of network topology

# Choose a set of soft-thresholding powers
powers = c(seq(4,10,by=1), seq(12,20, by=2))
# Initialize a list to hold the results of scale-free analysis
powerTables = vector(mode = "list", length = nSets)
# Call the network topology analysis function for each set in turn
for (set in 1:nSets)
  powerTables[[set]] = list(data = pickSoftThreshold(multiExpr[[set]]$data, powerVector=powers,
                                                     verbose = 2)[[2]]);
collectGarbage()
# Plot the results:
colors = c("black", "red","blue","green")
# Will plot these columns of the returned scale free analysis tables
plotCols = c(2,5,6,7)
colNames = c("Scale Free Topology Model Fit", "Mean connectivity", "Median connectivity",
             "Max connectivity");
# Get the minima and maxima of the plotted points
ylim = matrix(NA, nrow = 2, ncol = 4);
for (set in 1:nSets)
{
  for (col in 1:length(plotCols))
  {
    ylim[1, col] = min(ylim[1, col], powerTables[[set]]$data[, plotCols[col]], na.rm = TRUE);
    ylim[2, col] = max(ylim[2, col], powerTables[[set]]$data[, plotCols[col]], na.rm = TRUE);
  }
}

# Plot the quantities in the chosen columns vs. the soft thresholding power
sizeGrWindow(8, 6)
pdf(file = "./results/Plots/ALLEN_ADvsNorm_powerselection_connectivity_topology_model.pdf", wi = 12, he = 12);
par(mfcol = c(2,2));
par(mar = c(4.2, 4.2 , 2.2, 0.5))
cex1 = 0.7;
for (col in 1:length(plotCols)) for (set in 1:nSets)
{
  if (set==1)
  {
    plot(powerTables[[set]]$data[,1], -sign(powerTables[[set]]$data[,3])*powerTables[[set]]$data[,2],
         xlab="Soft Threshold (power)",ylab=colNames[col],type="n", ylim = ylim[, col],
         main = colNames[col]);
    addGrid();
  }
  if (col==1)
  {
    text(powerTables[[set]]$data[,1], -sign(powerTables[[set]]$data[,3])*powerTables[[set]]$data[,2],
         labels=powers,cex=cex1,col=colors[set]);
  } else
    text(powerTables[[set]]$data[,1], powerTables[[set]]$data[,plotCols[col]],
         labels=powers,cex=cex1,col=colors[set]);
  if (col==1)
  {
    legend("bottomright", legend = setLabels, col = colors, pch = 20) ;
  } else
    legend("topright", legend = setLabels, col = colors, pch = 20) ;
}
dev.off()

# Calculation of network adjacencies

## Network construction starts by calculating the adjacencies in the individual sets

softPower = 6;
nGenes = exprSize$nGenes;
nSamples = exprSize$nSamples;

## Initialize an appropriate array to hold the adjacencies

adjacencies = array(0, dim = c(nSets, nGenes, nGenes))
# Calculate adjacencies in each individual data set
for (set in 1:nSets)
  adjacencies[set, , ] = abs(cor(multiExpr[[set]]$data, use = "p"))^softPower;

# Calculation of Topological Overlap

## Initialize an appropriate array to hold the TOMs

TOM = array(0, dim = c(nSets, nGenes, nGenes));
# Calculate TOMs in each individual data set
for (set in 1:nSets)
  TOM[set, , ] = TOMsimilarity(adjacencies[set, , ]);

# Scaling of Topological Overlap Matrices to make them comparable across sets

## Define the reference percentile

scaleP = 0.95
# Set RNG seed for reproducibility of sampling
set.seed(12345)
# Sample sufficiently large number of TOM entries
nSamples = as.integer(1/(1-scaleP) * 1000);
# Choose the sampled TOM entries
scaleSample = sample(nGenes*(nGenes-1)/2, size = nSamples)
TOMScalingSamples = list();
# These are TOM values at reference percentile
scaleQuant = rep(1, nSets)

## Scaling powers to equalize reference TOM values

scalePowers = rep(1, nSets)
# Loop over sets
for (set in 1:nSets)
{
### Select the sampled TOM entries
TOMScalingSamples[[set]] = as.dist(TOM[set, , ])[scaleSample]
### Calculate the 95th percentile
scaleQuant[set] = quantile(TOMScalingSamples[[set]],
probs = scaleP, type = 8);
### Scale the male TOM
if (set>1)
{
scalePowers[set] = log(scaleQuant[1])/log(scaleQuant[set]);
TOM[set, ,] = TOM[set, ,]^scalePowers[set];
}
}

## For plotting, also scale the sampled TOM entries

scaledTOMSamples = list();
for (set in 1:nSets)
  scaledTOMSamples[[set]] = TOMScalingSamples[[set]]^scalePowers[set]

sizeGrWindow(6,6)
pdf(file = "./results/Plots/TOMScaling-QQPlot_TCx_HIP.pdf", wi = 6, he = 6);
# qq plot of the unscaled samples
qqUnscaled = qqplot(TOMScalingSamples[[1]], TOMScalingSamples[[2]], plot.it = TRUE, cex = 0.6, xlab = paste("TOM in", setLabels[1]), ylab = paste("TOM in", setLabels[2]), main = "Q-Q plot of TOM", pch = 20)
# qq plot of the scaled samples
qqScaled = qqplot(scaledTOMSamples[[1]], scaledTOMSamples[[2]], plot.it = FALSE)
points(qqScaled$x, qqScaled$y, col = "red", cex = 0.6, pch = 20);
abline(a=0, b=1, col = "blue")
legend("topleft", legend = c("Unscaled TOM", "Scaled TOM"), pch = 20, col = c("black", "red"))
dev.off()

sizeGrWindow(6,6)
pdf(file = "./results/Plots/TOMScaling-QQPlot_HIP_PCx.pdf", wi = 6, he = 6);
# qq plot of the unscaled samples
qqUnscaled = qqplot(TOMScalingSamples[[2]], TOMScalingSamples[[3]], plot.it = TRUE, cex = 0.6, xlab = paste("TOM in", setLabels[2]), ylab = paste("TOM in", setLabels[3]), main = "Q-Q plot of TOM", pch = 20)
# qq plot of the scaled samples
qqScaled = qqplot(scaledTOMSamples[[2]], scaledTOMSamples[[3]], plot.it = FALSE)
points(qqScaled$x, qqScaled$y, col = "red", cex = 0.6, pch = 20);
abline(a=0, b=1, col = "blue")
legend("topleft", legend = c("Unscaled TOM", "Scaled TOM"), pch = 20, col = c("black", "red"))
dev.off()

sizeGrWindow(6,6)
pdf(file = "./results/Plots/TOMScaling-QQPlot_TCx_PCx.pdf", wi = 6, he = 6);
# qq plot of the unscaled samples
qqUnscaled = qqplot(TOMScalingSamples[[1]], TOMScalingSamples[[3]], plot.it = TRUE, cex = 0.6, xlab = paste("TOM in", setLabels[1]), ylab = paste("TOM in", setLabels[3]), main = "Q-Q plot of TOM", pch = 20)
# qq plot of the scaled samples
qqScaled = qqplot(scaledTOMSamples[[1]], scaledTOMSamples[[3]], plot.it = FALSE)
points(qqScaled$x, qqScaled$y, col = "red", cex = 0.6, pch = 20);
abline(a=0, b=1, col = "blue")
legend("topleft", legend = c("Unscaled TOM", "Scaled TOM"), pch = 20, col = c("black", "red"))
dev.off()

# Calculation of consensus Topological Overlap

consensusTOM = pmin(TOM[1, , ],TOM[2, , ],TOM[3, , ]);

# Clustering and module identification

# Clustering
consTree = hclust(as.dist(1-consensusTOM), method = "average");
# We like large modules, so we set the minimum module size relatively high:
minModuleSize = 30;
# Module identification using dynamic tree cut:
unmergedLabels = cutreeDynamic(dendro = consTree, distM = 1-consensusTOM,
deepSplit = 2, cutHeight = 0.995,
minClusterSize = minModuleSize,
pamRespectsDendro = FALSE );
unmergedColors = labels2colors(unmergedLabels)

sizeGrWindow(8,6)
pdf(file = "./results/Plots/Dynamic_Tree_Cut.pdf", wi = 6, he = 6);
plotDendroAndColors(consTree, unmergedColors, "Dynamic Tree Cut",
dendroLabels = FALSE, hang = 0.03,
addGuide = TRUE, guideHang = 0.05)
dev.off()

# Merging of modules whose expression profiles are very similar

# Calculate module eigengenes
unmergedMEs = multiSetMEs(multiExpr, colors = NULL, universalColors = unmergedColors, excludeGrey = TRUE)
# Calculate consensus dissimilarity of consensus module eigengenes
consMEDiss = consensusMEDissimilarity(unmergedMEs);
# Cluster consensus modules
consMETree = hclust(as.dist(consMEDiss), method = "average");
# Plot the result
sizeGrWindow(7,6)
pdf(file = "./results/Plots/Consensus_clustering_of_consensus_module_eigengenes.pdf", wi = 6, he = 6);
par(mfrow = c(1,1))
plot(consMETree, main = "Consensus clustering of consensus module eigengenes",
xlab = "", sub = "")
abline(h=0.25, col = "red")
dev.off()

merge = mergeCloseModules(multiExpr, unmergedLabels, cutHeight = 0.25, verbose = 3)

# Numeric module labels
moduleLabels = merge$colors;
# Convert labels to colors
moduleColors = labels2colors(moduleLabels)
# Eigengenes of the new merged modules:
consMEs = merge$newMEs;
# remove "grey" module
consMEs = removeGreyME(consMEs, greyMEName = paste(moduleColor.getMEprefix(), "0", sep=""))

# plot the gene dendrogram again, this time with both the unmerged and the merged module colors
sizeGrWindow(9,6)
pdf(file = "./results/Plots/Consensus_clustering_of_consensus_module_eigengenes_new.pdf", wi = 6, he = 6);
plotDendroAndColors(consTree, cbind(unmergedColors, moduleColors),
c("Unmerged", "Merged"),
dendroLabels = FALSE, hang = 0.03,
addGuide = TRUE, guideHang = 0.05)
dev.off()

save(consMEs, moduleColors, moduleLabels, consTree, file = "./Consensus-NetworkConstruction-man.RData")

#===============================================================================
#
#  Code chunk 3: relate consensus modules to external trait info
#
#===============================================================================

Traits = vector(mode = "list", length = 3) #nSets

Traits[[1]] = list(data = as.data.frame(TCx_colsamples))
names(Traits[[1]]$data) = colnames(TCx_colsamples)
rownames(Traits[[1]]$data) = rownames(TCx_colsamples)
Traits[[2]] = list(data = as.data.frame(HIP_colsamples))
names(Traits[[2]]$data) = colnames(HIP_colsamples)
rownames(Traits[[2]]$data) = rownames(HIP_colsamples)
Traits[[3]] = list(data = as.data.frame(PCx_colsamples))
names(Traits[[3]]$data) = colnames(PCx_colsamples)
rownames(Traits[[3]]$data) = rownames(PCx_colsamples)

# Set up variables to contain the module-trait correlations

moduleTraitCor = list();
moduleTraitPvalue = list();
# Calculate the correlations
for (set in 1:nSets)
{
moduleTraitCor[[set]] = cor(consMEs[[set]]$data, Traits[[set]]$data, use = "p");
moduleTraitPvalue[[set]] = corPvalueFisher(moduleTraitCor[[set]], exprSize$nSamples[set]);
}
# Convert numerical lables to colors for labeling of modules in the plot
MEColors = labels2colors(as.numeric(substring(names(consMEs[[1]]$data), 3)));
MEColorNames = paste("ME", MEColors, sep="");

## Plot the module-trait relationship table for set number 1
set = 1
textMatrix = paste(signif(moduleTraitCor[[set]], 2), "\n(",
signif(moduleTraitPvalue[[set]], 1), ")", sep = "");
dim(textMatrix) = dim(moduleTraitCor[[set]])
sizeGrWindow(10,7)
pdf(file = "./results/Plots/ModuleTraitRelationships-TCx.pdf", wi = 6, he = 12);
par(mar = c(6, 8.8, 3, 2.2));
labeledHeatmap(Matrix = moduleTraitCor[[set]],
xLabels = names(Traits[[set]]$data),
yLabels = MEColorNames,
ySymbols = MEColorNames,
colorLabels = FALSE,
colors = greenWhiteRed(50),
textMatrix = textMatrix,
setStdMargins = FALSE,
cex.text = 0.5,
zlim = c(-1,1),
main = paste("Module--trait relationships in", setLabels[set]))
dev.off();

## Plot the module-trait relationship table for set number 2
set = 2
textMatrix = paste(signif(moduleTraitCor[[set]], 2), "\n(",
signif(moduleTraitPvalue[[set]], 1), ")", sep = "");
dim(textMatrix) = dim(moduleTraitCor[[set]])
sizeGrWindow(10,7)
pdf(file = "./results/Plots/ModuleTraitRelationships-HIP.pdf", wi = 6, he = 12);
par(mar = c(6, 8.8, 3, 2.2));
labeledHeatmap(Matrix = moduleTraitCor[[set]],
xLabels = names(Traits[[set]]$data),
yLabels = MEColorNames,
ySymbols = MEColorNames,
colorLabels = FALSE,
colors = greenWhiteRed(50),
textMatrix = textMatrix,
setStdMargins = FALSE,
cex.text = 0.5,
zlim = c(-1,1),
main = paste("Module--trait relationships in", setLabels[set]))
dev.off();

## Plot the module-trait relationship table for set number 3
set = 3
textMatrix = paste(signif(moduleTraitCor[[set]], 2), "\n(",
signif(moduleTraitPvalue[[set]], 1), ")", sep = "");
dim(textMatrix) = dim(moduleTraitCor[[set]])
sizeGrWindow(10,7)
pdf(file = "./results/Plots/ModuleTraitRelationships-PCx.pdf", wi = 6, he = 12);
par(mar = c(6, 8.8, 3, 2.2));
labeledHeatmap(Matrix = moduleTraitCor[[set]],
xLabels = names(Traits[[set]]$data),
yLabels = MEColorNames,
ySymbols = MEColorNames,
colorLabels = FALSE,
colors = greenWhiteRed(50),
textMatrix = textMatrix,
setStdMargins = FALSE,
cex.text = 0.5,
zlim = c(-1,1),
main = paste("Module--trait relationships in", setLabels[set]))
dev.off();

# Initialize matrices to hold the consensus correlation and p-value

consensusCor = matrix(NA, nrow(moduleTraitCor[[1]]), ncol(moduleTraitCor[[1]]));
consensusPvalue = matrix(NA, nrow(moduleTraitCor[[1]]), ncol(moduleTraitCor[[1]]));
# Find consensus negative correlations
negative = moduleTraitCor[[1]] < 0 & moduleTraitCor[[2]] < 0 & moduleTraitCor[[3]] < 0;
consensusCor[negative] = pmax(moduleTraitCor[[1]][negative], moduleTraitCor[[2]][negative], moduleTraitCor[[3]][negative]);
consensusPvalue[negative] = pmax(moduleTraitPvalue[[1]][negative], moduleTraitPvalue[[2]][negative], moduleTraitPvalue[[3]][negative]);
# Find consensus positive correlations
positive = moduleTraitCor[[1]] > 0 & moduleTraitCor[[2]] > 0 & moduleTraitCor[[3]] > 0;
consensusCor[positive] = pmin(moduleTraitCor[[1]][positive], moduleTraitCor[[2]][positive], moduleTraitCor[[3]][positive]);
consensusPvalue[positive] = pmax(moduleTraitPvalue[[1]][positive], moduleTraitPvalue[[2]][positive], moduleTraitPvalue[[3]][positive]);

textMatrix = paste(signif(consensusCor, 2), "\n(",
signif(consensusPvalue, 1), ")", sep = "");
dim(textMatrix) = dim(moduleTraitCor[[set]])
sizeGrWindow(10,7)
pdf(file = "./results/Plots/ModuleTraitRelationships-consensus.pdf", wi = 6, he = 12);
par(mar = c(6, 8.8, 3, 2.2));
labeledHeatmap(Matrix = consensusCor,
xLabels = names(Traits[[set]]$data),
yLabels = MEColorNames,
ySymbols = MEColorNames,
colorLabels = FALSE,
colors = greenWhiteRed(50),
textMatrix = textMatrix,
setStdMargins = FALSE,
cex.text = 0.5,
zlim = c(-1,1),
main = paste("Consensus module--trait relationships across\n",
paste(setLabels, collapse = " and ")))
dev.off();

save(moduleTraitCor, moduleTraitPvalue, consensusCor, consensusPvalue, file = "./Consensus-Network_module-trait-association.RData")

#===============================================================================
#
#  Code chunk 4: Exporting results of the network analysis
#
#===============================================================================

consMEs.unord = multiSetMEs(multiExpr, universalColors = moduleLabels, excludeGrey = TRUE)
GS = list();
kME = list();
for (set in 1:nSets)
{
GS[[set]] = corAndPvalue(multiExpr[[set]]$data, Traits[[set]]$data);
kME[[set]] = corAndPvalue(multiExpr[[set]]$data, consMEs.unord[[set]]$data);
}

GS.metaZ = (GS[[1]]$Z + GS[[2]]$Z + GS[[3]]$Z)/sqrt(3);
kME.metaZ = (kME[[1]]$Z + kME[[2]]$Z + kME[[3]]$Z)/sqrt(3);
GS.metaP = 3*pnorm(abs(GS.metaZ), lower.tail = FALSE);
kME.metaP = 3*pnorm(abs(kME.metaZ), lower.tail = FALSE);

GSmat = rbind(GS[[1]]$cor, GS[[2]]$cor, GS[[3]]$cor, GS[[1]]$p, GS[[2]]$p, GS[[3]]$p, GS.metaZ, GS.metaP);
nTraits = checkSets(Traits)$nGenes
traitNames = colnames(Traits[[1]]$data)
dim(GSmat) = c(nGenes, 8*nTraits)
rownames(GSmat) = names(multiExpr[[1]]$data);
colnames(GSmat) = spaste(
c("GS.set1.", "GS.set2.", "GS.set3.", "p.GS.set1.", "p.GS.set2.", "p.GS.set3.", "Z.GS.meta.", "p.GS.meta"),
rep(traitNames, rep(8, nTraits)))
# Same code for kME:
kMEmat = rbind(kME[[1]]$cor, kME[[2]]$cor, kME[[3]]$cor, kME[[1]]$p, kME[[2]]$p, kME[[3]]$p, kME.metaZ, kME.metaP);
MEnames = colnames(consMEs.unord[[1]]$data);
nMEs = checkSets(consMEs.unord)$nGenes
dim(kMEmat) = c(nGenes, 8*nMEs)
rownames(kMEmat) = names(multiExpr[[1]]$data);
colnames(kMEmat) = spaste(
c("kME.set1.", "kME.set2.", "kME.set3.", "p.kME.set1.", "p.kME.set2.", "p.kME.set3.", "Z.kME.meta.", "p.kME.meta"),
rep(MEnames, rep(8, nMEs)))

info = data.frame(EntrezID = names(multiExpr[[1]]$data),
ModuleLabel = moduleLabels,
ModuleColor = labels2colors(moduleLabels),
GSmat,
kMEmat);
write.csv(info, file = "./results/consensusAnalysis-CombinedNetworkResults.csv",
row.names = FALSE, quote = FALSE);

###
save.image("./WGCNA_brainregioncombined_enrichment_analysis.RData")
