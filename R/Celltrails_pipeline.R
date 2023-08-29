library(CellTrails)
library(ggplotify)
library(tidyverse)

# Import in vitro data and subset HCs & SCs
invitro_haircell <- readRDS("~/aggr/Invitro_22clusters.rds")
invitro_HC_SC = subset(invitro_haircell,idents = c(3,4,5,14,15))
data <- as.matrix(invitro_HC_SC@assays$RNA@data)
data_sub = data[apply(data,1,FUN = function(x) sum(x) != 0),]
HC_SC <- SingleCellExperiment(assays=list(logcounts=data_sub))

# Select trajectory features
trajFeatureNames(HC_SC) <- featureNames(HC_SC)
trajFeatureNames(HC_SC) <- filterTrajFeaturesByFF(HC_SC, threshold=0.7)
tfeat <- featureNames(HC_SC[trajFeatureNames(HC_SC), ])
trajFeatureNames(HC_SC) <- tfeat
showTrajInfo(HC_SC)

# Perform spectral embedding
se <- embedSamples(HC_SC)
d <- findSpectrum(se$eigenvalues, frac=100)
latentSpace(HC_SC) <- se$components[, d]

# Perform clustering
cl <- findStates(HC_SC, min_size=0.02, min_feat=5, max_pval=1e-4, min_fc=1.8)
states(HC_SC) <- cl

# Sample ordering
HC_SC <- connectStates(HC_SC, l=10)
HC_SC <- fitTrajectory(HC_SC)

# Export trajectory graph structure to graphml and then import layout computed from yEd
write.ygraphml(HC_SC, file="~/celltrails/HC_SC_8states.graphml", color_by="phenoName",name="state", node_label="state")
tl <- read.ygraphml(file="~/celltrails/HC_SC_8states_layout.graphml")

# Adjust layout and store to object
trajLayout(HC_SC, adjust=TRUE) <- tl

plotMap(HC_SC, color_by="phenoName", name="state")

saveRDS(HC_SC, "/home/rstudio/keith_data/Invitro_HC_SC_sce.rds")

