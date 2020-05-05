library(Seurat)

args = commandArgs(trailingOnly=T)

f_input = args[1]
f_output = args[2]
p_prefix = args[3]

v2 <- readRDS(f_input)
v3 <- UpdateSeuratObject(v2)

adt <- as.matrix(GetAssayData(v3, assay="CITE", slot="counts"))
rownames(adt) <- gsub(" ", "-", rownames(adt))

v3[["CITE"]] <- CreateAssayObject(counts = adt)
v3 <- NormalizeData(v3, assay="CITE", normalization.method = "CLR", margin = 2)
v3 <- ScaleData(v3, assay="CITE")


# =================
# Filter Artifacts
# =================
# 1. Scrublet
scrublet <- v3@meta.data$scrublet

# 2. S2
s2 <- v3@meta.data$S2

# 3. Remove cells with no ADT counts
ADT <- GetAssayData(v3, assay="CITE")
cs <- colSums(ADT)

# 4. Subset by cellnames
v3 <- subset(v3, cells=names(cs[(cs!=0)&(scrublet==0)&(s2==0)]))

# =================
# JointVis (Joint)
# =================
source("/home/FCAM/ycheng/projects/koralov/src/JointVis/utils_updated_v3.R")
require(methods); library(Rcpp); library(RcppArmadillo)
sourceCpp("/home/FCAM/ycheng/projects/koralov/src/JointVis/utils.cpp")

snn_type = "SNF"
data_type = c("GEX", "CITE")
cluster_method = "Spectral"
NF_method = "ANF"
iter = 2
p_Linf = 1
p_clusters = max(as.integer(v3@meta.data$res.0.8))+1
JointRes <- JointVisSeuratV3(v3, p_clusters,
                             snn_type, cluster_method, NF_method,
                             data = data_type, iter = iter, L = p_Linf)
output <- prep_result(JointRes)
v3 <- AddMetaData(v3, metadata=factor(output[["cluster"]], sort(as.numeric(levels(output[["cluster"]])))), col.name="SpectralJoint")
v3[["jointtsne"]] <- CreateDimReducObject(embeddings=output[["embeddings"]],
	                                  key="jointtsne_",
				          assay=DefaultAssay(v3))
Idents(v3) <- factor(output[["cluster"]],
	             sort(as.numeric(levels(output[["cluster"]]))))

# ===============
# JointVis (ADT)
# ===============
data_type = c("CITE")
ClusterADT <- JointVisSeuratV3(v3, p_clusters,
                               snn_type, cluster_method, NF_method,
                               data = data_type, iter = iter, L = p_Linf)
output_ADT <- prep_result(ClusterADT)

v3 <- AddMetaData(v3, metadata=factor(output_ADT[["cluster"]], sort(as.numeric(levels(output_ADT[["cluster"]])))), col.name="SpectralADT")
v3[["adttsne"]] <- CreateDimReducObject(embeddings=output_ADT[["embeddings"]],
	                                 key="adttsne_",
				         assay=DefaultAssay(v3))

# =================
# JointVis (GEX)
# =================
data_type = c("GEX")
ClusterGEX <- JointVisSeuratV3(v3, p_clusters,
                               snn_type, cluster_method, NF_method,
                               data = data_type, iter = iter, L = p_Linf)
output_GEX <- prep_result(ClusterGEX)
v3 <- AddMetaData(v3, metadata=factor(output_GEX[["cluster"]], sort(as.numeric(levels(output_GEX[["cluster"]])))), col.name="SpectralGEX")
v3[["gextsne"]] <- CreateDimReducObject(embeddings=output_GEX[["embeddings"]],
	                                key="gextsne_",
				        assay=DefaultAssay(v3))


# =================
# Run uMAP
# =================
library(umap)
run_umap <- function(output) {
  distmat <- 1-output[["SN"]]
  custom.config <- umap.defaults
  custom.config$input <- "dist"
  custom.config$random_state <- 386
  u1 <- umap(distmat, custom.config)
  rownames(u1$layout) <- rownames(output[["embeddings"]])
  colnames(u1$layout) <- c("umap_1", "umap_2")
  return(u1)
}

umap_joint <- run_umap(output)
umap_adt <- run_umap(output_ADT)
umap_gex <- run_umap(output_GEX)

v3[["jointumap"]] <- CreateDimReducObject(embeddings=umap_joint$layout,
		                          key="jointumap_",
					  assay=DefaultAssay(v3))
v3[["adtumap"]] <- CreateDimReducObject(embeddings=umap_adt$layout,
		                        key="adtumap_",
					assay=DefaultAssay(v3))
v3[["gexumap"]] <- CreateDimReducObject(embeddings=umap_gex$layout,
		                        key="gexumap_",
					assay=DefaultAssay(v3))

# ==============
# Save output
# ==============
saveRDS(output, paste(p_prefix, "_Joint_ANF_SNF_Spectral_", p_clusters, "_T_output.rds", sep=""))
saveRDS(output_ADT, paste(p_prefix, "_ADT_ANF_SNF_Spectral_", p_clusters, "_T_output.rds", sep=""))
saveRDS(output_GEX, paste(p_prefix, "_GEX_ANF_SNF_Spectral_", p_clusters, "_T_output.rds", sep=""))
saveRDS(v3, f_output)
