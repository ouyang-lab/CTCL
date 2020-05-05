path = "/home/FCAM/ycheng/projects/koralov"
setwd(path)

library(Seurat)
library(magrittr)
library(dplyr)
library(tibble)

args = commandArgs(trailingOnly=T)

sample = args[1]
sample_short = unlist(strsplit(sample, "_"))[1]
f_output = args[2]

# applicable from H19 onwards

if (sample_short %in% c("H5", "H6", "H7", "H8", "H9")) {
  sid = sample
} else {
  sid = unlist(strsplit(sample, "_"))[2]
}

if (sample_short %in% c("H5", "H6", "H7", "H8", "H9", "H11")) {
  cDNA = paste("data/cellranger/cDNA", paste(sample, "cDNA", sep="_"), 
               "outs", "filtered_gene_bc_matrices_h5_genes.csv.gz", sep="/")
} else {
  cDNA = paste("data/cellranger3/cDNA-HTO", paste(sample, "cDNA-HTO", sep="_"), 
               "outs", "filtered_feature_bc_matrix_genes.csv.gz", sep="/")
}

ADT = paste("data/ADT", paste(sample, "ADT.csv.gz", sep="_"), sep="/")
Tab = paste("data/VDJ", paste(sample, "TCRab.csv.gz", sep="_"), sep="/")
Tgd = paste("data/VDJ", paste(sample, "TCRgd.csv.gz", sep="_"), sep="/")
HTO = paste("data/HTO", paste(sample, "HTO.csv.gz", sep="_"), sep="/")

Cab = paste("data/VDJ", paste(sample, "clonotypeAB.csv.gz", sep="_"), sep="/")
Cgd = paste("data/VDJ", paste(sample, "clonotypeGD.csv.gz", sep="_"), sep="/")

print(cDNA)
print(ADT)
print(Tab)
print(Tgd)
print(HTO)
print(Cab)
print(Cgd)

# =========================================================================================
# Assign Sample Identifier and Remove HTO duplicates
# =========================================================================================
csv.cdna <- read.csv(cDNA, sep = ",", header = TRUE, row.names = 1)
colnames(csv.cdna) <- paste(sample_short, colnames(csv.cdna), sep = '_')

tmp <- (CreateSeuratObject(raw.data = csv.cdna, project=sid)
         %>% NormalizeData()
	 %>% ScaleData(display.progress = TRUE))
tmp@meta.data$stim <- sid

print("load cDNA")

# =========================================================================================
# Import HTO into Seurat Object
# =========================================================================================
csv.hto <- read.csv(HTO, sep = ",", header = TRUE, row.names = 1)
colnames(csv.hto) <- paste(sample_short, colnames(csv.hto), sep = '_')
rownames(csv.hto) <- paste0("HTO_", rownames(csv.hto))

tmp <- (tmp
          %>% SetAssayData(assay.type = "HTO", slot = "raw.data", new.data = csv.hto)
	  %>% NormalizeData(assay.type = "HTO", normalization.method = "genesCLR")
          %>% ScaleData(assay.type = "HTO", display.progress = TRUE)	  
	  %>% HTODemux(positive_quantile = 0.999))

print("load HTO")

# =========================================================================================
# Import TCRab into Seurat Object
# =========================================================================================
csv.TCRab <- read.csv(Tab, sep = ",", header = TRUE, row.names = 1)
colnames(csv.TCRab) <- paste(sample_short, colnames(csv.TCRab), sep = '_')
rownames(csv.TCRab) <- paste0("TCRab_", rownames(csv.TCRab))

tmp <- (tmp
          %>% SetAssayData(assay.type = "TCRab", slot = "raw.data", new.data = csv.TCRab)
	  %>% NormalizeData(assay.type = "TCRab", normalization.method = "genesCLR")
	  %>% ScaleData(assay.type = "TCRab", display.progress = TRUE))

print("load TCRab")

# =========================================================================================
# Import TCRgd into Seurat Object
# =========================================================================================
csv.TCRgd <- read.csv(Tgd, sep = ",", header = TRUE, row.names = 1)
colnames(csv.TCRgd) <- paste(sample_short, colnames(csv.TCRgd), sep = '_')
rownames(csv.TCRgd) <- paste0("TCRgd_", rownames(csv.TCRgd))

tmp <- (tmp
          %>% SetAssayData(assay.type = "TCRgd", slot = "raw.data", new.data = csv.TCRgd)
	  %>% NormalizeData(assay.type = "TCRgd", normalization.method = "genesCLR")
	  %>% ScaleData(assay.type = "TCRgd", display.progress = TRUE))

print("load TCRgd")

# =========================================================================================
# Import ADT into Seurat Object
# =========================================================================================
csv.adt <- read.csv(ADT, sep = ",", header = TRUE, row.names = 1)
colnames(csv.adt) <- paste(sample_short, colnames(csv.adt), sep = '_')
rownames(csv.adt) <- paste0("CITE_", rownames(csv.adt))

tmp <- (tmp
          %>% SetAssayData(assay.type = "CITE", slot = "raw.data", new.data = csv.adt)
	  %>% NormalizeData(assay.type = "CITE", normalization.method = "genesCLR")
	  %>% ScaleData(assay.type = "CITE", display.progress = TRUE))

print("load ADT")

# =========================================================================================
# Import clonotypeAB into Seurat Object
# =========================================================================================
csv.Cab <- read.csv(Cab, sep = ",", header = TRUE, row.names = 1)
colnames(csv.Cab) <- paste(sample_short, colnames(csv.Cab), sep = '_')
rownames(csv.Cab) <- paste0("clonotypeAB_", rownames(csv.Cab))

tmp <- (tmp
          %>% SetAssayData(assay.type = "clonotypeAB", slot = "raw.data", new.data = csv.Cab)
	  %>% NormalizeData(assay.type = "clonotypeAB", normalization.method = "genesCLR")
	  %>% ScaleData(assay.type = "clonotypeAB", display.progress = TRUE))

print("load clonotypeAB")

# =========================================================================================
# Import clonotypeGD into Seurat Object
# =========================================================================================
csv.Cgd <- read.csv(Cgd, sep = ",", header = TRUE, row.names = 1)
colnames(csv.Cgd) <- paste(sample_short, colnames(csv.Cgd), sep = '_')
rownames(csv.Cgd) <- paste0("clonotypeGD_", rownames(csv.Cgd))

tmp <- (tmp
          %>% SetAssayData(assay.type = "clonotypeGD", slot = "raw.data", new.data = csv.Cgd)
	  %>% NormalizeData(assay.type = "clonotypeGD", normalization.method = "genesCLR")
	  %>% ScaleData(assay.type = "clonotypeGD", display.progress = TRUE))

print("load clonotypeGD")


# =====
# Save
# =====
saveRDS(tmp, f_output)