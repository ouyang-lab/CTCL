# Pipeline and miscellaneous scripts for processing CTCL ECCITE-seq data
The pipeline was created using Snakemake

## Software Pre-requisites
  * [Cell Ranger] (v 3.0.1) - Process 10X 5' GEX and VDJ data for alignment, cell barcode identification, read de-duplication using UMI barcodes, tabulation gene counts and clonotype assembly
  * [CITE-seq-Count] (v 1.3.2) - Pre-processing of ADT and HTO data
  * [Seurat] (v 3.0.1) - Downstream single cell analyses
  * HTODemux - Addition function under the Seurat package to demultiplex pooled samples based on the their assigned HTO barcode
  * [Scrublet] (v 0.2.1) - Intra-sample doublet detection
  * JointVis (v 0.1) - Multi-modal joint analysis using network fusion techniques


  [Cell Ranger]: <https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/installation>
  [CITE-seq-Count]: <https://github.com/Hoohm/CITE-seq-Count>
  [Seurat]: <https://satijalab.org/seurat/install.html>
  [Scrublet]: <https://github.com/AllonKleinLab/scrublet>