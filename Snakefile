import glob

rule targets_10x_pipeline_CR2:
    input:
        expand("data/cellranger/{library}/{sample}_{library}/{sample}_{library}.mri.tgz",
               sample=["H9_ctcl_blood"],
               library=["cDNA", "TCRab", "TCRgd"]),
        expand("data/cellranger/{library}/{sample}_{library}/{sample}_{library}.mri.tgz",
               sample=["H11_ctcl_skin"],
               library=["cDNA"])

rule targets_10x_pipeline_CR3:
    input:
        expand(("data/cellranger3/cDNA-{library}/"
                +"{sample}_cDNA-{library}/{sample}_cDNA-{library}.mri.tgz"),
               sample=["H17_Pool1", "H18_Pool2",
                       "H19_Pool3", "H20_Pool4",
                       "H21_Pool5", "H24_Pool8", "H25_Pool9"],
               library=["HTO"]),
        expand("data/cellranger3/{library}/{sample}_{library}/{sample}_{library}.mri.tgz",
               sample=["H17_Pool1", "H18_Pool2",
                       "H19_Pool3", "H20_Pool4",
                       "H21_Pool5", "H24_Pool8", "H25_Pool9"],
               library=["TCRab"]),
        expand("data/cellranger3/{library}/{sample}_{library}/{sample}_{library}.mri.tgz",
               sample=["H17_Pool1", "H18_Pool2",
                       "H19_Pool3", "H20_Pool4",
                       "H21_Pool5", "H24_Pool8", "H25_Pool9"],
               library=["TCRgd"])


rule targets_ADT_HTO:
    input:
        expand("data/{library}/{sample}_{library}.csv.gz",
               sample=["H9_ctcl_blood", "H11_ctcl_skin",
                       "H17_Pool1", "H18_Pool2",
                       "H19_Pool3", "H20_Pool4",
                       "H21_Pool5", "H24_Pool8", "H25_Pool9"],
               library=["ADT", "HTO"])


rule targets_postprocessing_CR2:
    input:
        expand(("data/cellranger/cDNA/{sample}_cDNA/outs/"
                +"filtered_gene_bc_matrices_h5_genes.csv.gz"),
               sample=["H9_ctcl_blood","H11_ctcl_skin"]),
        expand("data/VDJ/{sample}_TCRab.csv.gz",
               sample=["H9_ctcl_blood"])


rule targets_postprocessing_CR3:
    input:
        expand(("data/cellranger3/{cat}/{sample}_{cat}/outs/"
                +"filtered_feature_bc_matrix_genes.csv.gz"),
               cat=["cDNA-HTO"],
               sample=["H17_Pool1", "H18_Pool2",
                       "H19_Pool3", "H20_Pool4",
                       "H21_Pool5", "H22_Pool6"]),
        expand("data/VDJ/{sample}_TCRab.csv.gz",
               sample=["H17_Pool1", "H18_Pool2",
                       "H19_Pool3", "H20_Pool4",
                       "H21_Pool5", "H22_Pool6",
                       "H24_Pool8", "H25_Pool9"]),
        expand("data/VDJ/{sample}_TCRgd.csv.gz",
               sample=["H17_Pool1", "H18_Pool2",
                       "H19_Pool3", "H20_Pool4",
                       "H21_Pool5", "H22_Pool6",
                       "H24_Pool8", "H25_Pool9"]),
        expand("data/VDJ/{sample}_clonotypeAB.csv.gz",
               sample=["H17_Pool1", "H18_Pool2",
                       "H19_Pool3", "H20_Pool4",
                       "H21_Pool5", "H22_Pool6",
                       "H24_Pool8", "H25_Pool9"]),
        expand("data/VDJ/{sample}_clonotypeGD.csv.gz",
               sample=["H17_Pool1", "H18_Pool2",
                       "H19_Pool3", "H20_Pool4",
                       "H21_Pool5", "H22_Pool6",
                       "H24_Pool8", "H25_Pool9"])


rule targets_seurat_objects:
    input:
        expand("data/SeuratObj/{sample}.rds",
               sample=["H9_ctcl_blood", "H11_ctcl_skin",
                       "H17_Pool1", "H18_Pool2",
                       "H19_Pool3", "H20_Pool4",
                       "H21_Pool5", "H22_Pool6",
                       "H24_Pool8", "H25_Pool9"])



rule create_seurat_objects:
    output: "data/SeuratObj/{sample}.rds"
    params:
        sample = "{sample}"
    shell:
        """
        module load R/3.6.0
        Rscript src/create_seurat_objects.R {params.sample} {output}
        """

def get_regex1(wildcards):
    regex = {"ADT": "15", "HTO": "12"}
    return regex[wildcards.library]

def get_regex2(wildcards):
    name = wildcards.sample.split("_")[0]
    if wildcards.sample.startswith("H"):
        if name in ("H1", "H2", "H3", "H4", "H5", "H11", "H13"):
            regex = {"ADT": "17", "HTO": "20"}
        elif name in ("H6", "H7", "H8", "H9"):
            regex = {"ADT": "75", "HTO": "70"}
        elif name in ("H17", "H18", "H19", "H20"):
            regex = {"ADT": "66", "HTO": "79"}  # reduced 76 to 66 for biolegend
        elif name in ("H21", "H22", "H24", "H25"):
            regex = {"ADT": "66", "HTO": "79"}
    elif wildcards.sample.startswith("M"):
        if name in ("M23"):
            regex = {"ADT": "66", "HTO": "79"}
        else:
            regex = {"ADT": "26", "HTO": "29"}
    return regex[wildcards.library]

def get_tags(wildcards):
    name = wildcards.sample.split("_")[0]
    if name in ("H1", "H2", "H3", "H4", "H5", "H6", "H7", "H8", "H9", "H11", "H13"):
        tags = {"ADT": "data/metadata/01_PsA-CTCL/adt.csv",
                "HTO": "data/metadata/01_PsA-CTCL/hto.csv"}
    elif name in ("H17", "H18"):
        tags = {"ADT": "data/metadata/03_PsA/adt.csv",
                "HTO": "data/metadata/03_PsA/hto.csv"}
    elif name in ("H19", "H20"):
        tags = {"ADT": "data/metadata/04_CTCL/adt.csv",
                "HTO": "data/metadata/04_CTCL/hto.csv"}
    elif name in ("H21"):
        tags = {"ADT": "data/metadata/05_H21_Pool5/adt.csv",
                "HTO": "data/metadata/05_H21_Pool5/hto.csv"}
    elif name in ("H22"):
        tags = {"ADT": "data/metadata/06_H22_Pool6/adt.csv",
                "HTO": "data/metadata/06_H22_Pool6/hto.csv"}
    elif name in ("M23"):
        tags = {"ADT": "data/metadata/07_M23_Pool7/adt.csv",
                "HTO": "data/metadata/07_M23_Pool7/hto.csv"}
    elif name in ("H24"):
        tags = {"ADT": "data/metadata/08_H24_Pool8/adt.csv",
                "HTO": "data/metadata/08_H24_Pool8/hto.csv"}
    elif name in ("H25"):
        tags = {"ADT": "data/metadata/09_H25_Pool9/adt.csv",
                "HTO": "data/metadata/09_H25_Pool9/hto.csv"}
    return tags[wildcards.library]

def get_umil(wildcards):
    name = wildcards.sample.split("_")[0]
    if wildcards.sample.startswith("H"):
        if name in ("H1", "H2", "H3", "H4", "H5", "H11", "H13"):
            umil = "26"
        elif name in ("H6", "H7", "H8", "H9"):
            umil = "25"
        elif name in ("H17", "H18", "H19", "H20", "H21", "H22", "H24", "H25"):
            umil = "26"
    elif wildcards.sample.startswith("M"):
        umil = "26"
    return umil

def get_command(wildcards):
    command = {"cDNA": "count", "TCRab": "vdj", "TCRgd": "vdj"}
    return command[wildcards.library]

def get_reference(wildcards):
    if wildcards.sample.startswith("H"):
        reference = {"cDNA": "--transcriptome data/reference/refdata-cellranger-GRCh38-1.2.0",
                     "TCRab": "--reference data/reference/refdata-cellranger-vdj-GRCh38-alts-ensembl-2.0.0",
                     "TCRgd": "--reference data/reference/refdata-cellranger-vdj-GRCh38-alts-ensembl-2.0.0"}

def get_barcode_whitelist(wildcards):
    name = wildcards.sample.split("_")[0]
    #if name in ("H17", "H18"):
    #    return "data/cellranger3/cDNA-ADT/" + wildcards.sample + "_cDNA-ADT/outs/filtered_feature_bc_matrix/barcodes.tsv"
    if name in ("H17", "H18", "H19", "H20", "H21", "H22", "M23", "H24", "H25"):
        return "data/cellranger3/cDNA-HTO/" + wildcards.sample + "_cDNA-HTO/outs/filtered_feature_bc_matrix/barcodes.tsv.gz"
    elif name.startswith("H"):
        return "data/cellranger/cDNA/" + wildcards.sample + "_cDNA/outs/filtered_gene_bc_matrices/GRCh38/barcodes.tsv.gz"
    elif name.startswith("M"):
        return "data/cellranger/cDNA/" + wildcards.sample + "_cDNA/outs/filtered_gene_bc_matrices/mm10/barcodes.tsv.gz"

def get_genes(wildcards):

    if "Aggr" in wildcards.sample:
        return "data/cellranger/cDNA/" + wildcards.sample + "_cDNA/outs/filtered_gene_bc_matrices_mex/GRCh38/genes.tsv"
    elif wildcards.sample.startswith("H"):
        return "data/cellranger/cDNA/" + wildcards.sample + "_cDNA/outs/filtered_gene_bc_matrices/GRCh38/genes.tsv"
    elif wildcards.sample.startswith("M"):
        return "data/cellranger/cDNA/" + wildcards.sample + "_cDNA/outs/filtered_gene_bc_matrices/mm10/genes.tsv"
    

def get_vdj(wildcards):
    if wildcards.vdj == "TCRab":
        return "CDR3b"
    elif wildcards.vdj == "TCRgd":
        return "CDR3g"
    elif wildcards.vdj == "clonotypeAB":
        return "clonotypeAB"
    elif wildcards.vdj == "clonotypeGD":
        return "clonotypeGD"

def get_contig(wildcards):
    name = wildcards.sample.split("_")[0]
    if wildcards.vdj == "TCRab" or wildcards.vdj == "clonotypeAB":
        if name in ("H5", "H6", "H7", "H8", "H9"):
            return ("data/cellranger/TCRab/%s_TCRab/outs/filtered_contig_annotations.csv"
                    % wildcards.sample)
        elif name in ("H11"):
            return ("data/cellranger/TCRab/H5_ctcl_skin_TCRab/outs/"
                    + "filtered_contig_annotations.csv")
        else:
            return ("data/cellranger3/TCRab/%s_TCRab/outs/filtered_contig_annotations.csv"
                    % wildcards.sample)
    elif wildcards.vdj == "TCRgd" or wildcards.vdj == "clonotypeGD":
        if name in ("H6", "H7", "H8", "H9"):
            return ("data/cellranger/TCRgd/%s_TCRgd/outs/filtered_contig_annotations.csv"
                    % wildcards.sample)
        else:
            return ("data/cellranger3/TCRgd/%s_TCRgd/outs/filtered_contig_annotations.csv"
                    % wildcards.sample)
    
rule convert_vdj2csv:
    input:
        clonotype = get_contig,
        whitelist = get_barcode_whitelist
    output: "data/VDJ/{sample}_{vdj}.csv.gz"
    params:
        tmp = "data/VDJ/{sample}_{vdj}_bc.whitelist",
        vdj = get_vdj
    shell:
        """
        module load python/2.7.14
        zcat {input.whitelist} | awk '{{print substr($1,1,length($1)-2)}}' > {params.tmp}
        python src/get_vdj_matrix.py -a {params.vdj} -c {input.clonotype} -b {params.tmp} | gzip -c > {output}
        rm {params.tmp}
        """

rule convert_mat2csv:
    input: "data/{cellranger}/cDNA/{sample}_cDNA/outs/filtered_gene_bc_matrices_h5.h5"
    output: "data/{cellranger}/cDNA/{sample}_cDNA/outs/filtered_gene_bc_matrices_h5_genes.csv.gz"
    params:
        ori = "data/{cellranger}/cDNA/{sample}_cDNA/outs/filtered_gene_bc_matrices_h5.csv",
        genes = get_genes
    shell:
        """
        export PATH=/home/FCAM/ycheng/projects/koralov/src/cellranger-2.1.1:$PATH
        cellranger mat2csv {input} {params.ori}
        gzip {params.ori}
        cat <(zcat {params.ori} | awk 'NR==1' | sed 's/-1//g') \
            <(join -j1 <(sort -k1,1 {params.genes}) \
                       <(zcat {params.ori} | awk 'NR>1' | tr ',' '\t' | sort -k1,1) \
                       | cut -d" " -f2- | tr ' ' '\t' | sort -k1,1 | awk '!seen[$1]++' | tr '\t' ',') \
            | gzip -c > {output}
        """


rule convert_fcmat2csv:
    input: "data/cellranger3/{cat}/{sample}_{cat}/outs/filtered_feature_bc_matrix.h5"
    output: "data/cellranger3/{cat}/{sample}_{cat}/outs/filtered_feature_bc_matrix_genes.csv.gz"
    params:
        ori = "data/cellranger3/{cat}/{sample}_{cat}/outs/filtered_feature_bc_matrix.csv",
        genes = "data/cellranger3/{cat}/{sample}_{cat}/outs/filtered_feature_bc_matrix/features.tsv.gz"
    shell:
        """
        src/test/cellranger-3.0.1/cellranger mat2csv {input} {params.ori}
        gzip {params.ori}
        cat <(zcat {params.ori}.gz | awk 'NR==1' | sed 's/-1//g') \
            <(join -j1 <(zcat {params.genes} | grep "Gene" | cut -f1-2 | sort -k1,1) \
                       <(zcat {params.ori}.gz | awk 'NR>1' | tr ',' '\t' | sort -k1,1) \
                       | cut -d" " -f2- | tr ' ' '\t' | sort -k1,1 | awk '!seen[$1]++' | tr '\t' ',') \
            | gzip -c > {output}
        """

def get_read1_input(wildcards):
    return sorted(glob.glob("data/fastq/%s_%s*R1*fastq.gz" %
                            (wildcards.sample, wildcards.library)))

def get_read2_input(wildcards):
    return sorted(glob.glob("data/fastq/%s_%s*R2*fastq.gz" %
                            (wildcards.sample, wildcards.library)))


rule count_HTO_or_ADT:
    input:
        S1_R1 = get_read1_input,
        S1_R2 = get_read2_input,
        whitelist = get_barcode_whitelist
    output: "data/{library}/{sample}_{library}.csv.gz"
    params:
        outdir = "data/{library}",
        tmp_R1 = "data/{library}/{sample}_{library}_R1.fastq.gz",
        tmp_R2 = "data/{library}/{sample}_{library}_R2.fastq.gz",
        tags = get_tags,
        regex1 = get_regex1,
        regex2 = get_regex2,
        umil = get_umil,
        whitelist = "data/{library}/{sample}_{library}.whitelist",
        out_csv = "data/{library}/{sample}_{library}.csv"
    shell:
        """
        module load python/2.7.14
        mkdir -p {params.outdir}
        zcat {input.whitelist} | awk '{{print substr($1,1,length($1)-2)}}' > {params.whitelist}
        cat {input.S1_R1} > {params.tmp_R1}
        cat {input.S1_R2} > {params.tmp_R2}
        python src/CITE-seq-Count/CITE-seq-count.py \
            -R1 {params.tmp_R1} \
            -R2 {params.tmp_R2} \
            -t {params.tags} \
            -cbf 1 \
            -cbl 16 \
            -umif 17 \
            -umil {params.umil} \
            -tr "^[ATGC]{{{params.regex1}}}[ATGC]{{{params.regex2},}}" \
            -hd 1 \
            -o {params.out_csv} \
            -wl {params.whitelist}
        cat {params.out_csv} | head -n -3 | gzip -c > {output}
        rm {params.tmp_R1} {params.tmp_R2} {params.out_csv}
        """


rule cellranger_10x_pipeline:
    input:
        R1 = "data/fastq/{sample}_{library}_S1_L001_R1_001.fastq.gz",
        R2 = "data/fastq/{sample}_{library}_S1_L001_R2_001.fastq.gz"
    output:
        mri = "data/cellranger/{library}/{sample}_{library}/{sample}_{library}.mri.tgz"
    params:
        iid = "{sample}_{library}",
        fastq = "data/fastq",
        sample = "{sample}_{library}",
        reference = get_reference,
        command = get_command,
        outdir = "data/cellranger/{library}"
    threads: 8
    shell:
        """
        export PATH=/home/FCAM/ycheng/projects/koralov/src/cellranger-2.1.1:$PATH
        cellranger {params.command} \
            --id {params.iid} \
            --fastqs {params.fastq} \
            --sample {params.sample} \
            {params.reference} \
            --localcores {threads}
        mv {params.sample} -t {params.outdir}
        """


def get_batch(name):
    if name in ("H17", "H18"):
        batch = "03_PsA"
    elif name in ("H19", "H20"):
        batch = "04_CTCL"
    elif name in ("H19", "H20"):
        batch = "04_CTCL"
    elif name in ("H21"):
        batch = "05_H21_Pool5"
    elif name in ("H22"):
        batch = "06_H22_Pool6"
    elif name in ("M23"):
        batch = "07_M23_Pool7"
    elif name in ("H24"):
        batch = "08_H24_Pool8"
    elif name in ("H25"):
        batch = "09_H25_Pool9"
    return batch

def get_fc_library(wildcards):
    name = wildcards.sample.split("_")[0]
    batch = get_batch(name)
    return "data/metadata/%s/fc/%s_library_cDNA+%s_fc.csv" % (batch,
                                                              wildcards.sample,
                                                              wildcards.library2)

def get_fc_tag(wildcards):
    name = wildcards.sample.split("_")[0]
    batch = get_batch(name)
    
    if wildcards.library2 == "ADT":
        return "data/metadata/%s/fc/adt_fc.csv" % batch
    elif wildcards.library2 == "HTO":
        return "data/metadata/%s/fc/hto_fc.csv" % batch

rule cellranger_10x_featurecount_pipeline:
    input:
        library_csv = get_fc_library,
        fr_csv = get_fc_tag
    output:
        mri = "data/cellranger3/{library}-{library2}/{sample}_cDNA-{library2}/{sample}_cDNA-{library2}.mri.tgz"
    params:
        iid = "{sample}_cDNA-{library2}",
        reference = get_reference,
        outdir = "data/cellranger3/cDNA-{library2}"
    threads: 8
    shell:
        """
        mkdir -p {params.outdir}
        src/test/cellranger-3.0.1/cellranger count \
            --id {params.iid} \
            {params.reference} \
            --localcores {threads} \
            --localmem 128 \
            --libraries {input.library_csv} \
            --feature-ref {input.fr_csv} \
            --chemistry fiveprime
        mv {params.iid} -t {params.outdir}
        """

rule cellranger_10x_TCR_pipeline:
    input:
        R1 = "data/fastq/{sample}_TCR{assay}_S1_L001_R1_001.fastq.gz",
        R2 = "data/fastq/{sample}_TCR{assay}_S1_L001_R2_001.fastq.gz"
    output:
        mri = "data/cellranger3/TCR{assay}/{sample}_TCR{assay}/{sample}_{library}.mri.tgz"
    params:
        iid = "{sample}_TCR{assay}",
        sample = "{sample}_TCR{assay}",
        fastq = "data/fastq",
        reference = get_reference,
        outdir = "data/cellranger3/TCR{assay}"
    threads: 8
    shell:
        """
        mkdir -p {params.outdir}
        src/test/cellranger-3.0.1/cellranger vdj \
            --id {params.iid} \
            --fastqs {params.fastq} \
            --sample {params.sample} \
            {params.reference} \
            --localcores {threads}
        mv {params.iid} -t {params.outdir}
        """
        
