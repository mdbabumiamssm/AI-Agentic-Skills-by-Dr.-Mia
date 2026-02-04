# Workflow Management Skills

Skills for building and managing bioinformatics pipelines.

## Sub-categories

| Category | Description | Skills |
|----------|-------------|--------|
| **workflow-management** | Pipeline tools | Snakemake, Nextflow, CWL, WDL |
| **workflows** | Complete pipelines | RNA-seq, variant calling, scRNA-seq |

## Key Tools

- **Snakemake** - Python-based workflow manager
- **Nextflow** - DSL for data-driven pipelines
- **CWL/cwltool** - Common Workflow Language
- **Cromwell** - WDL execution engine

## Pre-built Workflows

The `workflows` directory contains 35 ready-to-use pipelines:

- RNA-seq differential expression
- Whole genome/exome variant calling
- Single-cell RNA-seq analysis
- ChIP-seq peak calling
- ATAC-seq accessibility
- Metagenomics classification
- And many more...

## Example Snakemake

```python
# Snakefile
SAMPLES = ["sample1", "sample2", "sample3"]

rule all:
    input:
        expand("results/{sample}_aligned.bam", sample=SAMPLES)

rule align:
    input:
        r1 = "data/{sample}_R1.fastq.gz",
        r2 = "data/{sample}_R2.fastq.gz"
    output:
        "results/{sample}_aligned.bam"
    threads: 8
    shell:
        "bwa mem -t {threads} reference.fa {input.r1} {input.r2} | "
        "samtools sort -o {output}"
```

---
*Source: bioSkills collection - integrated into BioMedical Skills Library*
