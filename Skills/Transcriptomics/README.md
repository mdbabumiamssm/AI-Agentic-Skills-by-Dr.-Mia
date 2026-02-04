# Transcriptomics Skills

Comprehensive skills for RNA-seq analysis, from quantification to differential expression.

## Sub-categories

| Category | Description | Skills |
|----------|-------------|--------|
| **differential-expression** | DE analysis | DESeq2, edgeR, limma |
| **rna-quantification** | Expression quantification | featureCounts, Salmon, kallisto |
| **expression-matrix** | Count matrix handling | Gene ID mapping, normalization |
| **alternative-splicing** | Splicing analysis | rMATS, SUPPA2, LeafCutter |
| **ribo-seq** | Ribosome profiling | Translation efficiency, ORF detection |
| **small-rna-seq** | Small RNA analysis | miRNA/piRNA quantification |

## Key Tools

- **DESeq2/edgeR** - Differential expression
- **Salmon/kallisto** - Pseudo-alignment quantification
- **featureCounts** - Read counting
- **rMATS** - Splicing quantification

## Example DE Analysis

```r
library(DESeq2)

# Create DESeq2 object
dds <- DESeqDataSetFromMatrix(
    countData = counts,
    colData = metadata,
    design = ~ condition
)

# Run analysis
dds <- DESeq(dds)
results <- results(dds, contrast = c("condition", "treated", "control"))
```

---
*Source: bioSkills collection - integrated into BioMedical Skills Library*
