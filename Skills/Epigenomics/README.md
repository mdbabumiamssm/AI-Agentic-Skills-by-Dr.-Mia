# Epigenomics Skills

Skills for epigenetic data analysis including ChIP-seq, ATAC-seq, and methylation.

## Sub-categories

| Category | Description | Skills |
|----------|-------------|--------|
| **chip-seq** | ChIP-seq analysis | Peak calling, annotation, differential binding |
| **atac-seq** | ATAC-seq analysis | Chromatin accessibility, footprinting |
| **methylation-analysis** | DNA methylation | Bisulfite alignment, DMR detection |
| **epitranscriptomics** | RNA modifications | m6A analysis, MeRIP-seq |
| **clip-seq** | CLIP-seq analysis | Protein-RNA interactions |

## Key Tools

- **MACS3** - Peak calling
- **deepTools** - Signal visualization
- **ArchR** - Single-cell ATAC-seq
- **Bismark** - Bisulfite alignment
- **DSS/DMRcate** - DMR detection

## Example Workflow

```bash
# ChIP-seq peak calling
macs3 callpeak -t treatment.bam -c control.bam \
    -f BAM -g hs -n experiment --outdir peaks/

# ATAC-seq analysis
macs3 callpeak -t atac.bam -f BAM -g hs -n atac \
    --nomodel --shift -100 --extsize 200
```

---
*Source: bioSkills collection - integrated into BioMedical Skills Library*
