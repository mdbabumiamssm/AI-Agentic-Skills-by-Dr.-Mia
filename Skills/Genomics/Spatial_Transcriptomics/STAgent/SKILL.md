---
name: st_agent
description: Multimodal AI agent for spatial transcriptomics analysis (H&E + Gene Expression).
license: MIT
metadata:
  author: Artificial Intelligence Group
  version: "1.0.0"
compatibility:
  - system: Python 3.9+
  - gpu: Recommended for image analysis
allowed-tools:
  - run_shell_command
  - read_file
  - write_file
---

# STAgent (Spatial Transcriptomics Agent)

STAgent is a specialized multimodal agent designed to analyze Spatial Transcriptomics (ST) data. It bridges the gap between histology (visual tissue structure) and genomics (gene expression), allowing for spatially-aware analysis.

## When to Use This Skill

*   **Tissue Segmentation**: When you need to annotate tissue regions (e.g., tumor vs. stroma) based on H&E images.
*   **Spatially Variable Genes**: To identify genes with distinct spatial patterns.
*   **Cell-Cell Communication**: To infer interactions between neighboring cell types in a spatial context.
*   **Automated Analysis**: When you want to run a Scanpy/Squidpy workflow via natural language commands.

## Core Capabilities

1.  **Visual Reasoning**: Interprets H&E histology images to guide genomic analysis.
2.  **Code Generation**: Writes and executes Python code (using Scanpy, Squidpy) to answer queries.
3.  **Report Generation**: Produces summaries of spatial domains and gene markers.

## Workflow

1.  **Input Data**: The user provides an H&E image and a gene expression matrix (e.g., from 10x Visium or Xenium).
2.  **Intent Parsing**: The agent determines the goal (e.g., "Find tumor boundary").
3.  **Execution**: The agent runs the necessary analysis pipeline.
4.  **Output**: Returns a textual summary and generates visualization artifacts.

## Example Usage

**User**: "Analyze the tumor regions in sample X and find top spatial markers."

**Agent Action**:
```bash
python3 Skills/Genomics/Spatial_Transcriptomics/STAgent/st_agent_wrapper.py \
  --image_path "./data/sample_he.jpg" \
  --h5ad_path "./data/sample_counts.h5ad" \
  --task "identify_tumor_regions"
```

## Step-by-Step Guide

1.  **Prepare Data**: Ensure your spatial data is in standard formats (TIFF/JPEG for images, H5AD for expression).
2.  **Select Task**: Choose from `identify_tumor_regions`, `find_spatial_genes`, or `cell_cell_interactions`.
3.  **Run Agent**: Execute the wrapper script with appropriate paths.
4.  **Review Output**: Check the generated `analysis_report.md` and plots.
