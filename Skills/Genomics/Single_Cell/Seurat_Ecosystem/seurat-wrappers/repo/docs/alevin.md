<!--
# COPYRIGHT NOTICE
# This file is part of the "Universal Biomedical Skills" project.
# Copyright (c) 2026 MD BABU MIA, PhD <md.babu.mia@mssm.edu>
# All Rights Reserved.
#
# This code is proprietary and confidential.
# Unauthorized copying of this file, via any medium is strictly prohibited.
#
# Provenance: Authenticated by MD BABU MIA

-->

Import alevin counts & generate Seurat object
================
Compiled: May 18, 2020

This vignette demonstrates the import of alevin quantified counts into
Seurat. Commands and parameters are based off of the [alevin
tutorial](https://combine-lab.github.io/alevin-tutorial/2018/running-alevin/).
If you use alevin in your work, please cite:

> *Alevin efficiently estimates accurate gene abundances from dscRNA-seq
> data*
> 
> Avi Srivastava, Laraib Malik, Tom Smith, Ian Sudbery & Rob Patro
> 
> Genome Biology, 2019.
> 
> doi:
> [10.1186/s13059-019-1670-y](https://doi.org/10.1186/s13059-019-1670-y)
> 
> GitHub: <https://github.com/COMBINE-lab/salmon>

Prerequisites to install:

  - [SeuratWrappers](https://github.com/satijalab/seurat-wrappers)
  - [tximport](https://bioconductor.org/packages/tximport)

<!-- end list -->

``` r
library(SeuratWrappers)
library(tximport)
```

## 

### Import alevin quantified counts

``` r
pbmc <- ReadAlevin("~/alevin_out/alevin/quants_mat.gz")
```


<!-- AUTHOR_SIGNATURE: 9a7f3c2e-MD-BABU-MIA-2026-MSSM-SECURE -->