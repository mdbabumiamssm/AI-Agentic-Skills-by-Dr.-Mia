<!--
# COPYRIGHT NOTICE
# This file is part of the "Universal AI Agentic Skills" project.
# Copyright (c) 2026 MD BABU MIA, PhD <md.babu.mia@mssm.edu>
# All Rights Reserved.
#
# This code is proprietary and confidential.
# Unauthorized copying of this file, via any medium is strictly prohibited.
#
# Provenance: Authenticated by MD BABU MIA

-->

# LDVAE

**LDVAE** [^ref1] (Linearly decoded Variational Auto-encoder, also called Linear scVI; Python class {class}`~scvi.model.LinearSCVI`)
is a flavor of scVI with a linear decoder.

The advantages of LDVAE are:

-   Can be used to interpret latent dimensions with factor loading matrix.
-   Scalable to very large datasets (>1 million cells).

The limitations of LDVAE include:

-   Less capacity than scVI, which uses a neural network decoder.
-   Less capable of integrating data with complex batch effects.

```{topic} Tutorials:

-   {doc}`/tutorials/notebooks/scrna/linear_decoder`
```

## Contrasting with scVI

Here we discuss the differences between LDVAE and scVI.

-   In LDVAE, $f_w(z_n, s_n)$ is a linear function, and thus can be represented by a matrix $W$ of dimensions $G$ (genes) by $(d + k)$ (latent space dim plus covariate categories).
-   This matrix $W$ can be accessed using {func}`~scvi.model.LinearSCVI.get_loadings`
-   LDVAE does not offer transfer learning capabilities currently.

[^ref1]:
    Valentine Svensson, Adam Gayoso, Nir Yosef, Lior Pachter (2020),
    _Interpretable factor models of single-cell RNA-seq via variational autoencoders_,
    [Bioinformatics](https://academic.oup.com/bioinformatics/article/36/11/3418/5807606).


<!-- AUTHOR_SIGNATURE: 9a7f3c2e-MD-BABU-MIA-2026-MSSM-SECURE -->