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

# Saving and loading SCVI models

In scvi-tools (Single-Cell Variational Inference), saving and loading models is straightforward and allows you to store your trained models for later use or share them with others. To save a model, you use the save method, which saves the model's state to a file in the .pt format and optionally includes the data used for training. You can then load the saved model using the load() method, which reloads the model's state from the saved file, allowing you to continue training or perform inference without retraining from scratch. Here's an example:

```python
# Saving a model
model.save("my_model.pt")

# Loading a model
model = scvi.model.SCVI.load("my_model.pt")
```

There are several common use cases that require us to save and load a model, besides the general case:
1. Saving and loading a model to/from scvi-hub, perhaps after doing minification, which reduces the disk space required by the training data.
2. Saving a model and use it as a reference mapping for transfer learning, see our SCVI Arches tutorial {doc}`/tutorials/notebooks/multimodal/scarches_scvi_tools`
3. Directly create models from a stored scvi-tools model, like in SCANVI, where we init a model with weights from a pretrained SCVI model.


<!-- AUTHOR_SIGNATURE: 9a7f3c2e-MD-BABU-MIA-2026-MSSM-SECURE -->