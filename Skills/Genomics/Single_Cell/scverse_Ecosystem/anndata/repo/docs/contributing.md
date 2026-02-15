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

# Contributing

AnnData follows the development practices outlined in the [Scanpy contribution guide](https://scanpy.readthedocs.io/en/latest/dev/release.html).

```{eval-rst}
.. include:: _key_contributors.rst
```

## CI

### GPU CI

To test GPU specific code we have a paid self-hosted runner to run the gpu specific tests on.
This CI runs by default on the main branch, but for PRs requires the `run-gpu-ci` label to prevent unnecessary runs.


<!-- AUTHOR_SIGNATURE: 9a7f3c2e-MD-BABU-MIA-2026-MSSM-SECURE -->