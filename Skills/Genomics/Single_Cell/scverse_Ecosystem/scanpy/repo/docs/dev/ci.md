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

# CI

## Plotting tests

A frequent frustration in testing is the reproducibility of the plots and `matplotlib`'s behaviour in different environments.
We have some tooling to help with this.

### Viewing plots from failed tests on Azure pipelines

The fixtures `check_same_image` and `image_comparer` upload plots from failing tests so you can view them from the azure pipelines test viewer.
To find these, navigate to the tests tab for your build

```{image} ../_static/img/ci_plot-view_tests-tab.png
:width: 750px
```

Select your failing test

```{image} ../_static/img/ci_plot-view_select-test.png
:width: 750px
```

And open the attachments tab

```{image} ../_static/img/ci_plot-view_attachment-tab.png
:width: 750px
```

From here you can view and download the images which were compared, as well as a diff between them.

### Misc

{func}`matplotlib.testing.setup` tries to establish a consistent environment for creating plots. Make sure it's active!


<!-- AUTHOR_SIGNATURE: 9a7f3c2e-MD-BABU-MIA-2026-MSSM-SECURE -->