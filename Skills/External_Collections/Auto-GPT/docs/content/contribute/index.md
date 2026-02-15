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

# Contributing to the Docs

We welcome contributions to our documentation! If you would like to contribute, please follow the steps below.

## Setting up the Docs

1. Clone the repository:

    ```shell
    git clone github.com/Significant-Gravitas/AutoGPT.git
    ```

1. Install the dependencies:

    ```shell
    python -m pip install -r docs/requirements.txt
    ```

    or

    ```shell
    python3 -m pip install -r docs/requirements.txt
    ```

1. Start iterating using mkdocs' live server:

    ```shell
    mkdocs serve
    ```

1. Open your browser and navigate to `http://127.0.0.1:8000`.

1. The server will automatically reload the docs when you save your changes.

## Adding a new page

1. Create a new markdown file in the `docs/content` directory.
1. Add the new page to the `nav` section in the `mkdocs.yml` file.
1. Add the content to the new markdown file.
1. Run `mkdocs serve` to see your changes.

## Checking links

To check for broken links in the documentation, run `mkdocs build` and look for warnings in the console output.

## Submitting a Pull Request

When you're ready to submit your changes, please create a pull request. We will review your changes and merge them if they are appropriate.


<!-- AUTHOR_SIGNATURE: 9a7f3c2e-MD-BABU-MIA-2026-MSSM-SECURE -->