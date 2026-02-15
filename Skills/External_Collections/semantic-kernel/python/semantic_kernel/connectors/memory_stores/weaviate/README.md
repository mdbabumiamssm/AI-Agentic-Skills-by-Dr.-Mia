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

# Weaviate Memory Connector

[Weaviate](https://weaviate.io/developers/weaviate) is an open source vector database. Semantic Kernel provides a connector to allow you to store and retrieve information for you AI applications from a Weaviate database.

## Setup

There are a few ways you can deploy your Weaviate database:
- [Weaviate Cloud](https://weaviate.io/developers/weaviate/installation/weaviate-cloud-services)
- [Docker](https://weaviate.io/developers/weaviate/installation/docker-compose)
- [Embedded](https://weaviate.io/developers/weaviate/installation/embedded)
- Other cloud providers such as [Azure](https://azuremarketplace.microsoft.com/en-us/marketplace/apps/weaviatebv1686614539420.weaviate_1?tab=Overview), [AWS](https://weaviate.io/developers/weaviate/installation/aws-marketplace) or [GCP](https://weaviate.io/developers/weaviate/installation/gc-marketplace).

> Note that embedded mode is not supported on Windows yet: [GitHub issue](https://github.com/weaviate/weaviate/issues/3315) and it's still an experimental feature on Linux and MacOS.

## Using the Connector

Once the Weaviate database is up and running, and the environment variables are set, you can use the connector in your Semantic Kernel application. Please refer to this sample to see how to use the connector: [Complex Connector Sample](../../../../samples/concepts/memory/complex_memory.py)

<!-- AUTHOR_SIGNATURE: 9a7f3c2e-MD-BABU-MIA-2026-MSSM-SECURE -->