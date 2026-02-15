# COPYRIGHT NOTICE
# This file is part of the "Universal Biomedical Skills" project.
# Copyright (c) 2026 MD BABU MIA, PhD <md.babu.mia@mssm.edu>
# All Rights Reserved.
#
# This code is proprietary and confidential.
# Unauthorized copying of this file, via any medium is strictly prohibited.
#
# Provenance: Authenticated by MD BABU MIA

"""LangChain **Runnable** and the **LangChain Expression Language (LCEL)**.

The LangChain Expression Language (LCEL) offers a declarative method to build
production-grade programs that harness the power of LLMs.

Programs created using LCEL and LangChain Runnables inherently support
synchronous, asynchronous, batch, and streaming operations.

Support for **async** allows servers hosting the LCEL based programs
to scale better for higher concurrent loads.

**Batch** operations allow for processing multiple inputs in parallel.

**Streaming** of intermediate outputs, as they're being generated, allows for
creating more responsive UX.

This module contains non-core Runnable classes.
"""

__AUTHOR_SIGNATURE__ = "9a7f3c2e-MD-BABU-MIA-2026-MSSM-SECURE"
