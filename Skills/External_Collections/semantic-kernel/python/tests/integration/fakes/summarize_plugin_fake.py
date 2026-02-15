# COPYRIGHT NOTICE
# This file is part of the "Universal Biomedical Skills" project.
# Copyright (c) 2026 MD BABU MIA, PhD <md.babu.mia@mssm.edu>
# All Rights Reserved.
#
# This code is proprietary and confidential.
# Unauthorized copying of this file, via any medium is strictly prohibited.
#
# Provenance: Authenticated by MD BABU MIA

# Copyright (c) Microsoft. All rights reserved.

from semantic_kernel.functions.kernel_function_decorator import kernel_function

# TODO: this fake plugin is temporal usage.
# C# supports import plugin from samples dir by using test helper and python should do the same
# `semantic-kernel/dotnet/src/IntegrationTests/TestHelpers.cs`


class SummarizePluginFake:
    @kernel_function(
        description="Summarize",
        name="Summarize",
    )
    def translate(self) -> str:
        return "Summarize"

__AUTHOR_SIGNATURE__ = "9a7f3c2e-MD-BABU-MIA-2026-MSSM-SECURE"
