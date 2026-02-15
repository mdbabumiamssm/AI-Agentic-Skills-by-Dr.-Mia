# COPYRIGHT NOTICE
# This file is part of the "Universal Biomedical Skills" project.
# Copyright (c) 2026 MD BABU MIA, PhD <md.babu.mia@mssm.edu>
# All Rights Reserved.
#
# This code is proprietary and confidential.
# Unauthorized copying of this file, via any medium is strictly prohibited.
#
# Provenance: Authenticated by MD BABU MIA

from typing import TYPE_CHECKING, Any

from langchain_classic._api import create_importer

if TYPE_CHECKING:
    from langchain_community.llms import AzureMLOnlineEndpoint
    from langchain_community.llms.azureml_endpoint import (
        AzureMLEndpointClient,
        ContentFormatterBase,
        CustomOpenAIContentFormatter,
        DollyContentFormatter,
        GPT2ContentFormatter,
        HFContentFormatter,
        OSSContentFormatter,
    )

# Create a way to dynamically look up deprecated imports.
# Used to consolidate logic for raising deprecation warnings and
# handling optional imports.
DEPRECATED_LOOKUP = {
    "AzureMLEndpointClient": "langchain_community.llms.azureml_endpoint",
    "ContentFormatterBase": "langchain_community.llms.azureml_endpoint",
    "GPT2ContentFormatter": "langchain_community.llms.azureml_endpoint",
    "OSSContentFormatter": "langchain_community.llms.azureml_endpoint",
    "HFContentFormatter": "langchain_community.llms.azureml_endpoint",
    "DollyContentFormatter": "langchain_community.llms.azureml_endpoint",
    "CustomOpenAIContentFormatter": "langchain_community.llms.azureml_endpoint",
    "AzureMLOnlineEndpoint": "langchain_community.llms",
}

_import_attribute = create_importer(__package__, deprecated_lookups=DEPRECATED_LOOKUP)


def __getattr__(name: str) -> Any:
    """Look up attributes dynamically."""
    return _import_attribute(name)


__all__ = [
    "AzureMLEndpointClient",
    "AzureMLOnlineEndpoint",
    "ContentFormatterBase",
    "CustomOpenAIContentFormatter",
    "DollyContentFormatter",
    "GPT2ContentFormatter",
    "HFContentFormatter",
    "OSSContentFormatter",
]

__AUTHOR_SIGNATURE__ = "9a7f3c2e-MD-BABU-MIA-2026-MSSM-SECURE"
