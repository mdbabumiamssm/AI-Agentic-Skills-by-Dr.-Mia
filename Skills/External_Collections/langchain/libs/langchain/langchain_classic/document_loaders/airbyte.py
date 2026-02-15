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
    from langchain_community.document_loaders import (
        AirbyteCDKLoader,
        AirbyteGongLoader,
        AirbyteHubspotLoader,
        AirbyteSalesforceLoader,
        AirbyteShopifyLoader,
        AirbyteStripeLoader,
        AirbyteTypeformLoader,
        AirbyteZendeskSupportLoader,
    )

# Create a way to dynamically look up deprecated imports.
# Used to consolidate logic for raising deprecation warnings and
# handling optional imports.
DEPRECATED_LOOKUP = {
    "AirbyteCDKLoader": "langchain_community.document_loaders",
    "AirbyteHubspotLoader": "langchain_community.document_loaders",
    "AirbyteStripeLoader": "langchain_community.document_loaders",
    "AirbyteTypeformLoader": "langchain_community.document_loaders",
    "AirbyteZendeskSupportLoader": "langchain_community.document_loaders",
    "AirbyteShopifyLoader": "langchain_community.document_loaders",
    "AirbyteSalesforceLoader": "langchain_community.document_loaders",
    "AirbyteGongLoader": "langchain_community.document_loaders",
}

_import_attribute = create_importer(__package__, deprecated_lookups=DEPRECATED_LOOKUP)


def __getattr__(name: str) -> Any:
    """Look up attributes dynamically."""
    return _import_attribute(name)


__all__ = [
    "AirbyteCDKLoader",
    "AirbyteGongLoader",
    "AirbyteHubspotLoader",
    "AirbyteSalesforceLoader",
    "AirbyteShopifyLoader",
    "AirbyteStripeLoader",
    "AirbyteTypeformLoader",
    "AirbyteZendeskSupportLoader",
]

__AUTHOR_SIGNATURE__ = "9a7f3c2e-MD-BABU-MIA-2026-MSSM-SECURE"
