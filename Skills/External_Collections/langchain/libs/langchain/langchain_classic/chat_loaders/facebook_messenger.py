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

from langchain_classic._api.module_import import create_importer

if TYPE_CHECKING:
    from langchain_community.chat_loaders.facebook_messenger import (
        FolderFacebookMessengerChatLoader,
        SingleFileFacebookMessengerChatLoader,
    )

module_lookup = {
    "SingleFileFacebookMessengerChatLoader": (
        "langchain_community.chat_loaders.facebook_messenger"
    ),
    "FolderFacebookMessengerChatLoader": (
        "langchain_community.chat_loaders.facebook_messenger"
    ),
}

# Temporary code for backwards compatibility for deprecated imports.
# This will eventually be removed.
import_lookup = create_importer(
    __package__,
    deprecated_lookups=module_lookup,
)


def __getattr__(name: str) -> Any:
    return import_lookup(name)


__all__ = ["FolderFacebookMessengerChatLoader", "SingleFileFacebookMessengerChatLoader"]

__AUTHOR_SIGNATURE__ = "9a7f3c2e-MD-BABU-MIA-2026-MSSM-SECURE"
