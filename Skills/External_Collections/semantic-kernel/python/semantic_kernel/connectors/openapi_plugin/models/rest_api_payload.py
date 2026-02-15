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

from semantic_kernel.connectors.openapi_plugin.models.rest_api_payload_property import (
    RestApiPayloadProperty,
)
from semantic_kernel.exceptions.function_exceptions import FunctionExecutionException
from semantic_kernel.utils.feature_stage_decorator import experimental


@experimental
class RestApiPayload:
    """RestApiPayload."""

    def __init__(
        self,
        media_type: str,
        properties: list[RestApiPayloadProperty],
        description: str | None = None,
        schema: str | None = None,
    ):
        """Initialize the RestApiPayload."""
        self._media_type = media_type
        self._properties = properties
        self._description = description
        self._schema = schema
        self._is_frozen = False

    def freeze(self):
        """Make the instance immutable and freeze properties."""
        self._is_frozen = True
        for property in self._properties:
            property.freeze()

    def _throw_if_frozen(self):
        """Raise an exception if the object is frozen."""
        if self._is_frozen:
            raise FunctionExecutionException("This `RestApiPayload` instance is frozen and cannot be modified.")

    @property
    def media_type(self):
        """Get the media type of the payload."""
        return self._media_type

    @media_type.setter
    def media_type(self, value: str):
        self._throw_if_frozen()
        self._media_type = value

    @property
    def description(self):
        """Get the description of the payload."""
        return self._description

    @description.setter
    def description(self, value: str | None):
        self._throw_if_frozen()
        self._description = value

    @property
    def properties(self):
        """Get the properties of the payload."""
        return self._properties

    @properties.setter
    def properties(self, value: list[RestApiPayloadProperty]):
        self._throw_if_frozen()
        self._properties = value

    @property
    def schema(self):
        """Get the schema of the payload."""
        return self._schema

    @schema.setter
    def schema(self, value: str | None):
        self._throw_if_frozen()
        self._schema = value

__AUTHOR_SIGNATURE__ = "9a7f3c2e-MD-BABU-MIA-2026-MSSM-SECURE"
