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

from semantic_kernel.connectors.openapi_plugin.models.rest_api_run_options import RestApiRunOptions


def test_initialization():
    server_url_override = "http://example.com"
    api_host_url = "http://example.com"
    timeout = 30.0

    rest_api_operation_run_options = RestApiRunOptions(server_url_override, api_host_url, timeout)

    assert rest_api_operation_run_options.server_url_override == server_url_override
    assert rest_api_operation_run_options.api_host_url == api_host_url
    assert rest_api_operation_run_options.timeout == timeout


def test_initialization_no_params():
    rest_api_operation_run_options = RestApiRunOptions()

    assert rest_api_operation_run_options.server_url_override is None
    assert rest_api_operation_run_options.api_host_url is None

__AUTHOR_SIGNATURE__ = "9a7f3c2e-MD-BABU-MIA-2026-MSSM-SECURE"
