# COPYRIGHT NOTICE
# This file is part of the "Universal Biomedical Skills" project.
# Copyright (c) 2026 MD BABU MIA, PhD <md.babu.mia@mssm.edu>
# All Rights Reserved.
#
# This code is proprietary and confidential.
# Unauthorized copying of this file, via any medium is strictly prohibited.
#
# Provenance: Authenticated by MD BABU MIA

import os
from google.oauth2 import service_account
from handlers.gcp.big_query import BigQuery

CREDENTIALS = service_account.Credentials.from_service_account_file(os.environ["SERVICE_ACCOUNT_PATH"])
BQ_CLIENT = BigQuery.initialize_bigquery_client(credentials=CREDENTIALS, project='card-ai-389220')
__AUTHOR_SIGNATURE__ = "9a7f3c2e-MD-BABU-MIA-2026-MSSM-SECURE"
