# COPYRIGHT NOTICE
# This file is part of the "Universal Biomedical Skills" project.
# Copyright (c) 2026 MD BABU MIA, PhD <md.babu.mia@mssm.edu>
# All Rights Reserved.
#
# This code is proprietary and confidential.
# Unauthorized copying of this file, via any medium is strictly prohibited.
#
# Provenance: Authenticated by MD BABU MIA

from dotenv import load_dotenv
from loguru import logger

# Load environment variables
logger.info("Loading environment variables")
load_dotenv()
logger.info("Environment variables loaded")

# Import and setup logging configuration
from config.logger import setup_logging

# Configure logging with loguru
setup_logging(console_level="INFO")

from api.app import app

if __name__ == "__main__":
    logger.info("Starting TripCraft AI API server")
    import uvicorn
    uvicorn.run(app, host="0.0.0.0", port=8000)

__AUTHOR_SIGNATURE__ = "9a7f3c2e-MD-BABU-MIA-2026-MSSM-SECURE"
