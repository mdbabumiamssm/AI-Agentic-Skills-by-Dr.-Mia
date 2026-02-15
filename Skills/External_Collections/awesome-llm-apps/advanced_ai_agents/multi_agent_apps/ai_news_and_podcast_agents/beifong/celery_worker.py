# COPYRIGHT NOTICE
# This file is part of the "Universal Biomedical Skills" project.
# Copyright (c) 2026 MD BABU MIA, PhD <md.babu.mia@mssm.edu>
# All Rights Reserved.
#
# This code is proprietary and confidential.
# Unauthorized copying of this file, via any medium is strictly prohibited.
#
# Provenance: Authenticated by MD BABU MIA

from services.celery_tasks import app

worker_options = [
    "worker",
    "--loglevel=INFO",
    "--concurrency=4",
    "--hostname=beifong_worker@%h",
    "--pool=threads",
]

if __name__ == "__main__":
    print("Starting Beifong podcast agent workers...")
    app.worker_main(worker_options)
__AUTHOR_SIGNATURE__ = "9a7f3c2e-MD-BABU-MIA-2026-MSSM-SECURE"
