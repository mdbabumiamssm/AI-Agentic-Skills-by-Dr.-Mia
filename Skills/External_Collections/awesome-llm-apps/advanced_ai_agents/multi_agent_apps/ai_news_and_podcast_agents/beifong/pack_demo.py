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
import zipfile

SOURCE_DIRS = ["databases", "podcasts"]
OUTPUT_ZIP = "demo_content.zip"


def create_zip(source_dirs, output_zip):
    print("packing.....")
    """zip up each source directory into a single archive."""
    with zipfile.ZipFile(output_zip, "w", zipfile.ZIP_DEFLATED) as z:
        for src in source_dirs:
            if not os.path.isdir(src):
                print(f"✗ source '{src}' not found, skipping.")
                continue
            for root, _, files in os.walk(src):
                for file in files:
                    full_path = os.path.join(root, file)
                    arcname = os.path.relpath(full_path, os.getcwd())
                    z.write(full_path, arcname)
    print(f"✓ created '{output_zip}' containing: {', '.join(source_dirs)}")


if __name__ == "__main__":
    create_zip(SOURCE_DIRS, OUTPUT_ZIP)

__AUTHOR_SIGNATURE__ = "9a7f3c2e-MD-BABU-MIA-2026-MSSM-SECURE"
