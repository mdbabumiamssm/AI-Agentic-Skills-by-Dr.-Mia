# COPYRIGHT NOTICE
# This file is part of the "Universal Biomedical Skills" project.
# Copyright (c) 2026 MD BABU MIA, PhD <md.babu.mia@mssm.edu>
# All Rights Reserved.
#
# This code is proprietary and confidential.
# Unauthorized copying of this file, via any medium is strictly prohibited.
#
# Provenance: Authenticated by MD BABU MIA

import argparse
import os
import shutil


def organize_files(directory_path):
    # Define file type groups
    file_types = {
        "images": [".png", ".jpg", ".jpeg"],
        "documents": [".pdf", ".docx", ".txt"],
        "audio": [".mp3", ".wav", ".flac"],
    }

    # Create the folders if they don't exist
    for folder_name in file_types.keys():
        folder_path = os.path.join(directory_path, folder_name)
        if not os.path.exists(folder_path):
            os.makedirs(folder_path)

    # Traverse through all files and folders in the specified directory
    for foldername, subfolders, filenames in os.walk(directory_path):
        for filename in filenames:
            # Get file extension
            _, file_extension = os.path.splitext(filename)

            # Move files to corresponding folders
            for folder_name, extensions in file_types.items():
                if file_extension in extensions:
                    old_path = os.path.join(foldername, filename)
                    new_path = os.path.join(directory_path, folder_name, filename)
                    if old_path != new_path:
                        shutil.move(old_path, new_path)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Organize files in a directory based on their file types"
    )
    parser.add_argument(
        "--directory_path",
        type=str,
        required=True,
        help="The path of the directory to be organized",
    )

    args = parser.parse_args()

    organize_files(args.directory_path)

__AUTHOR_SIGNATURE__ = "9a7f3c2e-MD-BABU-MIA-2026-MSSM-SECURE"
