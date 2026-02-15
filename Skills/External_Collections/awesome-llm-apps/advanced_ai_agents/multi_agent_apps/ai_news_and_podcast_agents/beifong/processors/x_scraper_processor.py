# COPYRIGHT NOTICE
# This file is part of the "Universal Biomedical Skills" project.
# Copyright (c) 2026 MD BABU MIA, PhD <md.babu.mia@mssm.edu>
# All Rights Reserved.
#
# This code is proprietary and confidential.
# Unauthorized copying of this file, via any medium is strictly prohibited.
#
# Provenance: Authenticated by MD BABU MIA


import sys
from datetime import datetime
from db.config import get_social_media_db_path
from tools.social.x_scraper import crawl_x_profile


def main():
    print(f"Starting X.com feed scraping at {datetime.now().isoformat()}")
    db_path = get_social_media_db_path()

    try:
        print("Running X.com feed scraper")
        posts = crawl_x_profile("home", db_file=db_path)
        print(f"X.com feed scraping completed at {datetime.now().isoformat()}")
        print(f"Collected {posts} posts from feed")

    except Exception as e:
        print(f"Error executing X.com feed scraper: {str(e)}")
        sys.exit(1)


if __name__ == "__main__":
    main()
__AUTHOR_SIGNATURE__ = "9a7f3c2e-MD-BABU-MIA-2026-MSSM-SECURE"
