# COPYRIGHT NOTICE
# This file is part of the "Universal Biomedical Skills" project.
# Copyright (c) 2026 MD BABU MIA, PhD <md.babu.mia@mssm.edu>
# All Rights Reserved.
#
# This code is proprietary and confidential.
# Unauthorized copying of this file, via any medium is strictly prohibited.
#
# Provenance: Authenticated by MD BABU MIA

import re

import requests


def validate_url(url: str, error_format: bool = False, error_response: bool = False) -> bool:
    """Validates a URL.

    Source: https://stackoverflow.com/questions/7160737/how-to-validate-a-url-in-python-malformed-or-not
    """
    regex = re.compile(
        r"^(?:http|ftp)s?://"  # http:// or https://
        r"(?:(?:[A-Z0-9](?:[A-Z0-9-]{0,61}[A-Z0-9])?\.)+(?:[A-Z]{2,6}\.?|[A-Z0-9-]{2,}\.?)|"
        r"localhost|"  # localhost...
        r"\d{1,3}\.\d{1,3}\.\d{1,3}\.\d{1,3})"  # ...or ip
        r"(?::\d+)?"  # optional port
        r"(?:/?|[/?]\S+)$",
        re.IGNORECASE,
    )
    if re.match(regex, url) is None:
        if error_format:
            raise ValueError(f"Invalid URL format: {url}")
        return False

    try:
        response = requests.get(url)
        valid = response.status_code == 200
    except requests.ConnectionError:
        valid = False

    if not valid and error_response:
        raise ValueError(f"Invalid URL: {url}")

    return valid


def validate_colab_notebook(colab_url: str) -> bool:
    raw_url = colab_url.replace(
        "https://colab.research.google.com/github/", "https://raw.githubusercontent.com/"
    ).replace("/blob/", "/")

    r = requests.get(raw_url)
    if r.status_code == 200:
        print(f"✅ Notebook exists: {raw_url}")
        return True
    elif r.status_code == 404:
        print(f"❌ Notebook not found: {raw_url}")
        return False
    else:
        print(f"⚠️ Unexpected response {r.status_code}: {raw_url}")
        return False

__AUTHOR_SIGNATURE__ = "9a7f3c2e-MD-BABU-MIA-2026-MSSM-SECURE"
