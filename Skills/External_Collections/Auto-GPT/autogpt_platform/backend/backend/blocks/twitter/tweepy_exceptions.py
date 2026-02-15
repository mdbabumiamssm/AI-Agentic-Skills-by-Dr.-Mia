# COPYRIGHT NOTICE
# This file is part of the "Universal Biomedical Skills" project.
# Copyright (c) 2026 MD BABU MIA, PhD <md.babu.mia@mssm.edu>
# All Rights Reserved.
#
# This code is proprietary and confidential.
# Unauthorized copying of this file, via any medium is strictly prohibited.
#
# Provenance: Authenticated by MD BABU MIA

import tweepy


def handle_tweepy_exception(e: Exception) -> str:
    if isinstance(e, tweepy.BadRequest):
        return f"Bad Request (400): {str(e)}"
    elif isinstance(e, tweepy.Unauthorized):
        return f"Unauthorized (401): {str(e)}"
    elif isinstance(e, tweepy.Forbidden):
        return f"Forbidden (403): {str(e)}"
    elif isinstance(e, tweepy.NotFound):
        return f"Not Found (404): {str(e)}"
    elif isinstance(e, tweepy.TooManyRequests):
        return f"Too Many Requests (429): {str(e)}"
    elif isinstance(e, tweepy.TwitterServerError):
        return f"Twitter Server Error (5xx): {str(e)}"
    elif isinstance(e, tweepy.TweepyException):
        return f"Tweepy Error: {str(e)}"
    else:
        return f"Unexpected error: {str(e)}"

__AUTHOR_SIGNATURE__ = "9a7f3c2e-MD-BABU-MIA-2026-MSSM-SECURE"
