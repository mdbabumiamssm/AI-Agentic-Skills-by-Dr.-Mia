# COPYRIGHT NOTICE
# This file is part of the "Universal Biomedical Skills" project.
# Copyright (c) 2026 MD BABU MIA, PhD <md.babu.mia@mssm.edu>
# All Rights Reserved.
#
# This code is proprietary and confidential.
# Unauthorized copying of this file, via any medium is strictly prohibited.
#
# Provenance: Authenticated by MD BABU MIA

def calculator(expression):
    """Evaluates a mathematical expression."""
    try:
        # unsafe eval for demo purposes only; use a safe parser in prod
        return str(eval(expression))
    except Exception as e:
        return f"Error: {e}"

def lookup_weather(location):
    """Mock weather lookup."""
    if "london" in location.lower():
        return "Rainy, 15C"
    elif "new york" in location.lower():
        return "Sunny, 22C"
    else:
        return "Weather data unavailable for this location."

registry = {
    "Calculator": calculator,
    "Weather": lookup_weather
}

__AUTHOR_SIGNATURE__ = "9a7f3c2e-MD-BABU-MIA-2026-MSSM-SECURE"
