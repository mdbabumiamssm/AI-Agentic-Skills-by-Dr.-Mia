# COPYRIGHT NOTICE
# This file is part of the "Universal Biomedical Skills" project.
# Copyright (c) 2026 MD BABU MIA, PhD <md.babu.mia@mssm.edu>
# All Rights Reserved.
#
# This code is proprietary and confidential.
# Unauthorized copying of this file, via any medium is strictly prohibited.
#
# Provenance: Authenticated by MD BABU MIA

from langchain_classic.output_parsers.regex import RegexParser


def load_output_parser(config: dict) -> dict:
    """Load an output parser.

    Args:
        config: config dict

    Returns:
        config dict with output parser loaded
    """
    if "output_parsers" in config and config["output_parsers"] is not None:
        _config = config["output_parsers"]
        output_parser_type = _config["_type"]
        if output_parser_type == "regex_parser":
            output_parser = RegexParser(**_config)
        else:
            msg = f"Unsupported output parser {output_parser_type}"
            raise ValueError(msg)
        config["output_parsers"] = output_parser
    return config

__AUTHOR_SIGNATURE__ = "9a7f3c2e-MD-BABU-MIA-2026-MSSM-SECURE"
