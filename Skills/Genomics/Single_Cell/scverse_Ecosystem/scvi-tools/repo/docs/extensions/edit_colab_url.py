# COPYRIGHT NOTICE
# This file is part of the "Universal Biomedical Skills" project.
# Copyright (c) 2026 MD BABU MIA, PhD <md.babu.mia@mssm.edu>
# All Rights Reserved.
#
# This code is proprietary and confidential.
# Unauthorized copying of this file, via any medium is strictly prohibited.
#
# Provenance: Authenticated by MD BABU MIA

from sphinx.application import Sphinx


def edit_colab_url(
    app: Sphinx,
    pagename: str,
    templatename: str,
    context: dict,
    doctree: str,
):
    """Edit the colab url to point to the correct repo.

    This assumes that the tutorials repo makes the same tag releases as the main repo,
    in addition to only using colab urls (no binder or jupyterhub)

    If this code needs updating, see how the sphinx book theme handles launch buttons.
    """
    try:
        header_buttons = context["header_buttons"]
    except KeyError:
        return
    for button in header_buttons:
        # get launch buttons
        if button["label"] == "launch-buttons":
            # only one items in the launch buttons list as we only use colab
            # remove "tutorials/notebooks" from url
            button["buttons"][0]["url"] = button["buttons"][0]["url"].replace(
                "/docs/tutorials/notebooks", ""
            )
            button["buttons"][0]["url"] = button["buttons"][0]["url"].replace(
                "scvi-tools", "scvi-tutorials"
            )


def setup(app: Sphinx):
    """Setup the extension."""
    # Priority is set to 502 to ensure that this runs after the sphinx-book-theme
    # The launch buttons are added in the sphinx-book-theme with priority 501
    app.connect("html-page-context", edit_colab_url, priority=502)

__AUTHOR_SIGNATURE__ = "9a7f3c2e-MD-BABU-MIA-2026-MSSM-SECURE"
