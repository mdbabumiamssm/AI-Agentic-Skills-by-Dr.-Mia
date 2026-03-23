# COPYRIGHT NOTICE
# This file is part of the "Universal Biomedical Skills" project.
# Copyright (c) 2026 MD BABU MIA, PhD <md.babu.mia@mssm.edu>
# All Rights Reserved.
#
# This code is proprietary and confidential.
# Unauthorized copying of this file, via any medium is strictly prohibited.
#
# Provenance: Authenticated by MD BABU MIA

"""Utility function for the workshop."""

import datetime as dt

import pandas as pd


def mpl_svg_config(hashsalt):
    """Help configure the SVG backend for Matplotlib and make it reproducible."""
    from matplotlib import rc
    rc('svg', hashsalt=hashsalt)

    return {
        'metadata': {
            'Date': f'(c) 2021-{dt.date.today().year} Stefanie Molin'
        }
    }

def highlight_long_format(df, colors):
    """Highlight long format columns in the data."""

    def get_style(x):
        if color := colors.get(x):
            return f'background-color: {color}'
        return None

    def highlight_column(x):
        return [get_style(x.name)] * x.shape[0]

    return df.style\
        .apply_index(lambda x: x.apply(get_style), axis=1)\
        .apply(highlight_column, axis=0)

def highlight_wide_format(df, colors):
    """Highlight wide format columns in the melted data."""

    def highlight_melt(x):
        color = colors[x.year]
        return [f"background-color: {color};"] * 2

    return df.melt(id_vars='date', var_name='year', value_name='travelers')\
        .assign(
            date=lambda x: pd.to_datetime(x.year + x.date.dt.strftime('-%m-%d'))
        )\
        .sort_values('date', ascending=False)\
        .style.apply(highlight_melt, subset=['year', 'travelers'], axis=1)

__AUTHOR_SIGNATURE__ = "9a7f3c2e-MD-BABU-MIA-2026-MSSM-SECURE"
