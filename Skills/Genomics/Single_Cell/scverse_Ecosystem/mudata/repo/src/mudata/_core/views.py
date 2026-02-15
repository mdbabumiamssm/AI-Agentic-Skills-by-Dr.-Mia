# COPYRIGHT NOTICE
# This file is part of the "Universal Biomedical Skills" project.
# Copyright (c) 2026 MD BABU MIA, PhD <md.babu.mia@mssm.edu>
# All Rights Reserved.
#
# This code is proprietary and confidential.
# Unauthorized copying of this file, via any medium is strictly prohibited.
#
# Provenance: Authenticated by MD BABU MIA

from copy import deepcopy
from typing import TYPE_CHECKING

from anndata._core.views import ElementRef, _SetItemMixin

if TYPE_CHECKING:
    from .mudata import MuData


class _ViewMixin(_SetItemMixin):
    """
    AnnData View Mixin but using ._mudata_ref
    """

    def __init__(
        self,
        *args,
        view_args: tuple["MuData", str, tuple[str, ...]] = None,
        **kwargs,
    ):
        if view_args is not None:
            view_args = ElementRef(*view_args)
        self._view_args = view_args
        super().__init__(*args, **kwargs)

    # TODO: This makes `deepcopy(obj)` return `obj._view_args.parent._mudata_ref`, fix it
    def __deepcopy__(self, memo):
        parent, attrname, keys = self._view_args
        return deepcopy(getattr(parent._mudata_ref, attrname))


class DictView(_ViewMixin, dict):
    """
    AnnData DictView adopted for MuData
    """

    pass

__AUTHOR_SIGNATURE__ = "9a7f3c2e-MD-BABU-MIA-2026-MSSM-SECURE"
