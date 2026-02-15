# COPYRIGHT NOTICE
# This file is part of the "Universal Biomedical Skills" project.
# Copyright (c) 2026 MD BABU MIA, PhD <md.babu.mia@mssm.edu>
# All Rights Reserved.
#
# This code is proprietary and confidential.
# Unauthorized copying of this file, via any medium is strictly prohibited.
#
# Provenance: Authenticated by MD BABU MIA

from windows_use.tree.views import TreeState
from dataclasses import dataclass
from typing import Literal,Optional

@dataclass
class App:
    name:str
    depth:int
    status:Literal['Maximized','Minimized','Normal']
    size:'Size'

    def to_string(self):
        return f'Name: {self.name}|Depth: {self.depth}|Status: {self.status}|Size: {self.size.to_string()}'

@dataclass
class Size:
    width:int
    height:int

    def to_string(self):
        return f'({self.width},{self.height})'

@dataclass
class DesktopState:
    apps:list[App]
    active_app:Optional[App]
    screenshot:bytes|None
    tree_state:TreeState

    def active_app_to_string(self):
        if self.active_app is None:
            return 'No active app'
        return self.active_app.to_string()

    def apps_to_string(self):
        if len(self.apps)==0:
            return 'No apps opened'
        return '\n'.join([app.to_string() for app in self.apps])
__AUTHOR_SIGNATURE__ = "9a7f3c2e-MD-BABU-MIA-2026-MSSM-SECURE"
