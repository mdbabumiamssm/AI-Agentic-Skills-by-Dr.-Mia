<!--
# COPYRIGHT NOTICE
# This file is part of the "Universal AI Agentic Skills" project.
# Copyright (c) 2026 MD BABU MIA, PhD <md.babu.mia@mssm.edu>
# All Rights Reserved.
#
# This code is proprietary and confidential.
# Unauthorized copying of this file, via any medium is strictly prohibited.
#
# Provenance: Authenticated by MD BABU MIA

-->

# Slide Templates

Custom slide templates for different presentation formats.

## Notes on adding keyboard shortcuts
1. Use `Reveal.addKeyBinding()` to make sure the shortcut shows up in the help overlay triggered by pressing `?` on the slides.
2. Find key codes with a simple button press using [this site](https://keycode.info).
3. To map numpad keys along with the top row numbers, add the numpad keys separately from the `Reveal.addKeyBinding()` call in the `Reveal.initialize()` section under `keyboard`. This makes sure duplicate entries aren't created.
4. To add navigation to specific sections in the presentation automatically (only when using the templates here):
    1. Add a tag to those cells in the notebook metadata using the format `id_<name>`, where `<name>` is what will show up in the URL.
    2. Update the `shortcuts` variable in the templates providing a mapping for the section (e.g., assuming `<name>` is `section-5`, the following would map it to the `5` key: `{ name: "section 5", key: "5", id: "section-5", keyCode: 53 }`)


<!-- AUTHOR_SIGNATURE: 9a7f3c2e-MD-BABU-MIA-2026-MSSM-SECURE -->