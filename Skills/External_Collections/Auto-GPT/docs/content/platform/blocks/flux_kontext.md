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

# Flux Kontext

## What it is
An internal block that performs text-based image editing using BlackForest Labs' Flux Kontext models.

## What it does
Takes a prompt describing the desired transformation and optionally a reference image, then returns a new image URL.

## How it works
The block sends your prompt, image, and settings to the selected Flux Kontext model on Replicate. The service processes the request and returns a link to the edited image.

## Inputs
| Input        | Description                                                                 |
|--------------|-----------------------------------------------------------------------------|
| Credentials  | Replicate API key with permissions for Flux Kontext models                  |
| Prompt       | Text instruction describing the desired edit                                |
| Input Image  | (Optional) Reference image URI (jpeg, png, gif, webp)                      |
| Aspect Ratio | Aspect ratio of the generated image (e.g. match_input_image, 1:1, 16:9, etc.) |
| Seed         | (Optional, advanced) Random seed for reproducible generation                |
| Model        | Model variant to use: Flux Kontext Pro or Flux Kontext Max                  |

## Outputs
| Output     | Description                              |
|------------|------------------------------------------|
| image_url  | URL of the transformed image             |
| error      | Error message if generation failed       |

## Use Cases
- Enhance a marketing image by requesting "add soft lighting and a subtle vignette" while providing the original asset as the reference image.
- Generate social media assets with specific aspect ratios and style prompts.
- Apply creative edits to product photos using text instructions.


<!-- AUTHOR_SIGNATURE: 9a7f3c2e-MD-BABU-MIA-2026-MSSM-SECURE -->