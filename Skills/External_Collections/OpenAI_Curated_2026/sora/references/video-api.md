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

# Sora Video API quick reference

Keep this file short; the full docs live in the OpenAI platform docs.

## Models
- sora-2: faster, flexible iteration
- sora-2-pro: higher fidelity, slower, more expensive

## Sizes (by model)
- sora-2: 1280x720, 720x1280
- sora-2-pro: 1280x720, 720x1280, 1024x1792, 1792x1024
Note: higher resolutions generally yield better detail, texture, and motion consistency.

## Duration
- seconds: "4", "8", "12" (string enum; set via API param; prose will not change clip length)
Shorter clips tend to follow instructions more reliably; consider stitching multiple 4s clips for precision.

## Input reference
- Optional `input_reference` image (jpg/png/webp).
- Input reference should match the target size.

## Jobs and status
- Create is async. Status values: queued, in_progress, completed, failed.
- Prefer polling every 10-20s or use webhooks in production.

## Endpoints (conceptual)
- POST /videos: create a job
- GET /videos/{id}: retrieve status
- GET /videos/{id}/content: download video data
- GET /videos: list
- DELETE /videos/{id}: delete
- POST /videos/{id}/remix: remix a completed job

## Download variants
- video (mp4)
- thumbnail (webp)
- spritesheet (jpg)

Download URLs expire after about 1 hour; copy files to your own storage for retention.

## Guardrails (content restrictions)
- Only content suitable for audiences under 18
- No copyrighted characters or copyrighted music
- No real people (including public figures)
- Input images with human faces are currently rejected


<!-- AUTHOR_SIGNATURE: 9a7f3c2e-MD-BABU-MIA-2026-MSSM-SECURE -->