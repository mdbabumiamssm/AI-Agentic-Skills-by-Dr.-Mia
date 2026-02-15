<!--
# COPYRIGHT NOTICE
# This file is part of the "Universal Biomedical Skills" project.
# Copyright (c) 2026 MD BABU MIA, PhD <md.babu.mia@mssm.edu>
# All Rights Reserved.
#
# This code is proprietary and confidential.
# Unauthorized copying of this file, via any medium is strictly prohibited.
#
# Provenance: Authenticated by MD BABU MIA

-->

# gpt-4o-transcribe-diarize quick reference

- Input formats: mp3, mp4, mpeg, mpga, m4a, wav, webm.
- Max file size: 25 MB per request.
- response_format options: text, json, diarized_json.
- For audio longer than ~30 seconds, pass chunking_strategy (use "auto" to split into chunks).
- Known speakers: up to 4 references via extra_body known_speaker_names + known_speaker_references (data URLs).
- Prompting is not supported for gpt-4o-transcribe-diarize.


<!-- AUTHOR_SIGNATURE: 9a7f3c2e-MD-BABU-MIA-2026-MSSM-SECURE -->