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

# Find available voices for D-ID

1. **ElevenLabs**
   - Select any voice from the voice list: https://api.elevenlabs.io/v1/voices
   - Copy the voice_id
   - Use it as a string in the voice_id field in the CreateTalkingAvatarClip Block

2. **Microsoft Azure Voices**
    - Select any voice from the voice gallery: https://speech.microsoft.com/portal/voicegallery
    - Click on the "Sample code" tab on the right
    - Copy the voice name, for example: config.SpeechSynthesisVoiceName ="en-GB-AbbiNeural"
    - Use this string en-GB-AbbiNeural in the voice_id field in the CreateTalkingAvatarClip Block

3. **Amazon Polly Voices**
    - Select any voice from the voice list: https://docs.aws.amazon.com/polly/latest/dg/available-voices.html
    - Copy the voice name / ID
    - Use it as string in the voice_id field in the CreateTalkingAvatarClip Block

<!-- AUTHOR_SIGNATURE: 9a7f3c2e-MD-BABU-MIA-2026-MSSM-SECURE -->