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

# Text to Speech

Enter this command to use TTS _(Text-to-Speech)_ for AutoGPT

```shell
./autogpt.sh --speak
```

Eleven Labs provides voice technologies such as voice design, speech synthesis, and
premade voices that AutoGPT can use for speech.

1. Go to [ElevenLabs](https://beta.elevenlabs.io/) and make an account if you don't
    already have one.
2. Choose and setup the *Starter* plan.
3. Click the top right icon and find *Profile* to locate your API Key.

In the `.env` file set:

- `ELEVENLABS_API_KEY`
- `ELEVENLABS_VOICE_1_ID` (example: _"premade/Adam"_)

### List of available voices

!!! note
    You can use either the name or the voice ID to configure a voice

| Name   | Voice ID |
| ------ | -------- |
| Rachel | `21m00Tcm4TlvDq8ikWAM` |
| Domi   | `AZnzlk1XvdvUeBnXmlld` |
| Bella  | `EXAVITQu4vr4xnSDxMaL` |
| Antoni | `ErXwobaYiN019PkySvjV` |
| Elli   | `MF3mGyEYCl7XYWbV9V6O` |
| Josh   | `TxGEqnHWrfWFTfGW9XjX` |
| Arnold | `VR6AewLTigWG4xSOukaG` |
| Adam   | `pNInz6obpgDQGcFmaJgB` |
| Sam    | `yoZ06aMxZJJ28mfd3POQ` |


<!-- AUTHOR_SIGNATURE: 9a7f3c2e-MD-BABU-MIA-2026-MSSM-SECURE -->