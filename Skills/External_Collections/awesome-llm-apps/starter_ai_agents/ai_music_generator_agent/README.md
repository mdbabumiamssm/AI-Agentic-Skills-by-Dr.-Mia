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

## ModelsLab Music Generator

This is a Streamlit-based application that allows users to generate music using the ModelsLab API and OpenAI's GPT-4 model. Users can input a prompt describing the type of music they want to generate, and the application will generate a music track in MP3 format based on the given prompt.

## Features

- **Generate Music**: Enter a detailed prompt for music generation (genre, instruments, mood, etc.), and the app will generate a music track.
- **MP3 Output**: The generated music will be in MP3 format, available for listening or download.
- **User-Friendly Interface**: Simple and clean Streamlit UI for ease of use.
- **API Key Integration**: Requires both OpenAI and ModelsLab API keys to function. API keys are entered in the sidebar for authentication.

## Setup

### Requirements 

1. **API Keys**:
   - **OpenAI API Key**: Sign up at [OpenAI](https://platform.openai.com/api-keys) to obtain your API key.
   - **ModelsLab API Key**: Sign up at [ModelsLab](https://modelslab.com/dashboard/api-keys) to get your API key.

2. **Python 3.8+**: Ensure you have Python 3.8 or higher installed.

### Installation
1. Clone this repository:
   ```bash
   git clone https://github.com/Shubhamsaboo/awesome-llm-apps
   cd ai_agent_tutorials/ai_models_lab_music_generator_agent
   ```

2. Install the required Python packages:
   ```bash
   pip install -r requirements.txt
   ```
### Running the App

1. Start the Streamlit app:
   ```bash
   streamlit run models_lab_music_generator_agent.py
   ```

2. In the app interface:
   - Enter a music generation prompt
   - Click "Generate Music"
   - Play the music & Download it.

<!-- AUTHOR_SIGNATURE: 9a7f3c2e-MD-BABU-MIA-2026-MSSM-SECURE -->