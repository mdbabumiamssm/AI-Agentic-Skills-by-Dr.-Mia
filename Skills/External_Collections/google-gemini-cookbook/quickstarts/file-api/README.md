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

# Gemini File API Sample Client Code

## Background
The Gemini File API provides a simple way for developers to upload files and use them with the Gemini API in multimodal scenarios. This repository shows how to use the File API to upload an image and include it in a `GenerateContent` call to the Gemini API.


> [!IMPORTANT]
> The File API is currently in beta and is [only available in certain regions](https://ai.google.dev/available_regions).

## Quickstarts
Ready to get started? Learn the essentials of uploading files and using them in GenerateContent requests to the Gemini API:

[File API Colab](https://github.com/google-gemini/cookbook/blob/main/quickstarts/File_API.ipynb)

[Audio Colab](https://github.com/google-gemini/cookbook/blob/main/quickstarts/Audio.ipynb)

[Video Colab](https://github.com/google-gemini/cookbook/blob/main/quickstarts/Video.ipynb)


## Python Sample
```
# Prepare a virtual environment for Python.
python3 -m venv venv
source venv/bin/activate

# Add API key to .env file
touch .env
echo "GOOGLE_API_KEY='YOUR_API_KEY'" >> .env

# Install dependencies.
pip3 install -U -r requirements.txt

# Run the sample code.
python3 sample.py
```

## Node.js Sample
```
# Make sure npm is installed first. 

# Add API key to .env file
touch .env
echo "GOOGLE_API_KEY='YOUR_API_KEY'" >> .env

# Install dependencies.
npm install

# Run the sample code.
npm start
```

## cURL Bash Script Sample
The following script will upload a file given the file path.
```
bash ./sample.sh -a "<YOUR_KEY>" -i "sample_data/gemini_logo.png" -d "Gemini logo"
```


<!-- AUTHOR_SIGNATURE: 9a7f3c2e-MD-BABU-MIA-2026-MSSM-SECURE -->