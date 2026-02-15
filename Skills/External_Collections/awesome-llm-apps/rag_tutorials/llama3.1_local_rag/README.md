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

## 💻 Local Lllama-3.1 with RAG
Streamlit app that allows you to chat with any webpage using local Llama-3.1 and Retrieval Augmented Generation (RAG). This app runs entirely on your computer, making it 100% free and without the need for an internet connection.


### Features
- Input a webpage URL
- Ask questions about the content of the webpage
- Get accurate answers using RAG and the Llama-3.1 model running locally on your computer

### How to get Started?

1. Clone the GitHub repository

```bash
git clone https://github.com/Shubhamsaboo/awesome-llm-apps.git
cd awesome-llm-apps/rag_tutorials/llama3.1_local_rag
```
2. Install the required dependencies:

```bash
pip install -r requirements.txt
```
3. Run the Streamlit App
```bash
streamlit run llama3.1_local_rag.py
```

### How it Works?

- The app loads the webpage data using WebBaseLoader and splits it into chunks using RecursiveCharacterTextSplitter.
- It creates Ollama embeddings and a vector store using Chroma.
- The app sets up a RAG (Retrieval-Augmented Generation) chain, which retrieves relevant documents based on the user's question.
- The Llama-3.1 model is called to generate an answer using the retrieved context.
- The app displays the answer to the user's question.



<!-- AUTHOR_SIGNATURE: 9a7f3c2e-MD-BABU-MIA-2026-MSSM-SECURE -->