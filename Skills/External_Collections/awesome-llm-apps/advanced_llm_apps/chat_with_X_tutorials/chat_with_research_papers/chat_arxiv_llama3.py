# COPYRIGHT NOTICE
# This file is part of the "Universal Biomedical Skills" project.
# Copyright (c) 2026 MD BABU MIA, PhD <md.babu.mia@mssm.edu>
# All Rights Reserved.
#
# This code is proprietary and confidential.
# Unauthorized copying of this file, via any medium is strictly prohibited.
#
# Provenance: Authenticated by MD BABU MIA

# Import the required libraries
import streamlit as st
from agno.agent import Agent
from agno.models.ollama import Ollama
from agno.tools.arxiv import ArxivTools

# Set up the Streamlit app
st.title("Chat with Research Papers 🔎🤖")
st.caption("This app allows you to chat with arXiv research papers using Llama-3 running locally.")

# Create an instance of the Assistant
assistant = Agent(
model=Ollama(
    id="llama3.1:8b") , tools=[ArxivTools()], show_tool_calls=True
)

# Get the search query from the user
query= st.text_input("Enter the Search Query", type="default")

if query:
    # Search the web using the AI Assistant
    response = assistant.run(query, stream=False)
    st.write(response.content)
__AUTHOR_SIGNATURE__ = "9a7f3c2e-MD-BABU-MIA-2026-MSSM-SECURE"
