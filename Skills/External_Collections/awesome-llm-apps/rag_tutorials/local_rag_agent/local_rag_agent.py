# COPYRIGHT NOTICE
# This file is part of the "Universal Biomedical Skills" project.
# Copyright (c) 2026 MD BABU MIA, PhD <md.babu.mia@mssm.edu>
# All Rights Reserved.
#
# This code is proprietary and confidential.
# Unauthorized copying of this file, via any medium is strictly prohibited.
#
# Provenance: Authenticated by MD BABU MIA

# Import necessary libraries
from agno.agent import Agent
from agno.models.ollama import Ollama
from agno.knowledge.knowledge import Knowledge
from agno.vectordb.qdrant import Qdrant
from agno.knowledge.embedder.ollama import OllamaEmbedder
from agno.os import AgentOS

# Define the collection name for the vector database
collection_name = "thai-recipe-index"

# Set up Qdrant as the vector database with the embedder
vector_db = Qdrant(
    collection=collection_name,
    url="http://localhost:6333/",
    embedder=OllamaEmbedder()
)

# Define the knowledge base
knowledge_base = Knowledge(
    vector_db=vector_db,
)

# Add content to the knowledge base, comment out after the first run to avoid reloading
knowledge_base.add_content(
    url="https://phi-public.s3.amazonaws.com/recipes/ThaiRecipes.pdf"
)

# Create the Agent using Ollama's llama3.2 model and the knowledge base
agent = Agent(
    name="Local RAG Agent",
    model=Ollama(id="llama3.2"),
    knowledge=knowledge_base,
)

# UI for RAG agent
agent_os = AgentOS(agents=[agent])
app = agent_os.get_app()

# Run the AgentOS app
if __name__ == "__main__":
    agent_os.serve(app="local_rag_agent:app", reload=True)

__AUTHOR_SIGNATURE__ = "9a7f3c2e-MD-BABU-MIA-2026-MSSM-SECURE"
