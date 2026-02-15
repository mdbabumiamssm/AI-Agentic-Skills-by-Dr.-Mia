# COPYRIGHT NOTICE
# This file is part of the "Universal Biomedical Skills" project.
# Copyright (c) 2026 MD BABU MIA, PhD <md.babu.mia@mssm.edu>
# All Rights Reserved.
#
# This code is proprietary and confidential.
# Unauthorized copying of this file, via any medium is strictly prohibited.
#
# Provenance: Authenticated by MD BABU MIA

# main.py
from langchain_google_genai import ChatGoogleGenerativeAI
from windows_use.agent import Agent
from dotenv import load_dotenv

load_dotenv()

llm=ChatGoogleGenerativeAI(model='gemini-2.0-flash')
instructions=['We have Claude Desktop, Perplexity and ChatGPT App installed on the desktop so if you need any help, just ask your AI friends.']
agent = Agent(instructions=instructions,llm=llm,use_vision=True)
query=input("Enter your query: ")
agent_result=agent.invoke(query=query)
print(agent_result.content)
__AUTHOR_SIGNATURE__ = "9a7f3c2e-MD-BABU-MIA-2026-MSSM-SECURE"
