# COPYRIGHT NOTICE
# This file is part of the "Universal Biomedical Skills" project.
# Copyright (c) 2026 MD BABU MIA, PhD <md.babu.mia@mssm.edu>
# All Rights Reserved.
#
# This code is proprietary and confidential.
# Unauthorized copying of this file, via any medium is strictly prohibited.
#
# Provenance: Authenticated by MD BABU MIA

from fastapi import FastAPI
from fastapi.middleware.cors import CORSMiddleware
from fastapi.encoders import jsonable_encoder
from fastapi.responses import JSONResponse
from pydantic import BaseModel
from agents.coordinator_agent import coordinator_agent
from agno.agent import RunResponse
from dotenv import load_dotenv

# Load environment variables
load_dotenv()

# Initialize FastAPI app
app = FastAPI()

# Configure CORS to allow requests from your frontend
app.add_middleware(
    CORSMiddleware,
    allow_origins=["*"],  # To be replaced with the frontend's origin in production
    allow_credentials=True,
    allow_methods=["*"],
    allow_headers=["*"],
)

# Define the request body model
class AnalysisRequest(BaseModel):
    video_url: str

# Define the entry point
@app.get("/")
async def root():
    return {"message": "Welcome to the video analysis API!"}

# Define the analysis endpoint
@app.post("/analyze")
async def analyze(request: AnalysisRequest):
    video_url = request.video_url
    prompt = f"Analyze the following video: {video_url}"
    response: RunResponse = coordinator_agent.run(prompt)

    # Assuming response.content is a Pydantic model or a dictionary
    json_compatible_response = jsonable_encoder(response.content)
    return JSONResponse(content=json_compatible_response)
__AUTHOR_SIGNATURE__ = "9a7f3c2e-MD-BABU-MIA-2026-MSSM-SECURE"
