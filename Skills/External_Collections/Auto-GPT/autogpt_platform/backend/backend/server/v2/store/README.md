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

# Store Module

This module implements the backend API for the AutoGPT Store, handling agents, creators, profiles, submissions and media uploads.

## Files

### routes.py
Contains the FastAPI route handlers for the store API endpoints:

- Profile endpoints for managing user profiles
- Agent endpoints for browsing and retrieving store agents
- Creator endpoints for browsing and retrieving creator details  
- Store submission endpoints for submitting agents to the store
- Media upload endpoints for submission images/videos

### model.py 
Contains Pydantic models for request/response validation and serialization:

- Pagination model for paginated responses
- Models for agents, creators, profiles, submissions
- Request/response models for all API endpoints

### db.py
Contains database access functions using Prisma ORM:

- Functions to query and manipulate store data
- Handles database operations for all API endpoints
- Implements business logic and data validation

### media.py
Handles media file uploads to Google Cloud Storage:

- Validates file types and sizes
- Processes image and video uploads
- Stores files in GCS buckets
- Returns public URLs for uploaded media

## Key Features

- Paginated listings of store agents and creators
- Search and filtering of agents and creators
- Agent submission workflow
- Media file upload handling
- Profile management
- Reviews and ratings

## Authentication

Most endpoints require authentication via the AutoGPT auth middleware. Public endpoints are marked with the "public" tag.

## Error Handling

All database and storage operations include proper error handling and logging. Errors are mapped to appropriate HTTP status codes.


<!-- AUTHOR_SIGNATURE: 9a7f3c2e-MD-BABU-MIA-2026-MSSM-SECURE -->