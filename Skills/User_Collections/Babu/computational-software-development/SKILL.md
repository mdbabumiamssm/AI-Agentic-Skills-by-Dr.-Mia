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

---
name: computational-software-development
description: "Full-stack computational software development for biomedical and life science applications. Use when building AI-powered research platforms, clinical decision support systems, multi-LLM ensemble architectures, RAG pipelines, single-cell analysis tools, or biomedical web applications. Includes production-ready patterns for PyTorch deep learning, Flask web apps, citation verification systems, and scientific software distribution."
license: Proprietary
---

# Computational Software Development for Life Sciences

## Core Technologies & Architecture Patterns

### Technology Stack Overview

```python
# AI/ML Stack
import torch
import torch.nn as nn
from transformers import AutoModel, AutoTokenizer
import anthropic  # Claude API
import openai     # GPT API
import google.generativeai as genai  # Gemini API

# Scientific Computing
import scanpy as sc
import anndata as ad
import scvi
import pandas as pd
import numpy as np
from scipy import stats

# Web Framework
from flask import Flask, render_template, jsonify, request
import asyncio
import aiohttp

# Data Management
import faiss  # Vector database
from sentence_transformers import SentenceTransformer
import redis
from diskcache import Cache
```

---

## 1. Multi-LLM Ensemble Architecture

### Weighted Consensus System (BioMedAI Pattern)

```python
from dataclasses import dataclass
from typing import List, Dict, Optional
import asyncio
from abc import ABC, abstractmethod

@dataclass
class LLMResponse:
    content: str
    model: str
    confidence: float
    citations: List[Dict]
    latency_ms: float

class AIProvider(ABC):
    @abstractmethod
    async def query(self, prompt: str, context: str) -> LLMResponse:
        pass

class ClaudeProvider(AIProvider):
    def __init__(self, api_key: str, model: str = "claude-opus-4-5-20251101"):
        self.client = anthropic.Anthropic(api_key=api_key)
        self.model = model
        self.weight = 1.2  # Primary synthesizer weight

    async def query(self, prompt: str, context: str) -> LLMResponse:
        start = time.time()
        response = await asyncio.to_thread(
            self.client.messages.create,
            model=self.model,
            max_tokens=4096,
            messages=[{"role": "user", "content": f"{context}\n\n{prompt}"}]
        )
        return LLMResponse(
            content=response.content[0].text,
            model=self.model,
            confidence=0.0,  # Computed post-hoc
            citations=[],
            latency_ms=(time.time() - start) * 1000
        )

class MultiLLMEnsemble:
    """Parallel multi-LLM execution with weighted consensus."""

    def __init__(self, providers: List[AIProvider]):
        self.providers = providers
        self.weights = {p.__class__.__name__: getattr(p, 'weight', 1.0) for p in providers}

    async def query_parallel(self, prompt: str, context: str) -> Dict:
        """Execute all LLMs in parallel for <5s response time."""
        tasks = [provider.query(prompt, context) for provider in self.providers]
        responses = await asyncio.gather(*tasks, return_exceptions=True)

        valid_responses = [r for r in responses if isinstance(r, LLMResponse)]
        return self._synthesize_consensus(valid_responses)

    def _synthesize_consensus(self, responses: List[LLMResponse]) -> Dict:
        """Weighted synthesis of multi-model responses."""
        # Extract key claims from each response
        all_claims = []
        for r in responses:
            claims = self._extract_claims(r.content)
            weight = self.weights.get(r.model.split('-')[0].title() + 'Provider', 1.0)
            all_claims.extend([(c, weight) for c in claims])

        # Rank claims by weighted frequency
        claim_scores = {}
        for claim, weight in all_claims:
            claim_scores[claim] = claim_scores.get(claim, 0) + weight

        return {
            'consensus': sorted(claim_scores.items(), key=lambda x: -x[1]),
            'individual_responses': responses,
            'agreement_score': self._compute_agreement(responses)
        }
```

---

## 2. Hybrid RAG v3.0 Architecture

### Dense + Sparse Retrieval with Reciprocal Rank Fusion

```python
from sentence_transformers import SentenceTransformer
import faiss
from rank_bm25 import BM25Okapi
from typing import List, Tuple
import numpy as np

class HybridRAG:
    """Production-ready RAG with dense and sparse retrieval fusion."""

    def __init__(self, embedding_model: str = "intfloat/e5-large-v2"):
        self.embedder = SentenceTransformer(embedding_model)
        self.dense_index = None
        self.sparse_index = None
        self.documents = []

    def build_index(self, documents: List[Dict]):
        """Build dual dense + sparse indices."""
        self.documents = documents
        texts = [d['text'] for d in documents]

        # Dense index (FAISS)
        embeddings = self.embedder.encode(texts, normalize_embeddings=True)
        dimension = embeddings.shape[1]
        self.dense_index = faiss.IndexFlatIP(dimension)  # Inner product for cosine sim
        self.dense_index.add(embeddings.astype('float32'))

        # Sparse index (BM25)
        tokenized = [t.lower().split() for t in texts]
        self.sparse_index = BM25Okapi(tokenized)

    def retrieve(self, query: str, k: int = 10, alpha: float = 0.5) -> List[Dict]:
        """Reciprocal Rank Fusion of dense + sparse results."""
        # Dense retrieval
        query_emb = self.embedder.encode([query], normalize_embeddings=True)
        dense_scores, dense_ids = self.dense_index.search(query_emb.astype('float32'), k * 2)

        # Sparse retrieval
        sparse_scores = self.sparse_index.get_scores(query.lower().split())
        sparse_ids = np.argsort(sparse_scores)[-k * 2:][::-1]

        # Reciprocal Rank Fusion (RRF)
        rrf_k = 60  # Standard RRF parameter
        doc_scores = {}

        for rank, doc_id in enumerate(dense_ids[0]):
            doc_scores[doc_id] = doc_scores.get(doc_id, 0) + alpha / (rrf_k + rank + 1)

        for rank, doc_id in enumerate(sparse_ids):
            doc_scores[doc_id] = doc_scores.get(doc_id, 0) + (1 - alpha) / (rrf_k + rank + 1)

        # Return top-k fused results
        ranked = sorted(doc_scores.items(), key=lambda x: -x[1])[:k]
        return [{'document': self.documents[did], 'score': score} for did, score in ranked]
```

---

## 3. Citation Validation Pipeline

### Real-Time PubMed Citation Verification

```python
import re
import aiohttp
from typing import Optional, Dict, List

class CitationValidator:
    """3-layer citation validation: Format → Existence → Relevance."""

    PUBMED_API = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"

    async def validate_citation(self, citation: str, claim_context: str) -> Dict:
        """Full validation pipeline for a single citation."""
        result = {
            'original': citation,
            'format_valid': False,
            'exists': False,
            'relevant': False,
            'pmid': None,
            'title': None,
            'abstract': None
        }

        # Layer 1: Format validation
        pmid = self._extract_pmid(citation)
        if not pmid:
            pmid = await self._search_pubmed(citation)

        if not pmid:
            return result

        result['format_valid'] = True
        result['pmid'] = pmid

        # Layer 2: Existence check
        metadata = await self._fetch_pubmed_metadata(pmid)
        if not metadata:
            return result

        result['exists'] = True
        result['title'] = metadata.get('title')
        result['abstract'] = metadata.get('abstract')

        # Layer 3: Relevance scoring
        relevance = self._compute_relevance(claim_context, metadata)
        result['relevant'] = relevance > 0.3
        result['relevance_score'] = relevance

        return result

    async def _fetch_pubmed_metadata(self, pmid: str) -> Optional[Dict]:
        """Fetch article metadata from PubMed."""
        url = f"{self.PUBMED_API}/efetch.fcgi"
        params = {'db': 'pubmed', 'id': pmid, 'retmode': 'xml'}

        async with aiohttp.ClientSession() as session:
            async with session.get(url, params=params) as resp:
                if resp.status == 200:
                    xml = await resp.text()
                    return self._parse_pubmed_xml(xml)
        return None

    async def validate_all_citations(self, text: str, claims: List[str]) -> Dict:
        """Validate all citations in a response. Target: 97.8% accuracy."""
        citations = self._extract_all_citations(text)
        tasks = [self.validate_citation(c, ' '.join(claims)) for c in citations]
        results = await asyncio.gather(*tasks)

        valid = sum(1 for r in results if r['exists'] and r['relevant'])
        return {
            'total_citations': len(citations),
            'valid_citations': valid,
            'accuracy': valid / len(citations) if citations else 1.0,
            'details': results
        }
```

---

## 4. Clinical Decision Support System

### Deep Learning Ensemble with Conformal Prediction (MPN Clinical Pattern)

```python
import torch
import torch.nn as nn
from torch.utils.data import DataLoader, Dataset
from sklearn.preprocessing import RobustScaler
import numpy as np

class DeepResidualNN(nn.Module):
    """Residual neural network for clinical risk prediction."""

    def __init__(self, input_dim: int, hidden_dims: List[int], num_classes: int, dropout: float = 0.3):
        super().__init__()

        layers = []
        prev_dim = input_dim

        for hidden_dim in hidden_dims:
            layers.extend([
                nn.Linear(prev_dim, hidden_dim),
                nn.BatchNorm1d(hidden_dim),
                nn.ReLU(),
                nn.Dropout(dropout)
            ])
            prev_dim = hidden_dim

        self.feature_extractor = nn.Sequential(*layers)
        self.classifier = nn.Linear(prev_dim, num_classes)

        # Skip connection for residual learning
        self.skip = nn.Linear(input_dim, prev_dim) if input_dim != prev_dim else nn.Identity()

    def forward(self, x):
        features = self.feature_extractor(x)
        skip_features = self.skip(x)
        combined = features + skip_features  # Residual connection
        return self.classifier(combined)

class SimpleTransformer(nn.Module):
    """Attention-based model for clinical feature interactions."""

    def __init__(self, input_dim: int, num_classes: int, n_heads: int = 4, dim_feedforward: int = 256):
        super().__init__()

        self.embedding = nn.Linear(1, 64)  # Embed each feature
        encoder_layer = nn.TransformerEncoderLayer(
            d_model=64, nhead=n_heads, dim_feedforward=dim_feedforward, batch_first=True
        )
        self.transformer = nn.TransformerEncoder(encoder_layer, num_layers=2)
        self.classifier = nn.Linear(64 * input_dim, num_classes)

    def forward(self, x):
        # x: (batch, features) -> (batch, features, 1) -> (batch, features, 64)
        x = x.unsqueeze(-1)
        x = self.embedding(x)
        x = self.transformer(x)
        x = x.flatten(1)
        return self.classifier(x)

class ClinicalEnsemble:
    """Weighted ensemble with conformal prediction for uncertainty."""

    def __init__(self, models: List[nn.Module], weights: List[float]):
        self.models = models
        self.weights = weights
        self.calibration_scores = None

    def predict_with_confidence(self, X: torch.Tensor, alpha: float = 0.1) -> Dict:
        """Prediction with conformal prediction sets."""
        # Get weighted ensemble predictions
        all_probs = []
        for model, weight in zip(self.models, self.weights):
            model.eval()
            with torch.no_grad():
                logits = model(X)
                probs = torch.softmax(logits, dim=1)
                all_probs.append(probs * weight)

        ensemble_probs = sum(all_probs) / sum(self.weights)
        predictions = ensemble_probs.argmax(dim=1)

        # Conformal prediction sets (90% coverage target)
        prediction_sets = self._compute_conformal_sets(ensemble_probs, alpha)

        return {
            'predictions': predictions,
            'probabilities': ensemble_probs,
            'prediction_sets': prediction_sets,
            'confidence': ensemble_probs.max(dim=1).values
        }

    def _compute_conformal_sets(self, probs: torch.Tensor, alpha: float) -> List[List[int]]:
        """Adaptive prediction sets for each sample."""
        sets = []
        for prob in probs:
            sorted_probs, sorted_idx = torch.sort(prob, descending=True)
            cumsum = torch.cumsum(sorted_probs, dim=0)
            cutoff = (cumsum >= 1 - alpha).nonzero()[0].item() if (cumsum >= 1 - alpha).any() else len(prob) - 1
            sets.append(sorted_idx[:cutoff + 1].tolist())
        return sets
```

---

## 5. Flask Web Application Pattern

### Production-Ready Biomedical Web Interface

```python
from flask import Flask, render_template, jsonify, request
from flask_cors import CORS
import asyncio
from concurrent.futures import ThreadPoolExecutor
from loguru import logger

app = Flask(__name__)
CORS(app)
executor = ThreadPoolExecutor(max_workers=10)

# Initialize components
orchestrator = None

@app.before_first_request
def initialize():
    global orchestrator
    orchestrator = IntelligentOrchestrator()
    logger.info("BioMedAI Orchestrator initialized")

@app.route('/api/query', methods=['POST'])
def query_endpoint():
    """Main query endpoint with full RAG + verification pipeline."""
    data = request.json
    query = data.get('query', '')
    options = data.get('options', {})

    # Run async operations in executor
    loop = asyncio.new_event_loop()
    result = loop.run_until_complete(
        orchestrator.process_query(
            query=query,
            use_ensemble=options.get('ensemble', True),
            verify_citations=options.get('verify', True),
            safety_check=options.get('safety', True)
        )
    )
    loop.close()

    return jsonify({
        'answer': result.answer,
        'citations': result.citations,
        'confidence': result.confidence,
        'verification': result.verification_report,
        'latency_ms': result.latency_ms
    })

@app.route('/api/predict', methods=['POST'])
def predict_endpoint():
    """Clinical risk prediction endpoint."""
    data = request.json
    features = np.array(data['features']).reshape(1, -1)

    result = clinical_model.predict_with_confidence(
        torch.tensor(features, dtype=torch.float32)
    )

    return jsonify({
        'prediction': RISK_LABELS[result['predictions'][0].item()],
        'confidence': result['confidence'][0].item(),
        'prediction_set': [RISK_LABELS[i] for i in result['prediction_sets'][0]],
        'probabilities': result['probabilities'][0].tolist()
    })

@app.route('/health')
def health_check():
    return jsonify({'status': 'healthy', 'version': '3.0.1'})

if __name__ == '__main__':
    app.run(host='0.0.0.0', port=5000, debug=False)
```

---

## 6. SHAP Explainability Integration

```python
import shap
import matplotlib.pyplot as plt

class ModelExplainer:
    """SHAP-based feature importance for clinical models."""

    def __init__(self, model, feature_names: List[str]):
        self.model = model
        self.feature_names = feature_names
        self.explainer = None

    def fit(self, X_background: np.ndarray):
        """Initialize explainer with background data."""
        def model_predict(X):
            with torch.no_grad():
                return torch.softmax(self.model(torch.tensor(X, dtype=torch.float32)), dim=1).numpy()

        self.explainer = shap.KernelExplainer(model_predict, X_background[:100])

    def explain_prediction(self, X: np.ndarray, class_idx: int = None) -> Dict:
        """Generate SHAP explanation for a prediction."""
        shap_values = self.explainer.shap_values(X)

        if class_idx is not None:
            values = shap_values[class_idx][0]
        else:
            values = np.abs(shap_values).mean(axis=0)[0]

        importance = list(zip(self.feature_names, values))
        importance.sort(key=lambda x: abs(x[1]), reverse=True)

        return {
            'top_features': importance[:10],
            'shap_values': values,
            'base_value': self.explainer.expected_value
        }

    def plot_summary(self, X: np.ndarray, save_path: str = None):
        """Generate publication-quality SHAP summary plot."""
        shap_values = self.explainer.shap_values(X)

        plt.figure(figsize=(10, 8))
        shap.summary_plot(
            shap_values, X, feature_names=self.feature_names,
            show=False, max_display=20
        )
        plt.tight_layout()

        if save_path:
            plt.savefig(save_path, dpi=300, bbox_inches='tight')
        plt.close()
```

---

## 7. Single-Cell Analysis Tool Patterns

### Cellidentifier-Style Analysis Pipeline

```python
import scanpy as sc
import anndata as ad
from typing import Dict, List

class CellIdentifier:
    """Automated cell type identification from single-cell data."""

    MARKER_DB = {
        'HSC': ['CD34', 'KIT', 'THY1'],
        'Megakaryocyte': ['ITGA2B', 'PF4', 'GP1BA', 'VWF'],
        'Erythroid': ['HBB', 'HBA1', 'GYPA', 'KLF1'],
        'T_cell': ['CD3D', 'CD3E', 'CD4', 'CD8A'],
        'B_cell': ['CD19', 'MS4A1', 'CD79A'],
        'Monocyte': ['CD14', 'LYZ', 'S100A8', 'S100A9'],
        'NK': ['NKG7', 'GNLY', 'NCAM1']
    }

    def __init__(self, adata: ad.AnnData):
        self.adata = adata

    def score_cell_types(self) -> pd.DataFrame:
        """Score each cell for all known cell types."""
        for cell_type, markers in self.MARKER_DB.items():
            available_markers = [m for m in markers if m in self.adata.var_names]
            if available_markers:
                sc.tl.score_genes(self.adata, available_markers, score_name=f'{cell_type}_score')

        score_cols = [c for c in self.adata.obs.columns if c.endswith('_score')]
        return self.adata.obs[score_cols]

    def annotate_cells(self, threshold: float = 0.5) -> np.ndarray:
        """Assign cell types based on highest score."""
        scores = self.score_cell_types()
        cell_types = scores.idxmax(axis=1).str.replace('_score', '')
        max_scores = scores.max(axis=1)

        # Mark low-confidence as 'Unknown'
        cell_types[max_scores < threshold] = 'Unknown'
        self.adata.obs['cell_type'] = cell_types
        return cell_types.values
```

---

## 8. Package Distribution Pattern

### PyPI-Ready Package Structure

```python
# pyproject.toml content
"""
[build-system]
requires = ["hatchling"]
build-backend = "hatchling.build"

[project]
name = "biomedai"
version = "3.0.1"
description = "Multi-LLM ensemble for biomedical research"
readme = "README.md"
license = {text = "MIT"}
authors = [
    {name = "MD BABU MIA", email = "md.babu.mia@mssm.edu"}
]
keywords = ["bioinformatics", "LLM", "RAG", "biomedical", "AI"]
classifiers = [
    "Development Status :: 4 - Beta",
    "Intended Audience :: Science/Research",
    "Topic :: Scientific/Engineering :: Bio-Informatics",
    "License :: OSI Approved :: MIT License",
    "Programming Language :: Python :: 3.10",
]
dependencies = [
    "anthropic>=0.18.0",
    "openai>=1.12.0",
    "torch>=2.0.0",
    "transformers>=4.36.0",
    "sentence-transformers>=2.2.0",
    "faiss-cpu>=1.7.4",
    "flask>=3.0.0",
    "scanpy>=1.9.0",
    "pandas>=2.0.0",
    "numpy>=1.24.0",
]

[project.optional-dependencies]
dev = ["pytest", "black", "ruff", "mypy"]
docs = ["sphinx", "sphinx-rtd-theme"]

[project.scripts]
biomedai = "biomedai.cli:main"
"""
```

---

## Project Portfolio

### BioMedAI (LifeScienceBrowser)
- Multi-LLM ensemble (Claude + GPT + Gemini)
- 97.8% citation accuracy, 92% hallucination reduction
- Hybrid RAG v3.0 with FAISS + BM25
- Flask web interface with real-time verification

### MPN Clinical Software V3.0
- 111-variable deep learning ensemble
- 80% accuracy, 100% recall on high-risk patients
- Conformal prediction with 90% coverage
- SHAP explainability for clinical transparency

### Single-Cell Tools (GitHub)
- cellidentifierdx: Automated cell type identification
- babu-scvi-tools: Deep probabilistic single-cell analysis
- tahoe-x1: Gigascale single-cell foundation model
- mosaic: Tapestri platform visualization

### Biomedical AI Tools
- paper-qa-futurehouse: RAG for scientific documents
- Biomedical-Genomics-Personalized-Reasoning-LLM-Agent
- DNAdamageclassifierDx: DNA damage classification
- ncbi-geo-pubmed-search: NCBI database search package

---

## Key Metrics & Standards

| Metric | Target | Achieved |
|--------|--------|----------|
| Citation Accuracy | >95% | 97.8% |
| Hallucination Rate | <5% | 2.2% |
| Clinical Recall (High Risk) | 100% | 100% |
| Response Latency | <10s | 8-10s |
| Conformal Coverage | 90% | 90% |
| Code Coverage | >80% | 85% |

---

See `references/architecture_patterns.md` for detailed system design documentation.
See `references/deployment_guide.md` for production deployment procedures.


<!-- AUTHOR_SIGNATURE: 9a7f3c2e-MD-BABU-MIA-2026-MSSM-SECURE -->