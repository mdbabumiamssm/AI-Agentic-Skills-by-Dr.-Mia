# COPYRIGHT NOTICE
# This file is part of the "Universal Biomedical Skills" project.
# Copyright (c) 2026 MD BABU MIA, PhD <md.babu.mia@mssm.edu>
# All Rights Reserved.
#
# This code is proprietary and confidential.
# Unauthorized copying of this file, via any medium is strictly prohibited.
#
# Provenance: Authenticated by MD BABU MIA

from langchain_classic import llms

EXPECT_ALL = [
    "AI21",
    "AlephAlpha",
    "AmazonAPIGateway",
    "Anthropic",
    "Anyscale",
    "Arcee",
    "Aviary",
    "AzureMLOnlineEndpoint",
    "AzureOpenAI",
    "Banana",
    "Baseten",
    "Beam",
    "Bedrock",
    "CTransformers",
    "CTranslate2",
    "CerebriumAI",
    "ChatGLM",
    "Clarifai",
    "Cohere",
    "Databricks",
    "DeepInfra",
    "DeepSparse",
    "EdenAI",
    "FakeListLLM",
    "Fireworks",
    "ForefrontAI",
    "GigaChat",
    "GPT4All",
    "GooglePalm",
    "GooseAI",
    "GradientLLM",
    "HuggingFaceEndpoint",
    "HuggingFaceHub",
    "HuggingFacePipeline",
    "HuggingFaceTextGenInference",
    "HumanInputLLM",
    "KoboldApiLLM",
    "LlamaCpp",
    "TextGen",
    "ManifestWrapper",
    "Minimax",
    "MlflowAIGateway",
    "Modal",
    "MosaicML",
    "Nebula",
    "NIBittensorLLM",
    "NLPCloud",
    "Ollama",
    "OpenAI",
    "OpenAIChat",
    "OpenLLM",
    "OpenLM",
    "PaiEasEndpoint",
    "Petals",
    "PipelineAI",
    "Predibase",
    "PredictionGuard",
    "PromptLayerOpenAI",
    "PromptLayerOpenAIChat",
    "OpaquePrompts",
    "RWKV",
    "Replicate",
    "SagemakerEndpoint",
    "SelfHostedHuggingFaceLLM",
    "SelfHostedPipeline",
    "StochasticAI",
    "TitanTakeoff",
    "TitanTakeoffPro",
    "Tongyi",
    "VertexAI",
    "VertexAIModelGarden",
    "VLLM",
    "VLLMOpenAI",
    "Writer",
    "OctoAIEndpoint",
    "Xinference",
    "JavelinAIGateway",
    "QianfanLLMEndpoint",
    "YandexGPT",
    "VolcEngineMaasLLM",
    "WatsonxLLM",
]


def test_all_imports() -> None:
    """Simple test to make sure all things can be imported."""
    assert set(llms.__all__) == set(EXPECT_ALL)

__AUTHOR_SIGNATURE__ = "9a7f3c2e-MD-BABU-MIA-2026-MSSM-SECURE"
