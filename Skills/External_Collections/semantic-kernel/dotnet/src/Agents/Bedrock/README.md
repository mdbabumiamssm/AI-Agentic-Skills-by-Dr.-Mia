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

# Amazon Bedrock AI Agents in Semantic Kernel

## Overview

AWS Bedrock Agents is a managed service that allows users to stand up and run AI agents in the AWS cloud quickly.

## Tools/Functions

Bedrock Agents allow the use of tools via [action groups](https://docs.aws.amazon.com/bedrock/latest/userguide/agents-action-create.html).

The integration of Bedrock Agents with Semantic Kernel allows users to register kernel functions as tools in Bedrock Agents.

## Enable code interpretation

Bedrock Agents can write and execute code via a feature known as [code interpretation](https://docs.aws.amazon.com/bedrock/latest/userguide/agents-code-interpretation.html) similar to what OpenAI also offers.

## Enable user input

Bedrock Agents can [request user input](https://docs.aws.amazon.com/bedrock/latest/userguide/agents-user-input.html) in case of missing information to invoke a tool. When this is enabled, the agent will prompt the user for the missing information. When this is disabled, the agent will guess the missing information.

## Knowledge base

Bedrock Agents can leverage data saved on AWS to perform RAG tasks, this is referred to as the [knowledge base](https://docs.aws.amazon.com/bedrock/latest/userguide/agents-kb-add.html) in AWS.

## Multi-agent

Bedrock Agents support [multi-agent workflows](https://docs.aws.amazon.com/bedrock/latest/userguide/agents-multi-agent-collaboration.html) for more complex tasks. However, it employs a different pattern than what we have in Semantic Kernel, thus this is not supported in the current integration.


<!-- AUTHOR_SIGNATURE: 9a7f3c2e-MD-BABU-MIA-2026-MSSM-SECURE -->