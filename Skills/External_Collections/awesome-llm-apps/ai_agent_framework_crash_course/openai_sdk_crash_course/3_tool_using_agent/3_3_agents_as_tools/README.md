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

# Agents as Tools

Demonstrates advanced orchestration patterns where agents are used as tools by other agents.

## 🎯 What This Demonstrates

- **Agent.as_tool()**: Converting agents to tools for orchestration
- **Custom Agent Tools**: Using `@function_tool` with `Runner.run()`
- **Multi-Agent Workflows**: Coordinating multiple specialized agents
- **Custom Configuration**: Per-agent settings like max_turns and run_config

## 🚀 Quick Start

1. **Install OpenAI Agents SDK**:
   ```bash
   pip install openai-agents
   ```

2. **Set up environment**:
   ```bash
   cp ../env.example .env
   # Edit .env and add your OpenAI API key
   ```

3. **Run basic orchestration**:
   ```python
   from agents import Runner
   from agent import root_agent
   
   result = Runner.run_sync(root_agent, "Say 'Hello, how are you?' in Spanish.")
   print(result.final_output)
   ```

4. **Try advanced orchestration**:
   ```python
   from advanced_agent import advanced_orchestrator
   
   result = Runner.run_sync(advanced_orchestrator, "Research the benefits of AI in healthcare.")
   print(result.final_output)
   ```

## 💡 Key Concepts

### Basic Agent Tools (`agent.py`)
- **Agent.as_tool()**: Simple agent-to-tool conversion
- **Translation Orchestration**: Multiple language agents coordinated
- **Tool Naming**: Custom tool names and descriptions

### Advanced Agent Tools (`advanced_agent.py`)
- **@function_tool with Runner.run()**: Custom agent tool implementations
- **Custom Configuration**: Per-run settings (max_turns, temperature)
- **Research-Writing Pipeline**: Complex multi-stage workflows

## 🧪 Available Patterns

### Basic Orchestration
- Spanish translation agent
- French translation agent
- Orchestrator coordinates language tasks

### Advanced Orchestration  
- Research agent for information gathering
- Writing agent for content creation
- Custom tool functions with Runner configuration

## 🔗 Next Steps

- [Function Tools](../3_1_function_tools/README.md) - Custom function tools
- [Built-in Tools](../3_2_builtin_tools/README.md) - SDK provided tools
- [Tutorial 4: Running Agents](../../4_running_agents/README.md) - Advanced execution patterns


<!-- AUTHOR_SIGNATURE: 9a7f3c2e-MD-BABU-MIA-2026-MSSM-SECURE -->