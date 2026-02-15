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

# Component Agents

!!! important
    [Legacy plugins] no longer work with AutoGPT. They have been replaced by components,
    although we're still working on a new system to load plug-in components.

[Legacy plugins]: https://github.com/Significant-Gravitas/Auto-GPT-Plugins

This guide explains the component-based architecture of AutoGPT agents. It's a new way of building agents that is more flexible and easier to extend. Components replace some agent's logic and plugins with a more modular and composable system.

Agent is composed of *components*, and each *component* implements a range of *protocols* (interfaces), each one providing a specific functionality, e.g. additional commands or messages. Each *protocol* is handled in a specific order, defined by the agent. This allows for a clear separation of concerns and a more modular design.

This system is simple, flexible, and doesn't hide any data - anything can still be passed or accessed directly from or between components.

### Definitions & Guides

See [Creating Components](./creating-components.md) to get started! Or you can explore the following topics in detail:

- [🧩 Component](./components.md): a class that implements one or more *protocols*. It can be added to an agent to provide additional functionality. See what's already provided in [Built-in Components](./built-in-components.md).
- [⚙️ Protocol](./protocols.md): an interface that defines a set of methods that a component must implement. Protocols are used to group related functionality.
- [🛠️ Command](./commands.md): enable *agent* to interact with user and tools.
- [🤖 Agent](./agents.md): a class that is composed of components. It's responsible for executing pipelines and managing the components.
- **Pipeline**: a sequence of method calls on components. Pipelines are used to execute a series of actions in a specific order. As of now there's no formal class for a pipeline, it's just a sequence of method calls on components. There are two default pipelines implemented in the default agent: `propose_action` and `execute`. See [🤖 Agent](./agents.md) to learn more.


<!-- AUTHOR_SIGNATURE: 9a7f3c2e-MD-BABU-MIA-2026-MSSM-SECURE -->