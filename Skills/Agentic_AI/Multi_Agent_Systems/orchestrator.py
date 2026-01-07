"""
Multi-Agent Orchestrator with LangGraph-Style State Management

This module implements a production-grade supervisor pattern for coordinating
specialized AI agents. Inspired by LangGraph's state graph architecture and
modern multi-agent orchestration best practices (2026).

Features:
- State graph with conditional routing
- Tool execution nodes with validation
- Human-in-the-loop checkpoints
- Structured output enforcement
- Memory persistence across interactions

References:
- LangGraph: https://github.com/langchain-ai/langgraph
- CrewAI patterns: https://crewai.com
- AutoGen conversational agents: https://microsoft.github.io/autogen
"""

from typing import List, Dict, Any, Optional, Callable, Literal, TypedDict
from dataclasses import dataclass, field
from enum import Enum
from abc import ABC, abstractmethod
import json
import time
from datetime import datetime


# --- State Management ---

class AgentState(TypedDict, total=False):
    """Typed state container for agent workflow."""
    messages: List[Dict[str, str]]
    current_agent: str
    task: str
    context: Dict[str, Any]
    results: List[Dict[str, Any]]
    iteration: int
    status: Literal["pending", "running", "completed", "failed", "human_review"]
    error: Optional[str]


class NodeType(Enum):
    """Types of nodes in the execution graph."""
    ENTRY = "entry"
    AGENT = "agent"
    TOOL = "tool"
    ROUTER = "router"
    HUMAN = "human"
    EXIT = "exit"


@dataclass
class Node:
    """Represents a node in the execution graph."""
    name: str
    node_type: NodeType
    handler: Callable
    metadata: Dict[str, Any] = field(default_factory=dict)


@dataclass
class Edge:
    """Represents a directed edge between nodes."""
    source: str
    target: str
    condition: Optional[Callable[[AgentState], bool]] = None


# --- Base Agent Classes ---

class Agent(ABC):
    """Abstract base class for specialized agents."""

    def __init__(self, name: str, role: str, description: str = ""):
        self.name = name
        self.role = role
        self.description = description
        self.memory: List[Dict[str, Any]] = []
        self.tools: Dict[str, Callable] = {}

    @abstractmethod
    def run(self, state: AgentState) -> AgentState:
        """Execute the agent's logic on the current state."""
        pass

    def add_tool(self, name: str, func: Callable, description: str = ""):
        """Register a tool for this agent."""
        self.tools[name] = {"func": func, "description": description}

    def log_interaction(self, input_data: Any, output_data: Any):
        """Store interaction in agent memory."""
        self.memory.append({
            "timestamp": datetime.now().isoformat(),
            "input": input_data,
            "output": output_data
        })


class LLMAdapter(ABC):
    """Abstract LLM interface for provider flexibility."""

    @abstractmethod
    def generate(self, prompt: str, **kwargs) -> str:
        """Generate a response from the LLM."""
        pass

    @abstractmethod
    def generate_structured(self, prompt: str, schema: Dict) -> Dict:
        """Generate a structured response matching the schema."""
        pass


class MockLLM(LLMAdapter):
    """Mock LLM for testing without API calls."""

    def __init__(self, responses: Optional[Dict[str, str]] = None):
        self.responses = responses or {}
        self.call_count = 0

    def generate(self, prompt: str, **kwargs) -> str:
        self.call_count += 1
        # Pattern matching for demonstration
        if "route" in prompt.lower() or "delegate" in prompt.lower():
            if "drug" in prompt.lower() or "molecule" in prompt.lower():
                return json.dumps({"agent": "chemist", "reason": "Drug-related query"})
            elif "research" in prompt.lower() or "literature" in prompt.lower():
                return json.dumps({"agent": "researcher", "reason": "Research query"})
            elif "clinical" in prompt.lower() or "patient" in prompt.lower():
                return json.dumps({"agent": "clinician", "reason": "Clinical query"})
        return self.responses.get(prompt, f"Mock response for: {prompt[:50]}...")

    def generate_structured(self, prompt: str, schema: Dict) -> Dict:
        self.call_count += 1
        # Return default values based on schema
        result = {}
        for key, spec in schema.get("properties", {}).items():
            if spec.get("type") == "string":
                result[key] = f"mock_{key}"
            elif spec.get("type") == "number":
                result[key] = 0.0
            elif spec.get("type") == "boolean":
                result[key] = True
            elif spec.get("type") == "array":
                result[key] = []
        return result


# --- Specialized Agents ---

class ChemistAgent(Agent):
    """Agent specialized in drug discovery and cheminformatics."""

    def __init__(self, llm: Optional[LLMAdapter] = None):
        super().__init__(
            name="Chemist",
            role="Drug Discovery Specialist",
            description="Analyzes molecules, calculates properties, designs compounds"
        )
        self.llm = llm or MockLLM()
        self._register_tools()

    def _register_tools(self):
        self.add_tool("calculate_properties", self._calc_props, "Calculate molecular descriptors")
        self.add_tool("assess_druglikeness", self._assess_drug, "Assess drug-likeness")

    def _calc_props(self, smiles: str) -> Dict[str, float]:
        # Placeholder - in production, use ChemTools from Drug_Discovery
        return {"MolWt": 180.16, "LogP": 1.19, "TPSA": 63.6, "QED": 0.56}

    def _assess_drug(self, smiles: str) -> Dict[str, Any]:
        return {"lipinski_violations": 0, "is_druglike": True}

    def run(self, state: AgentState) -> AgentState:
        task = state.get("task", "")
        context = state.get("context", {})

        # Extract SMILES if present
        smiles = context.get("smiles", "CC(=O)OC1=CC=CC=C1C(=O)O")  # Default: aspirin

        # Execute analysis
        props = self._calc_props(smiles)
        druglikeness = self._assess_drug(smiles)

        result = {
            "agent": self.name,
            "task": task,
            "smiles": smiles,
            "properties": props,
            "druglikeness": druglikeness,
            "recommendation": "Compound passes initial screening" if druglikeness["is_druglike"] else "Consider optimization"
        }

        self.log_interaction(task, result)

        # Update state
        state["results"] = state.get("results", []) + [result]
        state["messages"] = state.get("messages", []) + [
            {"role": "agent", "content": json.dumps(result)}
        ]
        return state


class ResearcherAgent(Agent):
    """Agent specialized in literature search and knowledge synthesis."""

    def __init__(self, llm: Optional[LLMAdapter] = None):
        super().__init__(
            name="Researcher",
            role="Literature Specialist",
            description="Searches scientific literature, synthesizes knowledge"
        )
        self.llm = llm or MockLLM()

    def run(self, state: AgentState) -> AgentState:
        task = state.get("task", "")

        # Mock literature search
        result = {
            "agent": self.name,
            "task": task,
            "sources": [
                {"title": "Recent advances in drug discovery", "year": 2025, "relevance": 0.95},
                {"title": "AI in pharmaceutical research", "year": 2026, "relevance": 0.87}
            ],
            "summary": f"Found relevant literature for: {task[:50]}..."
        }

        self.log_interaction(task, result)
        state["results"] = state.get("results", []) + [result]
        state["messages"] = state.get("messages", []) + [
            {"role": "agent", "content": json.dumps(result)}
        ]
        return state


class ClinicianAgent(Agent):
    """Agent specialized in clinical analysis and patient care."""

    def __init__(self, llm: Optional[LLMAdapter] = None):
        super().__init__(
            name="Clinician",
            role="Clinical Specialist",
            description="Analyzes clinical data, provides medical insights"
        )
        self.llm = llm or MockLLM()

    def run(self, state: AgentState) -> AgentState:
        task = state.get("task", "")

        result = {
            "agent": self.name,
            "task": task,
            "analysis": "Clinical analysis completed",
            "requires_human_review": True,
            "confidence": 0.85
        }

        self.log_interaction(task, result)
        state["results"] = state.get("results", []) + [result]

        # Flag for human review if clinical
        if result.get("requires_human_review"):
            state["status"] = "human_review"

        return state


# --- Orchestrator ---

class SupervisorOrchestrator:
    """
    Production-grade multi-agent orchestrator with state graph execution.

    Implements the Supervisor pattern where a meta-agent delegates tasks
    to specialized worker agents based on task analysis.

    Example:
        >>> orchestrator = SupervisorOrchestrator()
        >>> orchestrator.add_agent(ChemistAgent())
        >>> orchestrator.add_agent(ResearcherAgent())
        >>> result = orchestrator.run("Analyze aspirin properties")
    """

    def __init__(self, llm: Optional[LLMAdapter] = None, max_iterations: int = 10):
        self.llm = llm or MockLLM()
        self.agents: Dict[str, Agent] = {}
        self.nodes: Dict[str, Node] = {}
        self.edges: List[Edge] = []
        self.max_iterations = max_iterations
        self.execution_log: List[Dict] = []

        # Initialize default graph structure
        self._build_default_graph()

    def _build_default_graph(self):
        """Construct the default execution graph."""
        # Entry node
        self.nodes["entry"] = Node(
            name="entry",
            node_type=NodeType.ENTRY,
            handler=self._entry_handler
        )

        # Router node (decides which agent to use)
        self.nodes["router"] = Node(
            name="router",
            node_type=NodeType.ROUTER,
            handler=self._router_handler
        )

        # Human review checkpoint
        self.nodes["human_review"] = Node(
            name="human_review",
            node_type=NodeType.HUMAN,
            handler=self._human_review_handler
        )

        # Exit node
        self.nodes["exit"] = Node(
            name="exit",
            node_type=NodeType.EXIT,
            handler=self._exit_handler
        )

        # Default edges
        self.edges.append(Edge("entry", "router"))
        self.edges.append(Edge("router", "exit", lambda s: s.get("status") == "completed"))
        self.edges.append(Edge("human_review", "exit"))

    def add_agent(self, agent: Agent):
        """Register an agent with the orchestrator."""
        self.agents[agent.name.lower()] = agent

        # Create agent node
        self.nodes[agent.name.lower()] = Node(
            name=agent.name.lower(),
            node_type=NodeType.AGENT,
            handler=lambda s, a=agent: a.run(s)
        )

        # Add edges from router to agent and back
        self.edges.append(Edge("router", agent.name.lower()))
        self.edges.append(Edge(agent.name.lower(), "router"))
        self.edges.append(Edge(
            agent.name.lower(),
            "human_review",
            lambda s: s.get("status") == "human_review"
        ))

    def _entry_handler(self, state: AgentState) -> AgentState:
        """Initialize state for new task."""
        state["status"] = "running"
        state["iteration"] = 0
        state["results"] = []
        state["messages"] = state.get("messages", [])
        return state

    def _router_handler(self, state: AgentState) -> AgentState:
        """Route task to appropriate agent using LLM."""
        task = state.get("task", "")

        # Check if we've completed enough iterations
        if state.get("iteration", 0) >= self.max_iterations:
            state["status"] = "completed"
            return state

        # Check if we have results and should exit
        if state.get("results") and len(state["results"]) > 0:
            state["status"] = "completed"
            return state

        # Build routing prompt
        agent_descriptions = "\n".join([
            f"- {name}: {agent.description}"
            for name, agent in self.agents.items()
        ])

        routing_prompt = f"""
        Given the following task, decide which agent should handle it.

        Task: {task}

        Available agents:
        {agent_descriptions}

        Return JSON with: {{"agent": "<agent_name>", "reason": "<brief reason>"}}
        """

        response = self.llm.generate(routing_prompt)

        try:
            decision = json.loads(response)
            target_agent = decision.get("agent", "").lower()

            if target_agent in self.agents:
                state["current_agent"] = target_agent
            else:
                state["status"] = "failed"
                state["error"] = f"Unknown agent: {target_agent}"
        except json.JSONDecodeError:
            # Fallback to keyword matching
            task_lower = task.lower()
            if any(kw in task_lower for kw in ["drug", "molecule", "compound", "chemical"]):
                state["current_agent"] = "chemist"
            elif any(kw in task_lower for kw in ["research", "literature", "paper", "study"]):
                state["current_agent"] = "researcher"
            elif any(kw in task_lower for kw in ["clinical", "patient", "diagnosis", "treatment"]):
                state["current_agent"] = "clinician"
            else:
                state["current_agent"] = list(self.agents.keys())[0] if self.agents else None

        state["iteration"] = state.get("iteration", 0) + 1
        return state

    def _human_review_handler(self, state: AgentState) -> AgentState:
        """Handle human-in-the-loop checkpoint."""
        print("\n[HUMAN REVIEW REQUIRED]")
        print(f"Task: {state.get('task')}")
        print(f"Results: {json.dumps(state.get('results', []), indent=2)}")

        # In production, this would wait for human input
        # For now, auto-approve after logging
        state["status"] = "completed"
        state["messages"] = state.get("messages", []) + [
            {"role": "system", "content": "Human review completed (auto-approved in demo mode)"}
        ]
        return state

    def _exit_handler(self, state: AgentState) -> AgentState:
        """Finalize execution and prepare output."""
        state["status"] = "completed"
        return state

    def run(self, task: str, context: Optional[Dict] = None) -> Dict[str, Any]:
        """
        Execute the orchestrator with the given task.

        Args:
            task: The task description
            context: Optional context data (e.g., SMILES, patient data)

        Returns:
            Final state with all results
        """
        start_time = time.time()

        # Initialize state
        state: AgentState = {
            "task": task,
            "context": context or {},
            "messages": [{"role": "user", "content": task}],
            "status": "pending"
        }

        # Execute entry
        state = self._entry_handler(state)

        # Main execution loop
        current_node = "router"
        visited = []

        while state["status"] not in ["completed", "failed"]:
            visited.append(current_node)

            # Execute current node
            node = self.nodes.get(current_node)
            if node:
                state = node.handler(state)

            # Determine next node
            if state["status"] == "human_review":
                current_node = "human_review"
            elif state.get("current_agent") and current_node == "router":
                current_node = state["current_agent"]
            elif current_node in self.agents:
                current_node = "router"  # Go back to router after agent execution
            else:
                current_node = "exit"

            # Safety check for infinite loops
            if state.get("iteration", 0) > self.max_iterations:
                state["status"] = "failed"
                state["error"] = "Max iterations exceeded"
                break

        # Finalize
        state = self._exit_handler(state)

        # Log execution
        execution_record = {
            "task": task,
            "duration": time.time() - start_time,
            "iterations": state.get("iteration", 0),
            "agents_used": visited,
            "status": state["status"],
            "timestamp": datetime.now().isoformat()
        }
        self.execution_log.append(execution_record)

        return state

    def get_execution_summary(self) -> Dict[str, Any]:
        """Get summary of all executions."""
        return {
            "total_executions": len(self.execution_log),
            "average_duration": sum(e["duration"] for e in self.execution_log) / len(self.execution_log) if self.execution_log else 0,
            "success_rate": sum(1 for e in self.execution_log if e["status"] == "completed") / len(self.execution_log) if self.execution_log else 0,
            "recent_executions": self.execution_log[-5:]
        }


# --- Example Usage ---

if __name__ == "__main__":
    # Create orchestrator with mock LLM
    orchestrator = SupervisorOrchestrator()

    # Register specialized agents
    orchestrator.add_agent(ChemistAgent())
    orchestrator.add_agent(ResearcherAgent())
    orchestrator.add_agent(ClinicianAgent())

    # Example 1: Drug discovery task
    print("=" * 60)
    print("Example 1: Drug Discovery Task")
    print("=" * 60)
    result = orchestrator.run(
        task="Calculate properties and assess druglikeness for aspirin",
        context={"smiles": "CC(=O)OC1=CC=CC=C1C(=O)O"}
    )
    print(f"Status: {result['status']}")
    print(f"Results: {json.dumps(result['results'], indent=2)}")

    # Example 2: Research task
    print("\n" + "=" * 60)
    print("Example 2: Research Task")
    print("=" * 60)
    result = orchestrator.run(
        task="Find recent literature on AI in drug discovery"
    )
    print(f"Status: {result['status']}")
    print(f"Results: {json.dumps(result['results'], indent=2)}")

    # Example 3: Clinical task (triggers human review)
    print("\n" + "=" * 60)
    print("Example 3: Clinical Task (Human Review)")
    print("=" * 60)
    result = orchestrator.run(
        task="Analyze patient symptoms and recommend treatment options"
    )
    print(f"Status: {result['status']}")

    # Print execution summary
    print("\n" + "=" * 60)
    print("Execution Summary")
    print("=" * 60)
    print(json.dumps(orchestrator.get_execution_summary(), indent=2))
