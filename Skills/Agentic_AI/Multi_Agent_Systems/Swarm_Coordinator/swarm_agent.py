import json
from typing import List, Dict, Any

class SwarmCoordinator:
    """
    A lightweight Multi-Agent Swarm Coordinator.
    This pattern allows a master agent to distribute tasks to specialized sub-agents,
    collect their responses, and synthesize a final output.
    This represents the 'Plan-and-Solve' or 'Debate' architecture trending in 2026.
    """

    def __init__(self, name: str = "Chief_Orchestrator"):
        self.name = name
        self.agents = {}

    def register_agent(self, agent_name: str, system_prompt: str):
        """Registers a sub-agent with a specific role."""
        self.agents[agent_name] = system_prompt
        print(f"[{self.name}] Registered sub-agent: {agent_name}")

    def run_swarm(self, task: str) -> Dict[str, Any]:
        """
        Orchestrates the execution of a task across the swarm.
        In a real deployment, this would asynchronously call LLMs for each agent.
        """
        print(f"\n[{self.name}] Initiating swarm execution for task: '{task}'")
        
        results = {}
        for agent_name, prompt in self.agents.items():
            # Mocking the LLM execution for the sub-agent
            print(f" -> Invoking {agent_name}...")
            results[agent_name] = self._mock_agent_response(agent_name, prompt, task)
            
        print(f"[{self.name}] Synthesizing final response...")
        final_synthesis = self._synthesize_results(task, results)
        
        return {
            "task": task,
            "sub_agent_results": results,
            "final_synthesis": final_synthesis
        }

    def _mock_agent_response(self, agent_name: str, role: str, task: str) -> str:
        """Simulates an LLM response based on the agent's role."""
        if "critic" in agent_name.lower():
            return f"Critique: The proposed approach to '{task}' lacks fallback mechanisms."
        elif "researcher" in agent_name.lower():
            return f"Research: Found 5 relevant papers addressing '{task}'."
        elif "coder" in agent_name.lower():
            return f"Code: Initialized repository structure for '{task}'."
        return f"Response to '{task}' from perspective of {role}."

    def _synthesize_results(self, task: str, results: Dict[str, str]) -> str:
        """Simulates a synthesizer LLM merging the sub-agent outputs."""
        synthesis = f"Final Executive Summary for: {task}\n"
        for agent, res in results.items():
            synthesis += f"- {agent}: {res}\n"
        return synthesis

if __name__ == "__main__":
    swarm = SwarmCoordinator()
    swarm.register_agent("AI_Researcher", "You are an expert AI researcher finding state-of-the-art papers.")
    swarm.register_agent("Code_Architect", "You are a senior software engineer designing scalable systems.")
    swarm.register_agent("Security_Critic", "You are a red-team hacker looking for vulnerabilities.")
    
    output = swarm.run_swarm("Design a self-improving AI agent framework")
    print("\n--- FINAL OUTPUT ---")
    print(output["final_synthesis"])
