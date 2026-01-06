from typing import List, Dict, Optional, Literal
from dataclasses import dataclass
import json

@dataclass
class AgentOutput:
    sender: str
    content: str
    next_action: Optional[str] = None

class BaseWorker:
    """
    A specialized worker agent (e.g., Coder, Researcher).
    """
    def __init__(self, name: str, description: str):
        self.name = name
        self.description = description

    def work(self, task: str) -> str:
        # Mock logic - replace with actual LLM call
        return f"[{self.name}]: I have processed the task '{task}' using my skills in {self.description}."

class SupervisorAgent:
    """
    The Orchestrator. It does not do the work; it delegates it.
    It outputs structured JSON to decide 'who acts next'.
    """
    def __init__(self, workers: List[BaseWorker]):
        self.workers = {w.name: w for w in workers}
        self.worker_descriptions = "\n".join([f"- {w.name}: {w.description}" for w in workers])

    def decide_next_step(self, objective: str, history: List[str]) -> Dict:
        """
        Determines the next worker or FINISH.
        Simulates an LLM call that returns JSON.
        """
        # Logic simulation based on keyword matching (Mocking the LLM router)
        last_msg = history[-1] if history else ""
        
        if "code" in objective.lower() and "Coder" in self.workers and "Coder" not in last_msg:
            return {"next": "Coder", "instruction": "Write the Python script for this."}
        elif "review" in objective.lower() and "Reviewer" in self.workers and "Reviewer" not in last_msg:
            return {"next": "Reviewer", "instruction": "Check the code for bugs."}
        else:
            return {"next": "FINISH", "instruction": "Task appears complete."}

    def run_mission(self, objective: str):
        print(f"--- Supervisor Mission: {objective} ---")
        history = []
        max_steps = 5
        
        for i in range(max_steps):
            decision = self.decide_next_step(objective, history)
            next_agent_name = decision["next"]
            instruction = decision["instruction"]
            
            print(f"\n[Supervisor]: Router decided -> {next_agent_name} ('{instruction}')")
            
            if next_agent_name == "FINISH":
                print("--- Mission Accomplished ---")
                break
                
            if next_agent_name in self.workers:
                worker = self.workers[next_agent_name]
                result = worker.work(instruction)
                print(result)
                history.append(f"{next_agent_name}: {result}")
            else:
                print(f"Error: Unknown worker {next_agent_name}")
                break

if __name__ == "__main__":
    # Define the "Digital Coworkers"
    coder = BaseWorker("Coder", "Writes Python software and scripts.")
    reviewer = BaseWorker("Reviewer", "Reviews code for security and performance.")
    
    # Initialize Manager
    supervisor = SupervisorAgent([coder, reviewer])
    
    # Run a mission
    supervisor.run_mission("I need a python script to parse a CSV, then review it.")
