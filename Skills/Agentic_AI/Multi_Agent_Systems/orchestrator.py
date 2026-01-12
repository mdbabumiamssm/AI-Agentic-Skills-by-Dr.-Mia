import asyncio
from typing import List, Dict
from dataclasses import dataclass

# Multi-Agent Orchestrator (Supervisor Pattern)
# Inspired by LangGraph / CrewAI
# Coordinates specialized agents to solve complex bio-tasks.

@dataclass
class AgentMessage:
    sender: str
    content: str
    required_action: str = None

class BaseAgent:
    def __init__(self, name: str, role: str):
        self.name = name
        self.role = role

    async def process(self, message: AgentMessage) -> AgentMessage:
        raise NotImplementedError

class ResearchAgent(BaseAgent):
    async def process(self, message: AgentMessage) -> AgentMessage:
        return AgentMessage(
            sender=self.name,
            content=f"Researched '{message.content}'. Found 3 relevant papers on PubMed.",
            required_action="review"
        )

class ReviewAgent(BaseAgent):
    async def process(self, message: AgentMessage) -> AgentMessage:
        return AgentMessage(
            sender=self.name,
            content=f"Reviewed findings. Paper 2 seems most relevant to the query.",
            required_action="complete"
        )

class Supervisor:
    def __init__(self):
        self.agents = {
            "researcher": ResearchAgent("Dr. Search", "Literature Search"),
            "reviewer": ReviewAgent("Dr. Crit", "Quality Control")
        }

    async def run_workflow(self, task: str):
        print(f"--- Starting Workflow: {task} ---")
        
        # Step 1: Assign to Researcher
        msg = AgentMessage(sender="User", content=task)
        response_1 = await self.agents["researcher"].process(msg)
        print(f"[{response_1.sender}]: {response_1.content}")

        # Step 2: Supervisor Decision Logic (Router)
        if response_1.required_action == "review":
            response_2 = await self.agents["reviewer"].process(response_1)
            print(f"[{response_2.sender}]: {response_2.content}")
            
        print("--- Workflow Complete ---")

if __name__ == "__main__":
    supervisor = Supervisor()
    asyncio.run(supervisor.run_workflow("Find latest EGFR inhibitors"))
