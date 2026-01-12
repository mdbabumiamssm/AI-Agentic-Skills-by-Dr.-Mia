from typing import Dict, Any, List, Callable
import json
import time

# Shared Context Blackboard
# A centralized memory store where agents publish tasks and data.
# Enables decoupled collaboration between autonomous agents.

class Blackboard:
    def __init__(self):
        self._memory: Dict[str, Any] = {}
        self._logs: List[str] = []

    def write(self, key: str, value: Any, author: str):
        """Publish data to the shared context."""
        self._memory[key] = {
            "value": value,
            "author": author,
            "timestamp": time.time()
        }
        self._log(f"📌 {author} published to [{key}]: {str(value)[:50]}...")

    def read(self, key: str) -> Any:
        """Retrieve data from the shared context."""
        data = self._memory.get(key)
        if data:
            self._log(f"👀 Access read on [{key}]")
            return data["value"]
        return None

    def get_full_context(self) -> Dict[str, Any]:
        """Retrieve the complete system state."""
        return {k: v["value"] for k, v in self._memory.items()}

    def _log(self, message: str):
        self._logs.append(message)
        print(f"[Blackboard] {message}")

# Global Instance
shared_context = Blackboard()
