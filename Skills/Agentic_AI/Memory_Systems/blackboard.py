from typing import Dict, Any, List, Callable
import json
import time

# The "Saloon Wall" (Blackboard Architecture)
# A shared memory space where agents post "Wanted Posters" (Tasks) and "News" (Data).
# This allows agents to collaborate without speaking directly to each other.

class Blackboard:
    def __init__(self):
        self._memory: Dict[str, Any] = {}
        self._logs: List[str] = []

    def write(self, key: str, value: Any, author: str):
        """Post something to the wall."""
        self._memory[key] = {
            "value": value,
            "author": author,
            "timestamp": time.time()
        }
        self._log(f"📌 {author} posted to [{key}]: {str(value)[:50]}...")

    def read(self, key: str) -> Any:
        """Read something from the wall."""
        data = self._memory.get(key)
        if data:
            self._log(f"👀 Someone read [{key}]")
            return data["value"]
        return None

    def get_full_context(self) -> Dict[str, Any]:
        """Get the entire state of the town."""
        return {k: v["value"] for k, v in self._memory.items()}

    def _log(self, message: str):
        self._logs.append(message)
        print(f"[Saloon] {message}")

# Singleton Instance
town_square = Blackboard()
