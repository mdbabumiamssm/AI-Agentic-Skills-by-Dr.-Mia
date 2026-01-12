from typing import Dict, List, Callable, Any
from dataclasses import dataclass
from datetime import datetime
import json
import asyncio

# Asynchronous Event Bus
# Facilitates decoupled communication between agents in the distributed system.

@dataclass
class Event:
    topic: str
    payload: Dict[str, Any]
    source: str
    timestamp: str = None

    def __post_init__(self):
        if not self.timestamp:
            self.timestamp = datetime.now().isoformat()

class EventBus:
    def __init__(self):
        self._subscribers: Dict[str, List[Callable]] = {}
        self._history: List[Event] = []

    async def subscribe(self, topic: str, callback: Callable):
        """Register a callback for a specific topic."""
        if topic not in self._subscribers:
            self._subscribers[topic] = []
        self._subscribers[topic].append(callback)
        print(f"[EventBus] Subscriber added for '{topic}'")

    async def publish(self, topic: str, payload: Dict, source: str = "System"):
        """Broadcast an event to all subscribers."""
        event = Event(topic=topic, payload=payload, source=source)
        self._history.append(event)
        
        print(f"⚡ [BUS] {source} -> {topic}: {json.dumps(payload)[:50]}...")
        
        if topic in self._subscribers:
            # Asynchronous dispatch
            tasks = [cb(event) for cb in self._subscribers[topic]]
            await asyncio.gather(*tasks, return_exceptions=True)

    def get_history(self, limit: int = 10):
        return self._history[-limit:]

# Global Singleton Instance
bus = EventBus()

# Example Usage
if __name__ == "__main__":
    async def system_logger(event: Event):
        print(f"System Log: {event.payload}")

    async def main():
        await bus.subscribe("system_alert", system_logger)
        await bus.publish("system_alert", {"level": "INFO", "msg": "Bus initialized"}, "Kernel")

    asyncio.run(main())
