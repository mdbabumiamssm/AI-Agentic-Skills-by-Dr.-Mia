import asyncio
from typing import Dict, List, Callable, Any
from dataclasses import dataclass
from datetime import datetime
import json

# The "Telegraph Line" of the Wild West
# Asynchronous Event Bus for Decoupled Agent Communication

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
        """Listen for messages on a specific channel."""
        if topic not in self._subscribers:
            self._subscribers[topic] = []
        self._subscribers[topic].append(callback)
        print(f"[EventBus] Subscribed to '{topic}'")

    async def publish(self, topic: str, payload: Dict, source: str = "System"):
        """Broadcast a message to the Wild West."""
        event = Event(topic=topic, payload=payload, source=source)
        self._history.append(event)
        
        print(f"⚡ [BUS] {source} -> {topic}: {json.dumps(payload)[:50]}...")
        
        if topic in self._subscribers:
            # Fire and forget (Wild West style)
            tasks = [cb(event) for cb in self._subscribers[topic]]
            await asyncio.gather(*tasks, return_exceptions=True)

    def get_history(self, limit: int = 10):
        return self._history[-limit:]

# Global Singleton Instance
bus = EventBus()

# Example Usage
if __name__ == "__main__":
    async def sheriff_logger(event: Event):
        print(f"🤠 Sheriff received: {event.payload}")

    async def main():
        await bus.subscribe("outlaw_sighted", sheriff_logger)
        await bus.publish("outlaw_sighted", {"name": "Billy the Kid", "location": "Saloon"}, "Bartender")

    asyncio.run(main())
