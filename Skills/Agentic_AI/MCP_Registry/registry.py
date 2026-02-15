from typing import List, Dict

class MCPRegistry:
    def __init__(self):
        self.active_servers = {}
        self.known_servers = {
            "github-mcp": {"command": "npx -y @modelcontextprotocol/server-github"},
            "postgres-mcp": {"command": "npx -y @modelcontextprotocol/server-postgres"},
            "filesystem-mcp": {"command": "npx -y @modelcontextprotocol/server-filesystem"}
        }

    def start_server(self, name: str):
        if name in self.known_servers:
            print(f"🔗 [MCP Registry] Starting {name}...")
            self.active_servers[name] = "Running"
            return True
        print(f"❌ [MCP Registry] Server {name} not found.")
        return False

    def get_tools(self, name: str) -> List[str]:
        if name == "github-mcp":
            return ["search_repositories", "read_file", "create_issue"]
        elif name == "filesystem-mcp":
            return ["ls", "cat", "grep"]
        return []

if __name__ == "__main__":
    reg = MCPRegistry()
    reg.start_server("github-mcp")
    print(reg.get_tools("github-mcp"))
