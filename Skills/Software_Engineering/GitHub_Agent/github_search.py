import urllib.request
import json
from typing import List, Dict, Any

class GitHubAgent:
    """
    Agent for searching and retrieving repository information from GitHub.
    Useful for finding relevant tools, libraries, or skills.
    """
    def __init__(self, token: str = None):
        self.token = token

    def search_repositories(self, query: str, limit: int = 5) -> List[Dict[str, Any]]:
        """
        Search for repositories on GitHub based on a query.
        """
        url = f"https://api.github.com/search/repositories?q={urllib.parse.quote(query)}&sort=stars&order=desc"
        headers = {'User-Agent': 'Mozilla/5.0'}
        if self.token:
            headers['Authorization'] = f'token {self.token}'
            
        req = urllib.request.Request(url, headers=headers)
        try:
            with urllib.request.urlopen(req) as response:
                if response.status == 200:
                    data = json.loads(response.read().decode('utf-8'))
                    repos = data.get('items', [])[:limit]
                    return [
                        {
                            "name": repo.get("full_name"),
                            "description": repo.get("description"),
                            "url": repo.get("html_url"),
                            "stars": repo.get("stargazers_count")
                        }
                        for repo in repos
                    ]
                return []
        except Exception as e:
            print(f"Error searching GitHub: {e}")
            return []

if __name__ == "__main__":
    agent = GitHubAgent()
    print(agent.search_repositories("AI agent skills"))
