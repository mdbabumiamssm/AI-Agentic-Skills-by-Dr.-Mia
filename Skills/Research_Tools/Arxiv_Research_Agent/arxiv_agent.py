import urllib.request
import urllib.parse
import xml.etree.ElementTree as ET
from typing import List, Dict, Any

class ArxivResearchAgent:
    """
    An agentic skill for querying the arXiv API to find scientific papers.
    Useful for research workflows, academic literature reviews, and summarizing state-of-the-art.
    """
    
    BASE_URL = "http://export.arxiv.org/api/query?"
    
    def search_papers(self, query: str, max_results: int = 5) -> List[Dict[str, Any]]:
        """
        Searches arXiv for papers matching the query.
        
        Args:
            query (str): The search query (e.g., 'machine learning healthcare').
            max_results (int): Maximum number of results to return.
            
        Returns:
            List of dictionaries containing paper metadata (title, authors, summary, link).
        """
        # Encode the query parameter
        search_query = urllib.parse.quote(f"all:{query}")
        url = f"{self.BASE_URL}search_query={search_query}&start=0&max_results={max_results}"
        
        try:
            with urllib.request.urlopen(url) as response:
                if response.status == 200:
                    xml_data = response.read()
                    return self._parse_arxiv_response(xml_data)
                else:
                    print(f"Error: arXiv API returned status code {response.status}")
                    return []
        except Exception as e:
            print(f"Error fetching data from arXiv: {e}")
            return []

    def _parse_arxiv_response(self, xml_data: bytes) -> List[Dict[str, Any]]:
        """
        Parses the Atom XML response from arXiv.
        """
        ns = {'atom': 'http://www.w3.org/2005/Atom'}
        root = ET.fromstring(xml_data)
        
        papers = []
        for entry in root.findall('atom:entry', ns):
            title = entry.find('atom:title', ns).text.strip()
            summary = entry.find('atom:summary', ns).text.strip()
            link = entry.find('atom:id', ns).text.strip()
            
            authors = []
            for author in entry.findall('atom:author', ns):
                name = author.find('atom:name', ns).text.strip()
                authors.append(name)
                
            papers.append({
                "title": title,
                "authors": authors,
                "summary": summary.replace('\n', ' '),
                "link": link
            })
            
        return papers

if __name__ == "__main__":
    agent = ArxivResearchAgent()
    print("Searching for 'Agentic AI workflows'...")
    results = agent.search_papers("Agentic AI workflows", max_results=2)
    for i, res in enumerate(results, 1):
        print(f"\n[{i}] {res['title']}")
        print(f"Authors: {', '.join(res['authors'])}")
        print(f"Link: {res['link']}")
