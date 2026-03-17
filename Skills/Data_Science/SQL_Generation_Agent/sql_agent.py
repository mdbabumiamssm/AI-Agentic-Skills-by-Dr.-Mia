from typing import Dict, Any, Optional
import json

class SQLGenerationAgent:
    """
    An agentic skill for generating SQL queries from natural language descriptions.
    In a full deployment, this would interface with an LLM (like Claude or GPT-4).
    This implementation provides a mock/template framework for the agentic workflow.
    """
    
    def __init__(self, model_name: str = "local-template"):
        self.model_name = model_name

    def generate_sql(self, natural_language_query: str, schema_context: Optional[str] = None) -> Dict[str, Any]:
        """
        Translates a natural language query into a SQL statement.
        
        Args:
            natural_language_query (str): The user's request (e.g., 'Get all users who signed up today').
            schema_context (str, optional): DDL or schema description to ground the query.
            
        Returns:
            Dict containing the generated SQL and execution metadata.
        """
        
        # In a real scenario, we would use self.adapter.generate(prompt)
        prompt = self._build_prompt(natural_language_query, schema_context)
        
        # Simulating LLM generation with a heuristic approach for demonstration
        simulated_sql = self._heuristic_sql_generator(natural_language_query)
        
        return {
            "status": "success",
            "model_used": self.model_name,
            "original_query": natural_language_query,
            "generated_sql": simulated_sql,
            "prompt_template": prompt
        }

    def _build_prompt(self, query: str, schema: Optional[str]) -> str:
        prompt = "You are an expert Data Engineer. Write a valid SQL query for the following request.\n"
        if schema:
            prompt += f"\nDatabase Schema Context:\n{schema}\n"
        prompt += f"\nUser Request: {query}\n"
        prompt += "\nOutput ONLY the SQL code."
        return prompt
        
    def _heuristic_sql_generator(self, query: str) -> str:
        """A simple mock generator for demonstration purposes."""
        lower_q = query.lower()
        if "user" in lower_q and "today" in lower_q:
            return "SELECT * FROM users WHERE DATE(created_at) = CURRENT_DATE;"
        elif "count" in lower_q or "how many" in lower_q:
            return "SELECT COUNT(*) FROM table_name;"
        else:
            return "SELECT * FROM my_table LIMIT 10; -- Mocked SQL"

if __name__ == "__main__":
    agent = SQLGenerationAgent()
    
    print("Test 1: Simple query")
    res1 = agent.generate_sql("Get all users who signed up today")
    print(json.dumps(res1, indent=2))
    
    print("\nTest 2: With schema")
    schema = "CREATE TABLE employees (id INT, name VARCHAR, salary DECIMAL);"
    res2 = agent.generate_sql("Find the top 5 highest paid employees", schema_context=schema)
    print(json.dumps(res2, indent=2))
