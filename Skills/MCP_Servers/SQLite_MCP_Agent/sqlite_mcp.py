import sqlite3
import json
from typing import Dict, Any

class SQLiteMCPAgent:
    """
    Model Context Protocol (MCP) compatible Agent for local SQLite databases.
    Allows LLMs to securely query, introspect, and analyze local SQLite files.
    """
    
    def __init__(self, db_path: str):
        self.db_path = db_path
        
    def execute_query(self, query: str) -> Dict[str, Any]:
        """
        Executes a read-only SQL query against the SQLite database.
        """
        if not query.strip().upper().startswith("SELECT"):
            return {"error": "Only SELECT queries are allowed for safety."}
            
        try:
            with sqlite3.connect(self.db_path) as conn:
                cursor = conn.cursor()
                cursor.execute(query)
                columns = [desc[0] for desc in cursor.description]
                rows = cursor.fetchall()
                
                # Format as a list of dicts for the LLM
                results = [dict(zip(columns, row)) for row in rows]
                
                return {
                    "status": "success",
                    "row_count": len(results),
                    "data": results
                }
        except Exception as e:
            return {"status": "error", "message": str(e)}

    def get_schema(self) -> Dict[str, Any]:
        """
        Introspects the database to return its schema (tables and columns).
        Useful for providing context to the LLM before it writes queries.
        """
        schema_query = "SELECT name FROM sqlite_master WHERE type='table';"
        try:
            with sqlite3.connect(self.db_path) as conn:
                cursor = conn.cursor()
                cursor.execute(schema_query)
                tables = cursor.fetchall()
                
                db_schema = {}
                for table in tables:
                    table_name = table[0]
                    cursor.execute(f"PRAGMA table_info({table_name});")
                    columns = cursor.fetchall()
                    # column format: (cid, name, type, notnull, dflt_value, pk)
                    db_schema[table_name] = [{"name": col[1], "type": col[2]} for col in columns]
                    
                return {
                    "status": "success",
                    "schema": db_schema
                }
        except Exception as e:
            return {"status": "error", "message": str(e)}

if __name__ == "__main__":
    # Quick self-test with an in-memory DB
    print("Testing SQLiteMCPAgent...")
    
    # Setup dummy db
    conn = sqlite3.connect("dummy.db")
    conn.execute("CREATE TABLE IF NOT EXISTS users (id INTEGER PRIMARY KEY, name TEXT, role TEXT);")
    conn.execute("INSERT OR IGNORE INTO users (id, name, role) VALUES (1, 'Alice', 'Admin'), (2, 'Bob', 'User');")
    conn.commit()
    conn.close()
    
    agent = SQLiteMCPAgent("dummy.db")
    
    print("\n--- Schema Introspection ---")
    print(json.dumps(agent.get_schema(), indent=2))
    
    print("\n--- Query Execution ---")
    print(json.dumps(agent.execute_query("SELECT * FROM users WHERE role = 'Admin';"), indent=2))
