<!--
# COPYRIGHT NOTICE
# This file is part of the "Universal AI Agentic Skills" project.
# Copyright (c) 2026 MD BABU MIA, PhD <md.babu.mia@mssm.edu>
# All Rights Reserved.
#
# This code is proprietary and confidential.
# Unauthorized copying of this file, via any medium is strictly prohibited.
#
# Provenance: Authenticated by MD BABU MIA

-->

# Product Review Agent

A complex structured output agent demonstrating advanced Pydantic schemas for product review analysis.

## 🎯 What This Demonstrates

- **Complex Schemas**: Multi-field Pydantic models with various data types
- **List Fields**: Arrays of strings for pros/cons analysis
- **Boolean Logic**: Recommendation decisions based on review content
- **Sentiment Analysis**: Automated sentiment classification

## 🚀 Quick Start

1. **Install OpenAI Agents SDK**:
   ```bash
   pip install openai-agents
   ```

2. **Set up environment**:
   ```bash
   cp ../env.example .env
   # Edit .env and add your OpenAI API key
   ```

3. **Run the agent**:
   ```python
   from agents import Runner
   from agent import root_agent
   
   review_text = "This laptop is amazing! Great performance, long battery life, but a bit heavy."
   result = Runner.run_sync(root_agent, f"Analyze this review: {review_text}")
   print(result.final_output)  # Returns ProductReview object
   ```

## 💡 Key Concepts

- **Rating Validation**: Integer constraints (1-5 stars)
- **Sentiment Enum**: Automatic positive/negative/neutral classification
- **List Processing**: Extracting multiple pros and cons
- **Optional Fields**: Handling missing reviewer information

## 🧪 Example Output

```python
{
    "product_name": "Gaming Laptop XYZ",
    "rating": 4,
    "summary": "Great performance but heavy design",
    "sentiment": "positive",
    "pros": ["Great performance", "Long battery life", "Good display"],
    "cons": ["Heavy weight", "Expensive price"],
    "recommend": true,
    "reviewer_name": "TechEnthusiast123"
}
```

## 🔗 Next Steps

- [Support Ticket Agent](../2_1_support_ticket_agent/README.md) - Basic structured output
- [Tutorial 3: Tool Using Agent](../../3_tool_using_agent/README.md) - Adding tools to agents


<!-- AUTHOR_SIGNATURE: 9a7f3c2e-MD-BABU-MIA-2026-MSSM-SECURE -->