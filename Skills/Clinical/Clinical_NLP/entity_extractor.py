"""
entity_extractor.py

A lightweight rule-based Clinical NLP extractor.
Identifies common medical terms using regex patterns.
"""

import argparse
import json
import re
import sys

class ClinicalNLP:
    def __init__(self):
        # Very simple regex patterns for demonstration
        self.patterns = {
            "PROBLEM": [
                r"(?i)diabetes\s*(type\s*\d)?",
                r"(?i)hypertension",
                r"(?i)pneumonia",
                r"(?i)chest pain",
                r"(?i)fracture",
                r"(?i)cancer"
            ],
            "MEDICATION": [
                r"(?i)metformin",
                r"(?i)lisinopril",
                r"(?i)aspirin",
                r"(?i)insulin",
                r"(?i)atorvastatin"
            ],
            "PROCEDURE": [
                r"(?i)x-ray",
                r"(?i)ct scan",
                r"(?i)mri",
                r"(?i)biopsy",
                r"(?i)surgery"
            ]
        }
        
        self.negation_triggers = ["no ", "denies ", "without ", "negative for "]

    def extract(self, text: str):
        entities = []
        
        for entity_type, regex_list in self.patterns.items():
            for pattern in regex_list:
                for match in re.finditer(pattern, text):
                    term = match.group()
                    start, end = match.span()
                    
                    # Check for negation (lookbehind window)
                    context_start = max(0, start - 20)
                    context = text[context_start:start].lower()
                    is_negated = any(trigger in context for trigger in self.negation_triggers)
                    
                    entities.append({
                        "text": term,
                        "type": entity_type,
                        "start": start,
                        "end": end,
                        "negated": is_negated
                    })
                    
        return sorted(entities, key=lambda x: x['start'])

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Clinical NLP Entity Extractor")
    parser.add_argument("--text", help="Clinical text to analyze")
    parser.add_argument("--file", help="Path to text file")
    parser.add_argument("--output", help="Path to save JSON output")
    
    args = parser.parse_args()
    
    input_text = ""
    if args.text:
        input_text = args.text
    elif args.file:
        try:
            with open(args.file, 'r') as f:
                input_text = f.read()
        except Exception as e:
            print(f"Error reading file: {e}", file=sys.stderr)
            sys.exit(1)
    else:
        print("Error: Must provide --text or --file", file=sys.stderr)
        sys.exit(1)
        
    nlp = ClinicalNLP()
    results = nlp.extract(input_text)
    
    print(json.dumps(results, indent=2))
    
    if args.output:
        with open(args.output, 'w') as f:
            json.dump(results, f, indent=2)
