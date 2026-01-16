"""
acmg_classifier.py

A rule-based classifier for ACMG Variant Interpretation.
Logic based on Richards et al. 2015.
"""

import argparse
import json
import sys

class ACMGClassifier:
    def __init__(self):
        # Weights for Pathogenicity
        self.p_weights = {
            "PVS1": 8,  # Very Strong
            "PS1": 4, "PS2": 4, "PS3": 4, "PS4": 4, # Strong
            "PM1": 2, "PM2": 2, "PM3": 2, "PM4": 2, "PM5": 2, "PM6": 2, # Moderate
            "PP1": 1, "PP2": 1, "PP3": 1, "PP4": 1, "PP5": 1 # Supporting
        }
        
        # Weights for Benign
        self.b_weights = {
            "BA1": 8, # Stand-alone
            "BS1": 4, "BS2": 4, "BS3": 4, "BS4": 4, # Strong
            "BP1": 1, "BP2": 1, "BP3": 1, "BP4": 1, "BP5": 1, "BP6": 1 # Supporting
        }

    def classify(self, evidence_codes: list):
        p_score = 0
        b_score = 0
        
        applied_codes = []
        
        for code in evidence_codes:
            code = code.upper().strip()
            if code in self.p_weights:
                p_score += self.p_weights[code]
                applied_codes.append({"code": code, "type": "Pathogenic", "weight": self.p_weights[code]})
            elif code in self.b_weights:
                b_score += self.b_weights[code]
                applied_codes.append({"code": code, "type": "Benign", "weight": self.b_weights[code]})
            else:
                print(f"Warning: Unknown code {code}", file=sys.stderr)

        # Simplified Scoring Logic (Richards et al is more complex logic-wise, this is a score-based approximation)
        # Pathogenic: PVS1 + (1 Strong OR 2 Mod etc...)
        # Here we use a sum threshold for demonstration
        
        verdict = "Uncertain Significance (VUS)"
        
        # Pathogenic logic (Approximation)
        if "PVS1" in evidence_codes:
            if p_score >= 10: verdict = "Pathogenic"
            elif p_score >= 9: verdict = "Likely Pathogenic"
        elif p_score >= 12: # e.g. 3 Strong
            verdict = "Pathogenic"
        elif p_score >= 6:
            verdict = "Likely Pathogenic"
            
        # Benign logic
        if "BA1" in evidence_codes:
            verdict = "Benign"
        elif b_score >= 8: # 2 Strong
            verdict = "Benign"
        elif b_score >= 4:
            verdict = "Likely Benign"
            
        return {
            "verdict": verdict,
            "pathogenicity_score": p_score,
            "benign_score": b_score,
            "evidence_applied": applied_codes
        }

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="ACMG Variant Classifier")
    parser.add_argument("--evidence", required=True, help="Comma-separated list of codes (e.g., PVS1,PM2)")
    parser.add_argument("--output", help="Path to save JSON output")
    
    args = parser.parse_args()
    
    codes = args.evidence.split(",")
    classifier = ACMGClassifier()
    result = classifier.classify(codes)
    
    print(json.dumps(result, indent=2))
    
    if args.output:
        with open(args.output, 'w') as f:
            json.dump(result, f, indent=2)
