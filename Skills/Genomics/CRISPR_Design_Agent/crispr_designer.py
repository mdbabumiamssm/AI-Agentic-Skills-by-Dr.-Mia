"""
CRISPR Design Agent (2026 Skills)

A specialized bioinformatics tool for designing gRNA (guide RNA) sequences.
Functionality:
1. Scan DNA for PAM sites (NGG).
2. Extract 20bp spacers (Protospacers).
3. Score efficiency (GC Content & Homopolymer checks).
4. Check for off-targets (Mocked).
"""

import re

class CRISPRDesigner:
    def __init__(self):
        self.pam_motif = "GG" # SpCas9 PAM is NGG
        self.spacer_length = 20

    def find_targets(self, dna_sequence: str):
        """
        Scans sequence for 'GG' and looks 21-23bp upstream.
        """
        dna = dna_sequence.upper()
        candidates = []
        
        # Iterate through sequence to find PAMs
        # Note: NGG means any base followed by GG.
        # We look for GG, then check index.
        
        for i in range(len(dna) - 1):
            if dna[i:i+2] == "GG":
                pam_index = i
                # NGG starts at pam_index-1. Spacer ends at pam_index-1.
                # Spacer start = (pam_index-1) - 20
                start = pam_index - 1 - self.spacer_length
                
                if start >= 0:
                    spacer = dna[start : start + self.spacer_length]
                    pam_full = dna[start + self.spacer_length : start + self.spacer_length + 3]
                    
                    score = self.calculate_score(spacer)
                    
                    candidates.append({
                        "location": start,
                        "spacer": spacer,
                        "pam": pam_full,
                        "gc_content": score["gc"],
                        "efficiency_score": score["efficiency"]
                    })
                    
        return sorted(candidates, key=lambda x: x["efficiency_score"], reverse=True)

    def calculate_score(self, spacer: str):
        # 1. GC Content (Ideal: 40-60%)
        g_count = spacer.count('G')
        c_count = spacer.count('C')
        gc_percent = (g_count + c_count) / len(spacer) * 100
        
        score = 0.0
        if 40 <= gc_percent <= 60:
            score += 50
        elif 30 <= gc_percent <= 80:
            score += 25
            
        # 2. Position Specific Weights (Mock Rule: G at pos 20 is good)
        if spacer[-1] == 'G':
            score += 10
            
        # 3. Penalize Poly-T (Termination signal)
        if "TTTT" in spacer:
            score -= 50
            
        return {"gc": round(gc_percent, 1), "efficiency": max(0, score)}

if __name__ == "__main__":
    designer = CRISPRDesigner()
    
    # Example Gene Fragment (TP53 exon mock)
    gene_seq = "ATGGAGGAGCCGCAGTCAGATCCTAGCGTCGAGCCCCCTCTGAGTCAGGAAACATTTTCAGACCTATGGAAACTGTGAGTGGATCCATTGGAAGGGC"
    
    print(f"Scanning DNA ({len(gene_seq)} bp) for CRISPR targets...\n")
    
    targets = designer.find_targets(gene_seq)
    
    print(f"{'Loc':<5} {'Spacer Sequence (20bp)':<22} {'PAM':<4} {'GC%':<5} {'Score'}")
    print("""---------------------"""")
    
    for t in targets:
        print(f"{t['location']:<5} {t['spacer']} {t['pam']} {t['gc_content']:<5} {t['efficiency_score']}")
        
    print(f"\nFound {len(targets)} potential guides.")
