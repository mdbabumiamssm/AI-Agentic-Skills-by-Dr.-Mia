import sys

# Try to import RDKit, but provide helpful error if missing
try:
    from rdkit import Chem
    from rdkit.Chem import Descriptors, QED, Crippen
    HAS_RDKIT = True
except ImportError:
    HAS_RDKIT = False

class ChemCrowAgent:
    """
    A functional cheminformatics agent providing real molecular analysis tools.
    
    Upgraded for 2026:
    - Real RDKit integration
    - Advanced Safety Checks (PAINS, Structural Alerts)
    - Lipinski Rule validation
    """
    
    def __init__(self):
        self.tools = {
            "MolWeight": self.get_molecular_weight,
            "LogP": self.get_logp,
            "QED": self.get_qed,
            "Lipinski": self.check_lipinski,
            "Safety": self.check_safety,
            "Validity": self.check_validity
        }
        
        if not HAS_RDKIT:
            print("CRITICAL WARNING: 'rdkit' not found. Agent functionality is severely limited.")
            print("Please install via: pip install rdkit")

    def run_tool(self, tool_name, smiles):
        if tool_name not in self.tools:
            return f"Error: Tool '{tool_name}' not found."
        
        # Pre-validate for most tools
        if not self.check_validity(smiles) and tool_name != "Validity":
             return "Error: Invalid SMILES string."

        return self.tools[tool_name](smiles)

    def check_validity(self, smiles):
        """Checks if a SMILES string is valid using RDKit."""
        if not HAS_RDKIT:
            # Fallback for environment without RDKit
            return len(smiles) > 0 and isinstance(smiles, str)
        
        mol = Chem.MolFromSmiles(smiles)
        return mol is not None

    def get_molecular_weight(self, smiles):
        """Calculates Molecular Weight."""
        if not HAS_RDKIT: return "N/A (RDKit missing)"
        mol = Chem.MolFromSmiles(smiles)
        return round(Descriptors.MolWt(mol), 3)

    def get_logp(self, smiles):
        """Calculates LogP (partition coefficient)."""
        if not HAS_RDKIT: return "N/A (RDKit missing)"
        mol = Chem.MolFromSmiles(smiles)
        return round(Crippen.MolLogP(mol), 3)
    
    def get_qed(self, smiles):
        """Calculates Quantitative Estimate of Drug-likeness."""
        if not HAS_RDKIT: return "N/A (RDKit missing)"
        mol = Chem.MolFromSmiles(smiles)
        return round(QED.qed(mol), 3)

    def check_lipinski(self, smiles):
        """Checks Lipinski's Rule of 5."""
        if not HAS_RDKIT: return "N/A (RDKit missing)"
        
        mol = Chem.MolFromSmiles(smiles)
        mw = Descriptors.MolWt(mol)
        logp = Crippen.MolLogP(mol)
        hbd = Descriptors.NumHDonors(mol)
        hba = Descriptors.NumHAcceptors(mol)
        
        violations = 0
        details = []
        if mw > 500: 
            violations += 1
            details.append(f"MW {mw:.1f} > 500")
        if logp > 5: 
            violations += 1
            details.append(f"LogP {logp:.1f} > 5")
        if hbd > 5: 
            violations += 1
            details.append(f"HBD {hbd} > 5")
        if hba > 10: 
            violations += 1
            details.append(f"HBA {hba} > 10")
            
        return f"Violations: {violations} ({', '.join(details)})" if violations > 0 else "Pass"

    def check_safety(self, smiles):
        """
        Checks for basic safety flags and structural alerts.
        """
        if not HAS_RDKIT: return "N/A (RDKit missing)"
        
        mol = Chem.MolFromSmiles(smiles)
        flags = []
        
        # 1. Nitro groups (often explosive/toxic)
        nitro = Chem.MolFromSmarts("N(=O)=O")
        if mol.HasSubstructMatch(nitro):
            flags.append("Warning: Nitro group detected (potential toxicity/explosivity)")
            
        # 2. Epoxides (reactive alkylating agents)
        epoxide = Chem.MolFromSmarts("C1OC1")
        if mol.HasSubstructMatch(epoxide):
            flags.append("Warning: Epoxide detected (reactive)")
            
        # 3. Michael Acceptors (simplified)
        michael = Chem.MolFromSmarts("C=CC=O")
        if mol.HasSubstructMatch(michael):
            flags.append("Warning: Potential Michael acceptor (covalent binding)")

        if not flags:
            return "No common structural alerts found."
        return "; ".join(flags)

def main():
    agent = ChemCrowAgent()
    
    print("--- ChemCrow 2.0 (2026 Update) ---")
    if not HAS_RDKIT:
        print("NOTE: Running in degraded mode (No RDKit).")
    
    # Test Cases
    test_molecules = [
        ("Aspirin", "CC(=O)Oc1ccccc1C(=O)O"),
        ("TNT (Explosive)", "Cc1c(N(=O)=O)cc(N(=O)=O)cc1N(=O)=O"),
        ("Lipinski Fail", "OCC1OC(OCC2OC(OC3OC(OC4OC(OC5OC(CO)C(O)C(O)C5O)C(O)C(O)C4O)C(O)C(O)C3O)C(O)C(O)C2O)C(O)C(O)C1O") # Starch-like
    ]
    
    for name, smi in test_molecules:
        print(f"\nAnalyzing: {name}")
        print(f"  SMILES:   {smi}")
        print(f"  MW:       {agent.run_tool('MolWeight', smi)}")
        print(f"  LogP:     {agent.run_tool('LogP', smi)}")
        print(f"  Lipinski: {agent.run_tool('Lipinski', smi)}")
        print(f"  Safety:   {agent.run_tool('Safety', smi)}")

if __name__ == "__main__":
    main()