# Systems Biology Skills

Skills for systems-level biological modeling and analysis.

## Sub-categories

| Category | Description | Skills |
|----------|-------------|--------|
| **systems-biology** | Metabolic modeling | FBA, metabolic reconstruction, network analysis |

## Key Capabilities

- **Flux Balance Analysis (FBA)** - Metabolic flux prediction
- **Metabolic Reconstruction** - Genome-scale models
- **Network Analysis** - Biological network properties
- **Dynamic Modeling** - ODE-based simulations
- **Constraint-Based Modeling** - Growth predictions

## Key Tools

- **COBRApy** - Constraint-based modeling
- **cobrapy** - Python COBRA implementation
- **CarveMe** - Automated model reconstruction
- **MEMOTE** - Model quality assessment
- **Escher** - Pathway visualization

## Example FBA

```python
import cobra
from cobra.io import read_sbml_model

# Load genome-scale model
model = read_sbml_model("e_coli_core.xml")

# Set objective (maximize growth)
model.objective = "BIOMASS_Ecoli_core_w_GAM"

# Run FBA
solution = model.optimize()

print(f"Growth rate: {solution.objective_value:.3f}")
print(f"Status: {solution.status}")

# Get flux distribution
fluxes = solution.fluxes
active_reactions = fluxes[abs(fluxes) > 1e-6]
```

## Gene Knockout Analysis

```python
# Single gene knockout
with model:
    model.genes.get_by_id("b0720").knock_out()
    ko_solution = model.optimize()
    print(f"Growth after KO: {ko_solution.objective_value:.3f}")

# Essential genes screen
from cobra.flux_analysis import single_gene_deletion

results = single_gene_deletion(model)
essential = results[results['growth'] < 0.01]
```

---
*Source: bioSkills collection - integrated into BioMedical Skills Library*
