<!--
# COPYRIGHT NOTICE
# This file is part of the "Universal Biomedical Skills" project.
# Copyright (c) 2026 MD BABU MIA, PhD <md.babu.mia@mssm.edu>
# All Rights Reserved.
#
# This code is proprietary and confidential.
# Unauthorized copying of this file, via any medium is strictly prohibited.
#
# Provenance: Authenticated by MD BABU MIA

-->

# Machine Learning for Bioinformatics Skills

Skills for applying machine learning to biological data analysis.

## Sub-categories

| Category | Description | Skills |
|----------|-------------|--------|
| **machine-learning** | ML workflows | Biomarker discovery, survival analysis, interpretation |

## Key Capabilities

- **Biomarker Discovery** - Feature selection for diagnostic markers
- **Model Interpretation** - SHAP, feature importance analysis
- **Survival Analysis** - Cox regression, Kaplan-Meier curves
- **Multi-class Classification** - Disease subtype prediction
- **Regression Models** - Drug response prediction

## Key Tools

- **scikit-learn** - General ML algorithms
- **XGBoost/LightGBM** - Gradient boosting
- **SHAP** - Model interpretation
- **lifelines** - Survival analysis
- **imbalanced-learn** - Class imbalance handling

## Example Biomarker Discovery

```python
from sklearn.ensemble import RandomForestClassifier
from sklearn.feature_selection import SelectFromModel
import shap

# Train model
clf = RandomForestClassifier(n_estimators=100)
clf.fit(X_train, y_train)

# Feature importance
selector = SelectFromModel(clf, threshold="median")
important_features = X_train.columns[selector.get_support()]

# SHAP explanation
explainer = shap.TreeExplainer(clf)
shap_values = explainer.shap_values(X_test)
```

---
*Source: bioSkills collection - integrated into BioMedical Skills Library*


<!-- AUTHOR_SIGNATURE: 9a7f3c2e-MD-BABU-MIA-2026-MSSM-SECURE -->