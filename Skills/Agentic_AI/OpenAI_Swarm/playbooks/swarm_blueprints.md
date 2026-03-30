# Swarm Blueprint Library (Stub)

| Blueprint | Agents | Notes |
|---|---|---|
| Hematology Crash Cart | IntakeAgent → DiagnosticsAgent → SafetyOfficer | TODO: Fill in labs + drug-SOP matrix.
| Oncology Prior Auth | Navigator → EvidenceSummarizer → AppealWriter | TODO: Map CPT/ICD from `Clinical/Prior_Authorization`.
| Self-Driving Lab Loop | HypothesisAgent → RobotProxy → DataReviewer | TODO: integrate with `Self_Driving_Labs` adapters.

> Populate this document with JSON/YAML agent specs derived from validated workflows.
