# XGBoost — A Numeric → Numeric Model

## The Fundamental Constraint

XGBoost is a **gradient-boosted decision-tree** ensemble.  Every split in every tree
evaluates a **numeric threshold** — `feature_j < value`.  This means:

| What XGBoost *accepts*       | What XGBoost *produces*      |
|------------------------------|------------------------------|
| A matrix of **numbers** (float/int) | A vector of **numbers** (probability or regression target) |

It cannot natively consume text, sequences, graphs, or molecular structures.

---

## Why This Matters for Drug–Target Interaction

Our pipeline predicts whether a small-molecule **drug** interacts with a **protein
target**.  The raw inputs are:

| Entity    | Raw Form                     | Example                    |
|-----------|------------------------------|----------------------------|
| Compound  | SMILES string                | `CC(=O)Oc1ccccc1C(=O)O`   |
| Protein   | Amino-acid sequence (FASTA)  | `MKWVTFISLLLLFSSA…`        |

Neither of these is numeric.

---

## How We Bridge the Gap

The project's **Phase 2 — Molecular Featurization** converts every entity into a
fixed-length numeric vector *before* XGBoost ever sees the data.

```
Compound  ──[RDKit Morgan FP]──►  2 048-d binary vector  (fp_0 … fp_2047)
Protein   ──[ESM-2 Embeddings]──►    320-d dense vector   (emb_0 … emb_319)
```

At training time these vectors are **concatenated** into a single feature row
(2 368 columns) paired with a binary label:

```
┌─────────────┬──────────────────┬───────┐
│ fp_0…fp_2047│ emb_0…emb_319    │ label │
├─────────────┼──────────────────┼───────┤
│ 0, 1, 0, …  │ -0.12, 0.44, …  │   1   │   ← known interaction
│ 1, 0, 1, …  │  0.31, -0.08, … │   0   │   ← no interaction
└─────────────┴──────────────────┴───────┘
```

XGBoost then does what it does best: learns non-linear numeric splits over these
2 368 features to output a **probability** of interaction.

---

## In This Codebase

| Step | Module | Output |
|------|--------|--------|
| Compound featurization | `src/features/fingerprint_generator.py` | 2 048-d Morgan fingerprint |
| Protein featurization  | `src/features/protein_embeddings.py`    | 320-d ESM-2 embedding |
| Feature concatenation  | `src/models/training_data.py`           | Combined DataFrame |
| Model training         | `src/models/xgboost_model.py`           | `XGBClassifier` (binary:logistic) |
| Inference              | `src/models/scoring.py` (`LearnedScorer`) | Ranked probabilities |

> **Key takeaway**: XGBoost never touches raw biology. The featurization layer
> translates molecular structure into numbers; XGBoost translates numbers into
> predictions. Separating these concerns is what makes the scoring layer
> pluggable — swap XGBoost for a GNN without touching featurization, or swap
> ESM-2 for ProtTrans without touching the model.
