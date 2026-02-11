# Hugging Face Transformers — Usage in This Project

## Confirmed: We Use Hugging Face

This project depends on the [Hugging Face `transformers`](https://huggingface.co/docs/transformers)
library (`^4.40`, declared in `pyproject.toml`).  It is the backbone of our
**protein embedding pipeline**.

### Where It Appears

| File | What It Does |
|------|--------------|
| `src/features/protein_embeddings.py` | Loads the **ESM-2** protein language model via `AutoTokenizer` + `AutoModel` |
| `pyproject.toml` | Declares `transformers = "^4.40"` and `torch = "^2.2"` |
| `poetry.lock` | Pins `huggingface-hub` and `tokenizers` as transitive dependencies |

### How It's Used (in 4 lines)

```python
from transformers import AutoTokenizer, AutoModel

tokenizer = AutoTokenizer.from_pretrained("facebook/esm2_t6_8M_UR50D")
model     = AutoModel.from_pretrained("facebook/esm2_t6_8M_UR50D")
# → tokenize protein sequences → forward pass → mean-pool → 320-d vector
```

The model name is **configurable** (`settings.esm2_model_name`), so swapping
ESM-2 for ProtTrans or Ankh requires a single config change — no code changes.

---

## What Are Transformers?

The **Transformer** is the neural-network architecture behind modern language
models (GPT, BERT, ESM).  Its core innovation is **self-attention**: every
position in a sequence can attend to every other position in a single step,
rather than processing tokens left-to-right.

### Architecture at a Glance

```
Input Tokens     →  Embedding Layer
                        ↓
              ┌─────────────────────┐
              │  Self-Attention     │ ← each token "looks at" all others
              │  + Feed-Forward     │
              │  + Layer Norm       │
              └────────┬────────────┘
                       ↓   × N layers
              Per-Token Hidden States
                       ↓
              Pooling / Task Head
```

**Key properties:**

| Property | Why It Matters |
|----------|---------------|
| **Parallelizable** | Unlike RNNs, all positions are processed simultaneously |
| **Long-range context** | Attention spans the full sequence in one step |
| **Pre-trainable** | Train once on billions of sequences, fine-tune or extract features cheaply |
| **Architecture-agnostic** | Works on text, protein sequences, molecular graphs, images |

---

## Why Transformers for Proteins?

Protein sequences are *biological language* — strings of amino-acid "letters"
(A, C, D, E, …) whose order determines 3D structure and function.

**Protein language models** like ESM-2 are trained on >250 million sequences via
**masked-language modelling** (predict hidden residues from context), learning:

- **Local motifs** — binding sites, catalytic triads
- **Long-range contacts** — residues far apart in sequence but close in 3D space
- **Evolutionary conservation** — functional constraints across species

The resulting per-residue embeddings capture **structural and functional
information** far richer than hand-crafted features.

### How We Use Them

```
Raw sequence       → Tokenizer  → Token IDs
Token IDs          → ESM-2      → Per-residue embeddings  (L × 320)
Per-residue embeddings → Mean-pool → Single 320-d vector
320-d vector       → XGBoost feature row (emb_0 … emb_319)
```

This 320-dimensional vector is the **protein half** of the feature vector that
XGBoost uses for interaction prediction.

---

## Why Hugging Face Specifically?

| Alternative | Why We Chose HuggingFace Instead |
|-------------|----------------------------------|
| Facebook `fair-esm` | Locks you into ESM models only |
| PyTorch Hub | No unified tokenizer API; manual padding/truncation |
| Custom download | Reinventing model management, caching, and serialization |

Hugging Face provides a **model-agnostic API** (`AutoTokenizer` / `AutoModel`)
that handles weight downloads, caching, tokenization, and GPU placement.
Swapping the underlying model is a one-line config change:

```bash
# .env
ESM2_MODEL_NAME=facebook/esm2_t6_8M_UR50D    # 8M params, fast
# ESM2_MODEL_NAME=facebook/esm2_t33_650M_UR50D  # 650M params, more accurate
# ESM2_MODEL_NAME=Rostlab/prot_t5_xl_uniref50   # ProtTrans alternative
```

> **Summary**: Hugging Face `transformers` is a first-class dependency in this
> project. It provides the protein-embedding backbone (ESM-2) that converts
> amino-acid sequences into the dense numeric vectors XGBoost consumes for
> interaction prediction.
