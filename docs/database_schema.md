# Database Schema

## Entity Relationship Diagram

```mermaid
erDiagram
    ingestion_runs {
        uuid run_id PK
        string source
        string status
        jsonb stats
        jsonb checksums
        timestamp started_at
        timestamp ended_at
    }

    proteins {
        uuid id PK
        string source
        string external_id
        string sequence
        jsonb metadata
    }

    compounds {
        uuid id PK
        string source
        string external_id
        string smiles
        jsonb metadata
    }

    interactions {
        uuid id PK
        uuid compound_id FK
        uuid protein_id FK
        string interaction_type
        int label
        float affinity_value
        string affinity_unit
        string assay_type
        float confidence
        string source
        string external_id
        timestamp ingested_at
    }

    interactions }|--|| compounds : "references"
    interactions }|--|| proteins : "references"
```

## Constraints

*   **Proteins**: `UNIQUE(source, external_id)`
*   **Compounds**: `UNIQUE(source, external_id)`
*   **Interactions**: `UNIQUE(source, external_id)`
    *   *Note: `external_id` refers to the source-provided record identifier.*
