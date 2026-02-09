# Security

This document defines the security posture of **Biological Molecular Interaction Intelligence**.

The platform uses **publicly available biological and molecular databases** for this project. However, the **same security controls are applied uniformly**, as if all data were proprietary research data. This ensures safe defaults, consistency, and future extensibility without changing the security model.

---

## 1. Design Assumptions

* Public datasets (e.g., UniProt, ChEMBL) are used for the current implementation.
* All molecular, protein, and interaction data is treated as **confidential research data by default**, regardless of source.
* The platform explicitly **does not support patient-linked, individual-derived, or clinical data**.

This design avoids special-case handling for public vs proprietary datasets and enforces a single, high-assurance security posture.

---

## 2. Data Classification

### 2.1 Molecular & Protein Research Data

Includes:

* Protein sequences
* Drug molecules (SMILES, SDF)
* Drug–protein and protein–protein interaction records
* Derived features, embeddings, and model artifacts

**Classification**: Confidential Research Data (by policy)

**PII / PHI**: ❌ No

**Regulatory Scope**:

* Not HIPAA
* Not PHI
* Not PII

Security requirements are driven by **data integrity, access control, and intellectual property protection**, not healthcare regulation.

---

## 3. Core Security Principles

### 3.1 Uniform Security Controls

* Identical controls for public and private datasets
* No security relaxation based on data source
* Predictable behavior across environments

---

### 3.2 Access Control

* No anonymous access to any system component

* Environment-scoped credentials for all services

* Role-based access separation:

  * ingestion
  * feature generation
  * model training
  * inference

* Least-privilege permissions enforced

---

### 3.3 Encryption

* Encryption in transit (TLS) for:

  * API traffic
  * database connections
  * object storage access

* Encryption at rest for:

  * PostgreSQL data volumes
  * Object storage (MinIO)

---

### 3.4 Auditability & Lineage

* All ingestion runs recorded with timestamps and checksums
* Dataset versions are immutable once ingested
* Object storage paths are versioned and content-addressed
* Model artifacts are stored with:

  * run identifiers
  * feature manifests
  * model signatures

---

### 3.5 Logging Policy

* No raw protein sequences, SMILES, or embeddings in logs

* Logs restricted to:

  * identifiers
  * record counts
  * execution metadata

* Debug logging disabled by default

---

## 4. Environment Isolation

* Development, test, and production environments are isolated
* Credentials are scoped per environment
* No cross-environment data reuse

---

## 5. Explicit Non-Goals

Biological Molecular Interaction Intelligence does **not**:

* Handle PHI or PII
* Accept patient identifiers
* Process individual-derived genomic or proteomic samples
* Claim healthcare regulatory compliance

Any system requiring such data must be handled separately.

---

## 6. Rationale

Using uniform security controls across public and proprietary datasets:

* Simplifies system design
* Prevents accidental data exposure
* Eliminates policy drift as data sources evolve
* Aligns with standard practices in molecular research platforms
