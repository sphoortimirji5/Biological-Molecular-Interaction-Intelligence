import uuid
from datetime import datetime
from typing import Optional, Dict, Any
from sqlalchemy import String, Text, TIMESTAMP, ForeignKey, Float, Integer, UniqueConstraint, Index
from sqlalchemy.dialects.postgresql import UUID, JSONB
from sqlalchemy.orm import Mapped, mapped_column, relationship
from sqlalchemy.sql import text, func
from src.db.database import Base

class IngestionRun(Base):
    __tablename__ = "ingestion_runs"

    run_id: Mapped[uuid.UUID] = mapped_column(UUID(as_uuid=True), primary_key=True, server_default=text("gen_random_uuid()"))
    status: Mapped[str] = mapped_column(String, nullable=False)  # STARTED, COMPLETED, FAILED
    source: Mapped[str] = mapped_column(String, nullable=False) # e.g. "chembl", "uniprot"
    stats: Mapped[Optional[Dict[str, Any]]] = mapped_column(JSONB, nullable=True)
    checksums: Mapped[Optional[Dict[str, Any]]] = mapped_column(JSONB, nullable=True)
    started_at: Mapped[datetime] = mapped_column(TIMESTAMP(timezone=True), server_default=func.now(), default=datetime.utcnow)
    ended_at: Mapped[Optional[datetime]] = mapped_column(TIMESTAMP(timezone=True), nullable=True)

class Protein(Base):
    __tablename__ = "proteins"

    id: Mapped[uuid.UUID] = mapped_column(UUID(as_uuid=True), primary_key=True, server_default=text("gen_random_uuid()"))
    source: Mapped[str] = mapped_column(String, nullable=False)
    external_id: Mapped[str] = mapped_column(String, nullable=False)
    sequence: Mapped[str] = mapped_column(Text, nullable=False) # Changed to Text
    metadata_: Mapped[Optional[Dict[str, Any]]] = mapped_column("metadata", JSONB, nullable=True)

    __table_args__ = (
        UniqueConstraint('source', 'external_id', name='uq_protein_source_external_id'),
    )

class Compound(Base):
    __tablename__ = "compounds"

    id: Mapped[uuid.UUID] = mapped_column(UUID(as_uuid=True), primary_key=True, server_default=text("gen_random_uuid()"))
    source: Mapped[str] = mapped_column(String, nullable=False)
    external_id: Mapped[str] = mapped_column(String, nullable=False)
    smiles: Mapped[str] = mapped_column(String, nullable=False)
    metadata_: Mapped[Optional[Dict[str, Any]]] = mapped_column("metadata", JSONB, nullable=True)

    __table_args__ = (
        UniqueConstraint('source', 'external_id', name='uq_compound_source_external_id'),
    )

class Interaction(Base):
    __tablename__ = "interactions"

    id: Mapped[uuid.UUID] = mapped_column(UUID(as_uuid=True), primary_key=True, server_default=text("gen_random_uuid()"))
    compound_id: Mapped[uuid.UUID] = mapped_column(ForeignKey("compounds.id"), index=True) # Added Index
    protein_id: Mapped[uuid.UUID] = mapped_column(ForeignKey("proteins.id"), index=True) # Added Index
    
    interaction_type: Mapped[str] = mapped_column(String, default="binds")
    label: Mapped[int] = mapped_column(Integer, default=1)  # 0/1 classification
    
    affinity_value: Mapped[Optional[float]] = mapped_column(Float, nullable=True)
    affinity_unit: Mapped[Optional[str]] = mapped_column(String, nullable=True)
    assay_type: Mapped[Optional[str]] = mapped_column(String, nullable=True)
    confidence: Mapped[Optional[float]] = mapped_column(Float, nullable=True)
    
    source: Mapped[str] = mapped_column(String, nullable=False)
    external_id: Mapped[str] = mapped_column(String, nullable=False) 
    ingested_at: Mapped[datetime] = mapped_column(TIMESTAMP(timezone=True), server_default=func.now(), default=datetime.utcnow)

    # Relationships
    compound = relationship("Compound")
    protein = relationship("Protein")

    __table_args__ = (
        UniqueConstraint('source', 'external_id', name='uq_interaction_source_external_id'),
    )
