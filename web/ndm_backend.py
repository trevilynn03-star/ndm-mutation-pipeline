from __future__ import annotations

import re
import shutil
import subprocess
import tempfile
from io import StringIO
from pathlib import Path

import pandas as pd
from Bio import AlignIO, SeqIO


BASE = Path(__file__).resolve().parent.parent

REFERENCE_MSA = BASE / "analysis/msa/NDM_MSA.fasta"
PRIMARY_METADATA = BASE / "analysis/reference/primary_sequences.tsv"

REFERENCE_ID = "CBZ39493.1"

EXPECTED_COLUMNS = [
    "accession",
    "variant",
    "reference_position",
    "reference_residue",
    "observed_residue",
    "mutation",
    "msa_column",
]


def get_accession(record) -> str:
    """Return a normalized accession from a FASTA record."""
    record_id = record.id.strip()

    if record_id.startswith("tr|"):
        parts = record_id.split("|")
        if len(parts) >= 2:
            return parts[1]

    return record_id


def extract_variant(description: str) -> str:
    """
    Extract an NDM variant from a sequence description.

    Handles descriptions such as:
    - NDM-1
    - blaNDM-16
    - beta-lactamase NDM-5
    """
    patterns = [
        r"\bbla(NDM-\d+)\b",
        r"\b(NDM-\d+)\b",
    ]

    for pattern in patterns:
        match = re.search(
            pattern,
            description,
            flags=re.IGNORECASE,
        )

        if match:
            value = match.group(1).upper()

            if not value.startswith("NDM-"):
                value = "NDM-" + value

            return value

    return "Unknown"


def load_variant_metadata() -> dict[str, str]:
    """
    Load authoritative accession → NDM variant mapping.

    primary_sequences.tsv contains:
        FASTA description<TAB>sequence_hash

    The accession and variant are therefore extracted from the
    FASTA-style description rather than relying on the MSA header.
    """
    if not PRIMARY_METADATA.exists():
        raise FileNotFoundError(
            f"Primary metadata file not found: {PRIMARY_METADATA}"
        )

    mapping = {}

    with PRIMARY_METADATA.open(
        "r",
        encoding="utf-8",
    ) as handle:

        for line in handle:
            line = line.strip()

            if not line:
                continue

            parts = line.split("\t")

            if len(parts) < 1:
                continue

            description = parts[0]

            first_token = description.split()[0]

            if first_token.startswith("tr|"):
                token_parts = first_token.split("|")

                if len(token_parts) >= 2:
                    accession = token_parts[1]
                else:
                    accession = first_token
            else:
                accession = first_token

            variant = extract_variant(description)

            mapping[accession] = variant

    return mapping


def get_variant_for_record(
    record,
    variant_metadata: dict[str, str],
) -> str:
    """
    Get the variant from authoritative metadata.

    Falls back to the FASTA description only for genuinely new
    uploaded sequences.
    """
    accession = get_accession(record)

    if accession in variant_metadata:
        return variant_metadata[accession]

    return extract_variant(record.description)


def load_locked_alignment():
    """Load the immutable 30-sequence reference MSA."""
    if not REFERENCE_MSA.exists():
        raise FileNotFoundError(
            f"Locked reference MSA not found: {REFERENCE_MSA}"
        )

    return AlignIO.read(
        REFERENCE_MSA,
        "fasta",
    )


def load_reference_record():
    """Return CBZ39493.1 from the locked MSA."""
    alignment = load_locked_alignment()

    for record in alignment:
        if get_accession(record) == REFERENCE_ID:
            return record

    raise ValueError(
        f"Reference {REFERENCE_ID} was not found in the locked MSA."
    )


def build_reference_coordinate_map(reference_record):
    """
    Map MSA columns to biological coordinates of CBZ39493.1.
    """
    coordinate_map = {}

    reference_position = 0

    for msa_column, residue in enumerate(
        str(reference_record.seq),
        start=1,
    ):
        if residue != "-":
            reference_position += 1

            coordinate_map[msa_column] = (
                reference_position,
                residue.upper(),
            )

    return coordinate_map


def run_mafft_add(
    uploaded_fasta: Path,
    output_alignment: Path,
) -> None:
    """
    Add genuinely new sequences to the immutable locked MSA.

    The original locked MSA is never modified.
    """
    mafft = shutil.which("mafft")

    if mafft is None:
        raise RuntimeError(
            "MAFFT was not found on PATH."
        )

    command = [
        mafft,
        "--add",
        str(uploaded_fasta),
        str(REFERENCE_MSA),
    ]

    process = subprocess.run(
        command,
        capture_output=True,
        text=True,
        check=False,
    )

    if process.returncode != 0:
        raise RuntimeError(
            "MAFFT alignment failed:\n"
            + process.stderr
        )

    if not process.stdout.strip():
        raise RuntimeError(
            "MAFFT completed but returned an empty alignment."
        )

    output_alignment.write_text(
        process.stdout,
        encoding="utf-8",
    )


def detect_mutations_for_accessions(
    alignment,
    uploaded_accessions: set[str],
    variant_metadata: dict[str, str],
) -> pd.DataFrame:
    """
    Detect mutations only for uploaded accessions.

    Mutations are measured against CBZ39493.1 using the locked
    reference-coordinate framework.
    """
    reference_record = None

    for record in alignment:
        if get_accession(record) == REFERENCE_ID:
            reference_record = record
            break

    if reference_record is None:
        raise ValueError(
            f"Reference {REFERENCE_ID} was not found in the alignment."
        )

    coordinate_map = build_reference_coordinate_map(
        reference_record
    )

    mutations = []

    for record in alignment:
        accession = get_accession(record)

        if accession not in uploaded_accessions:
            continue

        if accession == REFERENCE_ID:
            continue

        variant = get_variant_for_record(
            record,
            variant_metadata,
        )

        sequence = str(record.seq)

        for msa_column, (
            reference_position,
            reference_residue,
        ) in coordinate_map.items():

            observed_residue = sequence[
                msa_column - 1
            ].upper()

            # Ignore alignment gaps.
            if observed_residue == "-":
                continue

            # Ignore identical residues.
            if observed_residue == reference_residue:
                continue

            mutations.append(
                {
                    "accession": accession,
                    "variant": variant,
                    "reference_position": reference_position,
                    "reference_residue": reference_residue,
                    "observed_residue": observed_residue,
                    "mutation": (
                        f"{reference_residue}"
                        f"{reference_position}"
                        f"{observed_residue}"
                    ),
                    "msa_column": msa_column,
                }
            )

    return pd.DataFrame(
        mutations,
        columns=EXPECTED_COLUMNS,
    )


def analyse_fasta_bytes(fasta_bytes: bytes):
    """
    Analyse uploaded NDM protein FASTA sequences.

    Existing sequences from the locked dataset are taken directly
    from the locked MSA.

    New sequences are aligned using MAFFT --add.

    Only uploaded sequences contribute mutation results.
    """
    variant_metadata = load_variant_metadata()

    with tempfile.TemporaryDirectory(
        prefix="ndm_web_"
    ) as tmpdir:

        tmp = Path(tmpdir)

        uploaded_fasta = tmp / "uploaded.fasta"
        new_sequences_fasta = tmp / "new_sequences.fasta"
        aligned_fasta = tmp / "aligned.fasta"

        uploaded_fasta.write_bytes(fasta_bytes)

        records = list(
            SeqIO.parse(
                str(uploaded_fasta),
                "fasta",
            )
        )

        if not records:
            raise ValueError(
                "No FASTA sequences were found."
            )

        uploaded_accessions = {
            get_accession(record)
            for record in records
        }

        locked_alignment = load_locked_alignment()

        locked_records = {
            get_accession(record): record
            for record in locked_alignment
        }

        new_records = [
            record
            for record in records
            if get_accession(record)
            not in locked_records
        ]

        # Start with the immutable locked alignment.
        if new_records:
            SeqIO.write(
                new_records,
                str(new_sequences_fasta),
                "fasta",
            )

            run_mafft_add(
                new_sequences_fasta,
                aligned_fasta,
            )

            alignment = AlignIO.read(
                str(aligned_fasta),
                "fasta",
            )

        else:
            alignment = locked_alignment

        mutations = detect_mutations_for_accessions(
            alignment,
            uploaded_accessions,
            variant_metadata,
        )

        if mutations.empty:
            mutation_type_count = 0
            mutation_accession_count = 0
        else:
            mutation_type_count = (
                mutations["mutation"].nunique()
            )

            mutation_accession_count = (
                mutations["accession"].nunique()
            )

        summary = {
            "input_sequences": len(records),
            "mutation_observations": len(mutations),
            "mutation_bearing_accessions": (
                mutation_accession_count
            ),
            "unique_mutation_types": (
                mutation_type_count
            ),
            "new_sequences_aligned": len(new_records),
        }

        alignment_text = StringIO()

        AlignIO.write(
            alignment,
            alignment_text,
            "fasta",
        )

        return {
            "mutations": mutations,
            "alignment_text": alignment_text.getvalue(),
            "summary": summary,
        }
