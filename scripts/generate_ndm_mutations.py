#!/usr/bin/env python3

from pathlib import Path
from collections import Counter
import re
import sys

import pandas as pd
from Bio import AlignIO


# ============================================================
# CONFIGURATION
# ============================================================

BASE = Path(__file__).resolve().parent.parent

MSA_FILE = BASE / "analysis/msa/NDM_MSA.fasta"
EXISTING_FILE = BASE / "analysis/mutations/NDM_mutations.csv"

REFERENCE_ID = "CBZ39493.1"

OUTPUT_FILE = BASE / "analysis/mutations/NDM_mutations_REGENERATED.csv"
REPORT_FILE = BASE / "audit/mutation_generation/mutation_generation_dry_run.txt"


# ============================================================
# HELPERS
# ============================================================

def get_accession(record):
    """
    Normalize FASTA identifiers.

    Standard NCBI:
        AFB82585.1

    UniProt:
        tr|A0A0K1Z5K4|A0A0K1Z5K4_KLEPN
    """
    if record.id.startswith("tr|"):
        parts = record.id.split("|")
        if len(parts) >= 2:
            return parts[1]

    return record.id


def get_variant(description):
    """
    Extract NDM variant number from the complete FASTA header.

    This handles both NCBI-style and UniProt-style headers,
    including:
        NDM-17
        GN=blaNDM-16
    """
    match = re.search(r"NDM-\d+", description, re.IGNORECASE)

    if match:
        return match.group(0).upper()

    return "NOT_FOUND"


def build_reference_coordinate_map(reference_sequence):
    """
    Map each MSA column containing a reference residue to the
    corresponding ungapped biological reference position.

    Returns:
        {msa_column: (reference_position, reference_residue)}
    """

    mapping = {}
    reference_position = 0

    for msa_column, residue in enumerate(
        str(reference_sequence.seq),
        start=1
    ):
        if residue != "-":
            reference_position += 1
            mapping[msa_column] = (
                reference_position,
                residue
            )

    return mapping


def generate_mutations(alignment, reference):
    """
    Detect amino-acid substitutions relative to CBZ39493.1.

    Gaps in the observed sequence are ignored here because
    deletions are represented separately by the pipeline.
    """

    coordinate_map = build_reference_coordinate_map(reference)

    records = []

    for record in alignment:

        if record.id == reference.id:
            continue

        accession = get_accession(record)
        variant = get_variant(record.description)

        for msa_column, (reference_position, reference_residue) in coordinate_map.items():

            observed_residue = str(record.seq)[msa_column - 1]

            # Ignore alignment gaps.
            if observed_residue == "-":
                continue

            # Record substitutions only.
            if observed_residue != reference_residue:

                records.append({
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
                    "msa_column": msa_column
                })

    return pd.DataFrame(records)


def mutation_key(row):
    return (
        row["accession"],
        row["reference_position"],
        row["reference_residue"],
        row["observed_residue"],
        row["mutation"],
        row["msa_column"],
    )


# ============================================================
# MAIN
# ============================================================

def main():

    print("=" * 60)
    print("NDM MUTATION GENERATOR — DRY RUN")
    print("=" * 60)

    if not MSA_FILE.exists():
        raise SystemExit(f"ERROR: MSA not found: {MSA_FILE}")

    if not EXISTING_FILE.exists():
        raise SystemExit(
            f"ERROR: existing mutation file not found: {EXISTING_FILE}"
        )

    alignment = AlignIO.read(MSA_FILE, "fasta")
    existing = pd.read_csv(EXISTING_FILE)

    reference = None

    for record in alignment:
        if record.id == REFERENCE_ID:
            reference = record
            break

    if reference is None:
        raise SystemExit(
            f"ERROR: reference {REFERENCE_ID} not found in MSA"
        )

    generated = generate_mutations(alignment, reference)

    # --------------------------------------------------------
    # Basic checks
    # --------------------------------------------------------

    print()
    print("=== INPUT ===")
    print(f"MSA records:        {len(alignment)}")
    print(f"MSA length:         {alignment.get_alignment_length()}")
    print(f"Reference:          {REFERENCE_ID}")
    print(f"Reference length:   {len(reference.seq)}")

    print()
    print("=== REGENERATED DATA ===")
    print(f"Mutation rows:      {len(generated)}")
    print(f"Accessions:         {generated['accession'].nunique()}")
    print(f"Mutation types:     {generated['mutation'].nunique()}")

    print()
    print("=== EXISTING DATA ===")
    print(f"Mutation rows:      {len(existing)}")
    print(f"Accessions:         {existing['accession'].nunique()}")
    print(f"Mutation types:     {existing['mutation'].nunique()}")

    # --------------------------------------------------------
    # Exact mutation comparison
    # --------------------------------------------------------

    generated_counter = Counter(
        mutation_key(row)
        for _, row in generated.iterrows()
    )

    existing_counter = Counter(
        mutation_key(row)
        for _, row in existing.iterrows()
    )

    mutation_match = (
        generated_counter == existing_counter
    )

    print()
    print("=== MUTATION DETECTION ===")

    if mutation_match:
        print("PASS")
        print(
            "The independently regenerated mutation observations "
            "exactly match the existing dataset."
        )
    else:
        print("FAIL")
        print("Mutation observations differ.")

    # --------------------------------------------------------
    # Variant comparison
    # --------------------------------------------------------

    header_variants = {}

    for record in alignment:
        accession = get_accession(record)
        variant = get_variant(record.description)
        header_variants[accession] = variant

    variant_mismatches = []

    for accession in sorted(header_variants):

        generated_rows = generated[
            generated["accession"] == accession
        ]

        existing_rows = existing[
            existing["accession"] == accession
        ]

        generated_variant = (
            generated_rows["variant"].iloc[0]
            if len(generated_rows) > 0
            else header_variants[accession]
        )

        existing_variants = sorted(
            existing_rows["variant"].dropna().unique().tolist()
        )

        if (
            len(existing_variants) > 0
            and existing_variants != [generated_variant]
        ):
            variant_mismatches.append({
                "accession": accession,
                "expected": generated_variant,
                "existing": ",".join(existing_variants)
            })

    print()
    print("=== VARIANT ANNOTATION ===")

    if variant_mismatches:

        for item in variant_mismatches:
            print(
                f"{item['accession']}: "
                f"expected={item['expected']} "
                f"existing={item['existing']}"
            )

    else:
        print("No variant-label mismatches found.")

    # --------------------------------------------------------
    # Critical NDM-16 check
    # --------------------------------------------------------

    ndm16_generated = generated[
        generated["accession"] == "A0A0K1Z5K4"
    ]

    print()
    print("=== NDM-16 CRITICAL CHECK ===")

    if len(ndm16_generated) == 0:
        print("FAIL: A0A0K1Z5K4 has no mutation record.")
    else:
        print(
            ndm16_generated[
                ["accession", "variant", "mutation", "reference_position", "msa_column"]
            ].to_string(index=False)
        )

        if set(ndm16_generated["variant"]) == {"NDM-16"}:
            print("PASS: A0A0K1Z5K4 correctly resolves to NDM-16.")
        else:
            print("FAIL: NDM-16 variant extraction failed.")

    # --------------------------------------------------------
    # Unknown check
    # --------------------------------------------------------

    unknown_existing = int(
        (existing["variant"] == "Unknown").sum()
    )

    unknown_generated = int(
        (generated["variant"] == "Unknown").sum()
    )

    print()
    print("=== UNKNOWN VARIANT CHECK ===")
    print(f"Existing Unknown records:    {unknown_existing}")
    print(f"Regenerated Unknown records: {unknown_generated}")

    # --------------------------------------------------------
    # Save report only
    # --------------------------------------------------------

    REPORT_FILE.parent.mkdir(parents=True, exist_ok=True)

    with open(REPORT_FILE, "w") as report:

        report.write(
            "NDM MUTATION GENERATION DRY-RUN REPORT\n"
        )
        report.write("=" * 50 + "\n\n")

        report.write(f"MSA records: {len(alignment)}\n")
        report.write(
            f"MSA length: {alignment.get_alignment_length()}\n"
        )
        report.write(
            f"Reference: {REFERENCE_ID}\n"
        )
        report.write(
            f"Reference length: {len(reference.seq)}\n\n"
        )

        report.write(
            f"Existing mutation rows: {len(existing)}\n"
        )
        report.write(
            f"Regenerated mutation rows: {len(generated)}\n"
        )
        report.write(
            f"Exact mutation match: {mutation_match}\n"
        )
        report.write(
            f"Existing Unknown: {unknown_existing}\n"
        )
        report.write(
            f"Regenerated Unknown: {unknown_generated}\n"
        )

        report.write("\nVariant mismatches:\n")

        if variant_mismatches:
            for item in variant_mismatches:
                report.write(
                    f"{item['accession']} "
                    f"expected={item['expected']} "
                    f"existing={item['existing']}\n"
                )
        else:
            report.write("None\n")

    # --------------------------------------------------------
    # IMPORTANT:
    # Do not overwrite the current mutation file.
    # --------------------------------------------------------

    generated.to_csv(
        OUTPUT_FILE,
        index=False
    )

    print()
    print("=== DRY-RUN OUTPUT ===")
    print(f"Generated comparison file:")
    print(OUTPUT_FILE)
    print(f"Audit report:")
    print(REPORT_FILE)

    print()
    print("IMPORTANT:")
    print("The existing NDM_mutations.csv was NOT modified.")

    # --------------------------------------------------------
    # Exit status
    # --------------------------------------------------------

    if not mutation_match:
        sys.exit(1)

    if unknown_generated != 0:
        sys.exit(1)

    if not (
        len(ndm16_generated) > 0
        and set(ndm16_generated["variant"]) == {"NDM-16"}
    ):
        sys.exit(1)

    print()
    print("DRY-RUN VALIDATION: PASS")


if __name__ == "__main__":
    main()
