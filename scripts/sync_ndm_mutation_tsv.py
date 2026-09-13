#!/usr/bin/env python3

from pathlib import Path
import pandas as pd

BASE = Path(__file__).resolve().parent.parent

CSV_FILE = BASE / "analysis/mutations/NDM_mutations.csv"
TSV_FILE = BASE / "analysis/mutations/NDM_mutations.tsv"
REPORT_FILE = BASE / "audit/mutation_generation/mutation_tsv_sync_report.txt"

EXPECTED_COLUMNS = [
    "accession",
    "variant",
    "reference_position",
    "reference_residue",
    "observed_residue",
    "mutation",
    "msa_column",
]

df = pd.read_csv(CSV_FILE)

if list(df.columns) != EXPECTED_COLUMNS:
    raise SystemExit(
        "FAIL: Unexpected mutation CSV columns:\n"
        + str(list(df.columns))
    )

df.to_csv(TSV_FILE, sep="\t", index=False)

# Re-read the written TSV for independent verification.
tsv = pd.read_csv(TSV_FILE, sep="\t")

checks = {}

checks["row_count"] = len(df) == 106 and len(tsv) == 106
checks["accession_count"] = (
    df["accession"].nunique() == 27
    and tsv["accession"].nunique() == 27
)
checks["mutation_type_count"] = (
    df["mutation"].nunique() == 82
    and tsv["mutation"].nunique() == 82
)
checks["unknown_count"] = (
    (df["variant"] == "Unknown").sum() == 0
    and (tsv["variant"] == "Unknown").sum() == 0
)

ndm16 = tsv[tsv["accession"] == "A0A0K1Z5K4"]

checks["ndm16_record"] = (
    len(ndm16) == 1
    and ndm16.iloc[0]["variant"] == "NDM-16"
    and ndm16.iloc[0]["mutation"] == "R264H"
    and int(ndm16.iloc[0]["reference_position"]) == 264
    and int(ndm16.iloc[0]["msa_column"]) == 293
)

# Compare complete tabular content.
csv_normalized = df.astype(str).reset_index(drop=True)
tsv_normalized = tsv.astype(str).reset_index(drop=True)

checks["csv_tsv_identical"] = csv_normalized.equals(tsv_normalized)

report = []
report.append("NDM MUTATION TSV SYNCHRONIZATION REPORT")
report.append("=" * 60)
report.append("")
report.append(f"CSV source: {CSV_FILE}")
report.append(f"TSV output: {TSV_FILE}")
report.append("")
report.append(f"Rows: {len(tsv)}")
report.append(f"Accessions: {tsv['accession'].nunique()}")
report.append(f"Mutation types: {tsv['mutation'].nunique()}")
report.append(
    f"Unknown variants: {(tsv['variant'] == 'Unknown').sum()}"
)
report.append("")
report.append("NDM-16 record:")
report.append(ndm16.to_string(index=False))
report.append("")
report.append("CHECKS")
report.append("-" * 60)

for name, result in checks.items():
    report.append(f"{name}: {'PASS' if result else 'FAIL'}")

all_pass = all(checks.values())

report.append("")
report.append(
    "FINAL VERDICT: "
    + ("PASS — TSV synchronized exactly with authoritative CSV."
       if all_pass
       else "FAIL — synchronization check failed.")
)

REPORT_FILE.write_text("\n".join(report) + "\n")

if not all_pass:
    print("\n".join(report))
    raise SystemExit(1)

print("\n".join(report))
