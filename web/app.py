from io import StringIO

import pandas as pd
import streamlit as st
from Bio import SeqIO

from ndm_backend import analyse_fasta_bytes


# ============================================================
# PAGE CONFIGURATION
# ============================================================

st.set_page_config(
    page_title="NDM Mutation Analysis",
    page_icon="🧬",
    layout="wide",
)


# ============================================================
# TITLE
# ============================================================

st.title("NDM Mutation Analysis Pipeline")

st.write(
    "Automated analysis of amino-acid mutations in "
    "NDM-family β-lactamase protein sequences."
)

st.info(
    "Automated analysis of amino-acid mutations in NDM-family β-lactamase protein sequences."
)


# ============================================================
# FASTA VALIDATION
# ============================================================

ALLOWED_RESIDUES = set(
    "ACDEFGHIKLMNPQRSTVWYBXZJUO"
)


def validate_fasta(uploaded_bytes):
    """
    Validate uploaded FASTA data before mutation analysis.
    """

    errors = []
    records = []

    # --------------------------------------------------------
    # UTF-8 validation
    # --------------------------------------------------------

    try:
        fasta_text = uploaded_bytes.decode("utf-8")
    except UnicodeDecodeError:
        return {
            "valid": False,
            "errors": ["File is not valid UTF-8 text."],
            "records": [],
            "count": 0,
            "min_length": None,
            "max_length": None,
            "mean_length": None,
        }

    # --------------------------------------------------------
    # FASTA header check
    # --------------------------------------------------------

    if not fasta_text.strip().startswith(">"):
        errors.append(
            "File does not appear to be a FASTA file."
        )
        return {
            "valid": False,
            "errors": errors,
            "records": [],
            "count": 0,
            "min_length": None,
            "max_length": None,
            "mean_length": None,
        }

    # --------------------------------------------------------
    # Parse sequences
    # --------------------------------------------------------

    try:
        records = list(
            SeqIO.parse(
                StringIO(fasta_text),
                "fasta",
            )
        )
    except Exception as exc:
        errors.append(
            f"FASTA parsing failed: {exc}"
        )

    if not records:
        errors.append(
            "No FASTA sequences were detected."
        )

    # --------------------------------------------------------
    # Duplicate identifiers
    # --------------------------------------------------------

    identifiers = [
        record.id
        for record in records
    ]

    duplicate_ids = sorted(
        {
            identifier
            for identifier in identifiers
            if identifiers.count(identifier) > 1
        }
    )

    if duplicate_ids:
        errors.append(
            "Duplicate sequence identifiers detected: "
            + ", ".join(duplicate_ids)
        )

    # --------------------------------------------------------
    # Sequence validation
    # --------------------------------------------------------

    lengths = []

    for record in records:

        sequence = str(record.seq).upper()

        if not sequence:
            errors.append(
                f"Empty sequence detected: {record.id}"
            )
            continue

        lengths.append(len(sequence))

        invalid_residues = sorted(
            set(sequence) - ALLOWED_RESIDUES
        )

        if invalid_residues:
            errors.append(
                f"{record.id}: invalid residue(s): "
                + ", ".join(invalid_residues)
            )

    # --------------------------------------------------------
    # Metrics
    # --------------------------------------------------------

    if lengths:
        min_length = min(lengths)
        max_length = max(lengths)
        mean_length = sum(lengths) / len(lengths)
    else:
        min_length = None
        max_length = None
        mean_length = None

    return {
        "valid": len(errors) == 0,
        "errors": errors,
        "records": records,
        "count": len(records),
        "min_length": min_length,
        "max_length": max_length,
        "mean_length": mean_length,
    }


# ============================================================
# UPLOAD
# ============================================================

st.subheader("1. Upload NDM protein sequences")

uploaded_file = st.file_uploader(
    "Upload a FASTA file containing NDM protein sequences",
    type=["fasta", "fa", "faa"],
)


if uploaded_file is None:

    st.info(
        "Upload an NDM protein FASTA file to begin."
    )

else:

    uploaded_bytes = uploaded_file.getvalue()

    st.write(
        f"**File:** {uploaded_file.name}"
    )

    st.write(
        f"**File size:** {len(uploaded_bytes):,} bytes"
    )

    # --------------------------------------------------------
    # Validate
    # --------------------------------------------------------

    st.subheader("2. FASTA input check")

    validation = validate_fasta(uploaded_bytes)

    if validation["valid"]:

        st.success(
            "FASTA validation passed."
        )

        metric1, metric2, metric3 = st.columns(3)

        metric1.metric(
            "Sequences",
            validation["count"],
        )

        metric2.metric(
            "Minimum length",
            f"{validation['min_length']} aa",
        )

        metric3.metric(
            "Maximum length",
            f"{validation['max_length']} aa",
        )

        st.write(
            f"Mean sequence length: "
            f"{validation['mean_length']:.1f} aa"
        )

        st.write(
            "✓ FASTA structure detected"
        )

        st.write(
            "✓ Sequences successfully parsed"
        )

        st.write(
            "✓ Protein residue validation passed"
        )

        st.write(
            "✓ No duplicate sequence identifiers"
        )

        st.write(
            "✓ No empty sequences"
        )

    else:

        st.error(
            "FASTA validation failed."
        )

        for error in validation["errors"]:
            st.error(error)

        st.stop()

    # ========================================================
    # ALIGNMENT AND MUTATION DETECTION
    # ========================================================

    st.divider()

    st.subheader("3. Multiple sequence alignment and mutation detection")

    st.write(
        "The validated NDM mutation-analysis backend aligns the "
        "uploaded sequences against the locked NDM reference framework "
        "and detects amino-acid substitutions relative to the reference."
    )

    run_analysis = st.button(
        "Run Alignment and Mutation Detection",
        type="primary",
    )

    if run_analysis:

        with st.spinner(
            "Running multiple sequence alignment and mutation detection..."
        ):

            try:
                result = analyse_fasta_bytes(
                    uploaded_bytes
                )

            except Exception as exc:

                st.error(
                    "Alignment and mutation detection failed."
                )

                st.exception(exc)

                st.stop()

        # The mutation table must exist before validation or display.
        mutations = result["mutations"].copy()
        summary = result["summary"]

        st.success(
            "Multiple sequence alignment and mutation detection completed successfully."
        )

        # ----------------------------------------------------
        # Tool validation
        # ----------------------------------------------------

        st.subheader("4. Mutation results")

        st.caption(
            "Summary and detailed mutation results generated from the "
            "uploaded NDM protein sequences."
        )

        st.markdown("### Processing summary")

        summary_data = {
            "Metric": [
                "Input sequences",
                "Mutation observations",
                "Mutation-bearing sequences",
                "Unique mutation types",
            ],
            "Value": [
                summary["input_sequences"],
                summary["mutation_observations"],
                summary["mutation_bearing_accessions"],
                summary["unique_mutation_types"],
            ],
        }

        st.dataframe(
            pd.DataFrame(summary_data),
            width="stretch",
            hide_index=True,
        )

        st.markdown("### Detected mutations")

        if mutations.empty:

            st.info(
                "No amino-acid substitutions were detected in the "
                "uploaded sequences relative to the selected reference."
            )

        else:

            display_columns = [
                "accession",
                "variant",
                "reference_position",
                "reference_residue",
                "observed_residue",
                "mutation",
                "msa_column",
            ]

            mutations = mutations[
                display_columns
            ]

            st.dataframe(
                mutations,
                width="stretch",
                hide_index=True,
                column_config={
                    "accession": st.column_config.TextColumn(
                        "Accession"
                    ),
                    "variant": st.column_config.TextColumn(
                        "NDM variant"
                    ),
                    "reference_position": st.column_config.NumberColumn(
                        "Reference position"
                    ),
                    "reference_residue": st.column_config.TextColumn(
                        "Reference residue"
                    ),
                    "observed_residue": st.column_config.TextColumn(
                        "Observed residue"
                    ),
                    "mutation": st.column_config.TextColumn(
                        "Mutation"
                    ),
                    "msa_column": st.column_config.NumberColumn(
                        "MSA column"
                    ),
                },
            )

            csv_data = mutations.to_csv(
                index=False
            ).encode("utf-8")

            st.download_button(
                label="Download mutation results (CSV)",
                data=csv_data,
                file_name="NDM_mutation_results.csv",
                mime="text/csv",
            )

            # ------------------------------------------------
            # Mutation frequency
            # ------------------------------------------------

            st.subheader("5. Mutation frequency")

            st.caption(
                "Live frequency summary for the currently uploaded sequences. "
                "Formal dataset-level frequency analysis and visualisation are "
                "provided in the R-based analysis below."
            )

            mutation_counts = (
                mutations["mutation"]
                .value_counts()
                .rename_axis("mutation")
                .reset_index(name="count")
            )

            mutation_counts["percentage"] = (
                mutation_counts["count"]
                / summary["input_sequences"]
                * 100
            )

            mutation_counts["percentage"] = (
                mutation_counts["percentage"]
                .round(2)
            )

            mutation_counts["percentage"] = (
                mutation_counts["percentage"].map(lambda x: f"{x:.2f}%")
            )

            mutation_counts = mutation_counts.head(10)

            st.dataframe(
                mutation_counts,
                width='stretch',
                hide_index=True,
            )

        # ----------------------------------------------------
        # NDM-16 check
        # ----------------------------------------------------

        st.subheader("6. Tool validation")

        st.caption(
            "The mutation-detection tool was checked against selected "
            "characterised NDM variants with published mutation records."
        )

        validation_cases = [
            {
                "accession": "AFB82585.1",
                "variant": "NDM-4",
                "expected": ["M154L"],
            },
            {
                "accession": "AXH80118.1",
                "variant": "NDM-5",
                "expected": ["V88L", "M154L", "A263P"],
            },
            {
                "accession": "AOT73359.1",
                "variant": "NDM-17",
                "expected": ["V88L", "M154L", "E170K"],
            },
        ]

        validation_rows = []

        for case in validation_cases:
            accession = case["accession"]
            expected = set(case["expected"])

            observed_rows = mutations[
                mutations["accession"] == accession
            ]

            if observed_rows.empty:
                validation_rows.append(
                    {
                        "Accession": accession,
                        "NDM variant": case["variant"],
                        "Expected mutation(s)": ", ".join(case["expected"]),
                        "Detected mutation(s)": "Not present",
                        "Result": "NOT EVALUATED",
                    }
                )
                continue

            observed = set(
                observed_rows["mutation"].astype(str).tolist()
            )

            missing = expected - observed

            if not missing:
                validation_result = "PASS"
            else:
                validation_result = "FAIL"

            validation_rows.append(
                {
                    "Accession": accession,
                    "NDM variant": case["variant"],
                    "Expected mutation(s)": ", ".join(case["expected"]),
                    "Detected mutation(s)": ", ".join(sorted(observed)),
                    "Result": validation_result,
                }
            )

        validation_df = pd.DataFrame(validation_rows)

        st.dataframe(
            validation_df,
            width="stretch",
            hide_index=True,
        )

        evaluated = validation_df[
            validation_df["Result"] != "NOT EVALUATED"
        ]

        failed = evaluated[
            evaluated["Result"] == "FAIL"
        ]

        if evaluated.empty:
            st.info(
                "No validation cases were present in this upload."
            )
        elif failed.empty:
            st.success(
                "Validation PASS — all expected characterised "
                "mutations were detected in the evaluated cases."
            )
        else:
            st.error(
                f"Validation FAIL — {len(failed)} evaluated case(s) "
                "did not contain all expected mutations."
            )

        st.caption(
            "PASS means the published expected mutation(s) were detected. "
            "Additional substitutions do not constitute a failure because "
            "the pipeline reports all observed substitutions relative to "
            "the selected reference sequence."
        )

        # ----------------------------------------------------
        # Mutation results
        # ----------------------------------------------------

        st.subheader("7. R-based visual analysis")

        st.caption(
            "Validated visual summaries generated from the locked "
            "NDM mutation-analysis dataset."
        )

        st.markdown("### 7.1 Variant composition")

        st.image(
            "results/figures/figure1_variant_composition.png",
            width=700,
        )

        st.caption(
            "What this shows: the number of sequences represented by "
            "each NDM variant in the analysed dataset."
        )

        st.markdown("### 7.2 Most frequently detected mutations")

        st.image(
            "results/figures/figure2_top_mutations.png",
            width=700,
        )

        st.caption(
            "What this shows: the most frequently observed amino-acid "
            "substitutions and their frequencies within the dataset."
        )

        st.markdown("### 7.3 Mutation landscape")

        st.image(
            "results/figures/figure3_mutation_landscape.png",
            width=700,
        )

        st.caption(
            "What this shows: the positions and types of amino-acid "
            "substitutions observed relative to the reference sequence."
        )

        st.markdown("### 7.4 Mutation burden")

        st.image(
            "results/figures/figure4_mutation_burden.png",
            width=700,
        )

        st.caption(
            "What this shows: the number of detected mutations carried "
            "by each sequence in the analysed dataset."
        )

        st.markdown("### 7.5 Variant–mutation heatmap")

        st.image(
            "results/figures/figure5_variant_mutation_heatmap.png",
            width=700,
        )

        st.caption(
            "What this shows: the distribution of selected mutations "
            "across the NDM variants represented in the dataset."
        )

        st.markdown("### 7.6 Mutation co-occurrence")

        st.image(
            "results/figures/figure6_mutation_cooccurrence.png",
            width=550,
        )

        st.caption(
            "What this shows: mutation pairs that were observed together "
            "within the same sequences."
        )

        st.markdown("### 7.7 Integrated mutation summary")

        st.image(
            "results/figures/figure7_ndm_mutation_summary.png",
            width=700,
        )

        st.caption(
            "What this shows: an integrated summary of the main mutation "
            "patterns and dataset-level statistics."
        )

        # ----------------------------------------------------
        # Completion note
        # ----------------------------------------------------

        st.divider()

        st.success(
            "Analysis complete. Results were generated "
            "using the validated NDM mutation-analysis backend."
        )

        st.caption(
            "NDM mutation analysis completed using the validated NDM-specific analysis workflow. "
            "Results should be interpreted within the "
            "curated study dataset and reference framework."
        )
