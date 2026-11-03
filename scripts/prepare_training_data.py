#!/usr/bin/env python3
"""Convert IEDB/TESLA public data into NeoFlux training format.

Usage:
    # Convert IEDB export (T-cell assay data)
    python scripts/prepare_training_data.py \\
        --input-iedb iedb_export.tsv \\
        --output training_data.tsv

    # Convert TESLA supplementary table
    python scripts/prepare_training_data.py \\
        --input-tesla tesla_data.tsv \\
        --output training_data.tsv

    # Compute missing features using NeoFlux predictors
    python scripts/prepare_training_data.py \\
        --input-iedb iedb_export.tsv \\
        --output training_data.tsv \\
        --compute-features

Input formats:
    IEDB: TSV with columns [peptide, allele, assay_type, response, ic50]
    TESLA: TSV with columns [patient, hla, peptide, mutation, immunogenicity, ic50]
"""

import argparse
import csv
import sys
from pathlib import Path

import numpy as np
import pandas as pd


def parse_iedb(path: Path) -> pd.DataFrame:
    rows = []
    with open(path) as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            peptide = (row.get("peptide") or row.get("Description") or "").strip()
            allele = (row.get("allele") or row.get("Allele") or "").strip()
            response = row.get("response") or row.get("Response") or row.get("Qualitative_Measure") or ""
            ic50_str = row.get("ic50") or row.get("IC50") or row.get("Quantitative_measure") or ""
            assay = row.get("assay_type") or row.get("Assay_Group") or ""

            if not peptide or not allele:
                continue

            is_positive = str(response).lower() in ("positive", "1", "true", "+", "yes")
            ic50 = float(ic50_str) if ic50_str and ic50_str not in ("NULL", "") else None
            if ic50 is not None and ic50 <= 0:
                ic50 = None

            rows.append({
                "peptide": peptide,
                "allele": allele,
                "label": 1 if is_positive else 0,
                "ic50": ic50,
                "assay": assay,
                "source": "iedb",
            })
    return pd.DataFrame(rows)


_GENERIC_ALLELES = frozenset({
    "hla class i", "hla class ii", "h2 class i", "h2 class ii",
    "h2-d class i", "h2-d class ii", "h2-b class i", "h2-b class ii",
    "h2-k class i", "h2-k class ii", "h2-l class i", "h2-l class ii",
    "hla-i", "hla-ii", "mhc class i", "mhc class ii",
    "rt1 class i", "rt1 class ii",
})


def _is_specific_allele(allele: str) -> bool:
    allele_clean = allele.strip().strip('"').lower()
    if not allele_clean or allele_clean in _GENERIC_ALLELES:
        return False
    if allele_clean.startswith(("hla-", "h2-", "h-2", "mamu-", "patr-", "gogo-")):
        return "*" in allele_clean or allele_clean.count(":") >= 1
    return False


def parse_iedb_v3(path: Path) -> pd.DataFrame:
    rows = []
    with open(path) as f:
        header_line = None
        for i, line in enumerate(f):
            if i == 1:
                header_line = line
                break
        if header_line is None:
            return pd.DataFrame(rows)

        def get(row_parts: list[str], idx_1based: int) -> str:
            idx = idx_1based - 1
            if 0 <= idx < len(row_parts):
                return row_parts[idx].strip().strip('"')
            return ""

        for line in f:
            line = line.rstrip("\n")
            if not line.strip():
                continue
            parts = line.split("\t")
            if len(parts) < 146:
                continue

            peptide = get(parts, 12)
            allele = get(parts, 142)
            mhc_class = get(parts, 146).upper()
            response = get(parts, 123)
            ic50_str = get(parts, 125)
            assay = get(parts, 119)

            if not peptide or not allele:
                continue
            if not _is_specific_allele(allele):
                continue

            is_positive = str(response).lower() in ("positive", "1", "true", "+", "yes")
            ic50 = None
            if ic50_str:
                try:
                    val = float(ic50_str)
                    if val > 0:
                        ic50 = val
                except ValueError:
                    pass

            rows.append({
                "peptide": peptide,
                "allele": allele,
                "label": 1 if is_positive else 0,
                "ic50": ic50,
                "assay": assay,
                "mhc_class": mhc_class,
                "source": "iedb_v3",
            })
    return pd.DataFrame(rows)


def parse_tesla(path: Path) -> pd.DataFrame:
    rows = []
    with open(path) as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            peptide = (row.get("peptide") or row.get("Peptide") or row.get("mutant_peptide") or "").strip()
            allele = (row.get("hla") or row.get("HLA") or row.get("allele") or "").strip()
            wt_peptide = (row.get("wildtype_peptide") or row.get("wt_peptide") or "").strip()
            immunogenicity = row.get("immunogenicity") or row.get("Immunogenicity") or row.get("response") or ""
            ic50_str = row.get("ic50") or row.get("IC50") or row.get("binding_affinity") or ""
            patient = (row.get("patient") or row.get("Patient") or row.get("sample_id") or "").strip()

            if not peptide or not allele:
                continue

            is_positive = str(immunogenicity).lower() in ("positive", "1", "true", "+", "yes", "responder")
            ic50 = float(ic50_str) if ic50_str and ic50_str not in ("NULL", "") else None
            if ic50 is not None and ic50 <= 0:
                ic50 = None

            rows.append({
                "peptide": peptide,
                "wildtype_peptide": wt_peptide,
                "allele": allele,
                "label": 1 if is_positive else 0,
                "ic50": ic50,
                "patient": patient,
                "source": "tesla",
            })
    return pd.DataFrame(rows)


def compute_mhcflurry_features(df: pd.DataFrame) -> pd.DataFrame:
    """Compute MHCflurry binding/presentation features per (peptide, allele) pair.

    Predictions are done per supported allele over the batch of peptides seen
    with that allele (one presentation-predictor call per allele). Pairs whose
    allele is unsupported by MHCflurry, or whose peptide length is outside the
    class I range (8-15), are left as NaN so the gradient-boosted models handle
    them natively as missing rather than dropping the record.
    """
    try:
        from mhcflurry import Class1PresentationPredictor
        presentation_predictor = Class1PresentationPredictor.load()
    except Exception as e:
        print(f"Warning: MHCflurry not available ({e}), skipping binding prediction", file=sys.stderr)
        return df

    supported = set(presentation_predictor.supported_alleles)

    # Normalise length to the class I window MHCflurry accepts.
    lengths = df["peptide"].str.len()
    valid_len = (lengths >= 8) & (lengths <= 15)

    # Build the set of unique (peptide, allele) pairs we can actually score.
    scorable = df[valid_len & df["allele"].isin(supported)]
    pairs = scorable[["peptide", "allele"]].drop_duplicates()
    n_alleles = pairs["allele"].nunique()
    print(
        f"  Predicting binding for {len(pairs)} unique pairs across "
        f"{n_alleles} supported alleles...",
        file=sys.stderr,
    )

    lookup: dict[tuple[str, str], dict[str, float]] = {}
    for processed, (allele, grp) in enumerate(pairs.groupby("allele"), start=1):
        peptides = grp["peptide"].tolist()
        try:
            out = presentation_predictor.predict(
                peptides=peptides,
                alleles={"s": [allele]},
                verbose=0,
            )
        except Exception as e:
            print(f"  skip allele {allele}: {e}", file=sys.stderr)
            continue
        for _, row in out.iterrows():
            lookup[(row["peptide"], allele)] = {
                "ic50": float(row["affinity"]),
                "percent_rank": float(row["presentation_percentile"]),
                "presentation_score": float(row["presentation_score"]),
            }
        if processed % 50 == 0:
            print(f"    {processed}/{n_alleles} alleles done", file=sys.stderr)

    ic50s, ranks, pres_scores = [], [], []
    for pep, allele in zip(df["peptide"], df["allele"]):
        hit = lookup.get((pep, allele))
        if hit is None:
            ic50s.append(None)
            ranks.append(None)
            pres_scores.append(None)
        else:
            ic50s.append(hit["ic50"])
            ranks.append(hit["percent_rank"])
            pres_scores.append(hit["presentation_score"])

    df["mhcflurry_ic50"] = ic50s
    df["mhcflurry_percent_rank"] = ranks
    df["mhcflurry_presentation_score"] = pres_scores

    use_ic50 = df["ic50"].isna() & df["mhcflurry_ic50"].notna()
    df.loc[use_ic50, "ic50"] = df.loc[use_ic50, "mhcflurry_ic50"]

    return df


def compute_processing_scores(df: pd.DataFrame) -> pd.DataFrame:
    from neoflux.antigen_processing import (
        predict_proteasome_cleavage,
        predict_stability_half_life,
        predict_tap_transport,
    )

    cleavage = []
    tap = []
    stability = []
    for peptide in df["peptide"]:
        cleavage.append(predict_proteasome_cleavage(peptide))
        tap.append(predict_tap_transport(peptide))
        ic50_val = None
        stability.append(predict_stability_half_life(peptide, ic50_val))

    df["proteasome_score"] = cleavage
    df["tap_score"] = tap
    df["stability_score"] = stability
    return df


def build_feature_matrix(df: pd.DataFrame) -> pd.DataFrame:
    from neoflux.features import FEATURE_COLUMNS
    from neoflux.utils import compute_hydrophobicity

    matrix = pd.DataFrame(index=df.index)

    peptide_lengths = df["peptide"].str.len()
    hydrophobicity = df["peptide"].apply(compute_hydrophobicity)

    if "mhcflurry_ic50" in df.columns:
        matrix["best_ic50"] = df["mhcflurry_ic50"].combine_first(df["ic50"])
    else:
        matrix["best_ic50"] = df.get("ic50")
    matrix["best_percent_rank"] = df.get("mhcflurry_percent_rank")
    matrix["is_strong_binder"] = (matrix["best_ic50"] < 50).astype(float)
    matrix["is_weak_binder"] = ((matrix["best_ic50"] >= 50) & (matrix["best_ic50"] < 500)).astype(float)
    matrix["presentation_score"] = df.get("mhcflurry_presentation_score")
    matrix["pmhc_half_life"] = df.get("stability_score")
    matrix["binding_promiscuity"] = 0.0

    matrix["gene_tpm"] = None
    matrix["transcript_tpm"] = None
    matrix["vaf"] = None
    matrix["rna_vaf"] = None
    matrix["mutant_expression"] = None

    matrix["is_clonal"] = 0.0
    matrix["cancer_cell_fraction"] = None

    matrix["is_snv"] = 1.0
    matrix["is_indel"] = 0.0
    matrix["is_frameshift"] = 0.0
    matrix["is_fusion"] = 0.0
    matrix["is_splice"] = 0.0
    matrix["is_ncorf"] = 0.0

    matrix["peptide_length"] = peptide_lengths.astype(float)
    matrix["hydrophobicity"] = hydrophobicity
    matrix["self_similarity"] = None
    matrix["agretopicity"] = None

    if {"proteasome_score", "tap_score"}.issubset(df.columns):
        matrix["proteasome_score"] = df["proteasome_score"]
        matrix["tap_score"] = df["tap_score"]
    else:
        matrix["proteasome_score"] = None
        matrix["tap_score"] = None

    matrix["presenting_allele_loh"] = 0.0
    matrix["b2m_mutated"] = 0.0

    matrix["ms_detected"] = None
    matrix["pepquery_score"] = None
    matrix["spectral_count"] = None

    matrix["is_driver_gene"] = 0.0
    matrix["label"] = df["label"].values

    missing = [c for c in FEATURE_COLUMNS if c not in matrix.columns]
    if missing:
        for col in missing:
            matrix[col] = None
    return matrix[FEATURE_COLUMNS + ["label"]]


def main():
    parser = argparse.ArgumentParser(description="Prepare IEDB/TESLA training data for NeoFlux")
    parser.add_argument("--input-iedb", help="IEDB export TSV (simple column header format)")
    parser.add_argument("--input-iedb-v3", help="IEDB v3 TSV export (2-row header, 161 columns)")
    parser.add_argument("--input-tesla", help="TESLA supplementary table TSV")
    parser.add_argument("--output", required=True)
    parser.add_argument("--compute-features", action="store_true", help="Run MHCflurry and processing predictors")
    args = parser.parse_args()

    if args.input_iedb_v3:
        print(f"Reading IEDB v3 data from {args.input_iedb_v3}...", file=sys.stderr)
        df = parse_iedb_v3(Path(args.input_iedb_v3))
    elif args.input_iedb:
        print(f"Reading IEDB data from {args.input_iedb}...", file=sys.stderr)
        df = parse_iedb(Path(args.input_iedb))
    elif args.input_tesla:
        print(f"Reading TESLA data from {args.input_tesla}...", file=sys.stderr)
        df = parse_tesla(Path(args.input_tesla))
    else:
        parser.error("Provide either --input-iedb, --input-iedb-v3, or --input-tesla")

    if df.empty:
        print("Error: No valid records found in input", file=sys.stderr)
        sys.exit(1)

    print(f"  Found {len(df)} records ({df['label'].sum()} positive, {(1 - df['label']).sum()} negative)", file=sys.stderr)

    if args.compute_features:
        print("Computing MHCflurry binding predictions...", file=sys.stderr)
        df = compute_mhcflurry_features(df)
        print("Computing antigen processing scores...", file=sys.stderr)
        df = compute_processing_scores(df)

    print("Building feature matrix...", file=sys.stderr)
    matrix = build_feature_matrix(df)
    matrix.insert(0, "peptide", df["peptide"].values)
    matrix.insert(1, "allele", df["allele"].values)
    if "patient" in df.columns:
        matrix.insert(0, "patient", df["patient"].values)

    matrix.to_csv(args.output, sep="\t", index=False)
    print(f"Wrote {len(matrix)} training records to {args.output}", file=sys.stderr)
    print(f"  Features: {len(matrix.columns)} columns, {len([c for c in matrix.columns if c != 'label'])} predictors", file=sys.stderr)


if __name__ == "__main__":
    main()
