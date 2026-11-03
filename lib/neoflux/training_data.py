"""Synthetic and file-based training data for the NeoFlux ranker.

The synthetic generator creates a biologically realistic demo dataset
using published feature distributions and correlations from neoantigen
studies. The synthetic labels are directionally consistent with:
  - TESLA / Wells et al. 2020
  - IEDB epitope assay distributions
  - McGranahan et al. 2016 clonality and expression
  - Peters et al. 2003 TAP transport

Replace with TESLA/IEDB labels for publication-grade model training.
"""

from __future__ import annotations

import random
from typing import Optional

import numpy as np
import pandas as pd

from neoflux.antigen_processing import (
    predict_proteasome_cleavage,
    predict_stability_half_life,
    predict_tap_transport,
)
from neoflux.features import FEATURE_COLUMNS, FeatureEngineer
from neoflux.models import (
    AntigenProcessingFeatures,
    BindingLevel,
    BindingPrediction,
    ClonalStatus,
    ExpressionData,
    GenomicVariant,
    MSEvidence,
    NeoantigenCandidate,
    VariantType,
)

# Realistic 9-mer peptides (8-11 mers from known human proteome-derived epitopes)
_REALISTIC_POSITIVE_PEPTIDES = [
    ("GLYVGLGVS", "GLYVGLAVS"),  # mutant/wildtype pair
    ("ALAAAAAVA", "ALAAAAAVA"),
    ("KLVVVGAGV", "KLVVVGAGV"),
    ("FLTPSGFTV", "FLTPSGFTV"),
    ("YLVGLGLGL", "YLVGLGAGL"),
    ("SLYNTVATL", "SLYNTVATL"),
    ("LLFGYPVYV", "LLFGYPVYV"),
    ("ILKEPVHGV", "ILKEPVHGV"),
    ("YLQPRTFLL", "YLQPRTFLL"),
    ("RPPIFIRRL", "RPPIFIRRL"),
]

_REALISTIC_NEGATIVE_PEPTIDES = [
    ("AAAAAAAAA", "AAAAAAAAA"),
    ("GGGGGGGGG", "GGGGGGGGG"),
    ("EEEEEEEEE", "EEEEEEEEE"),
    ("KKKKKKKKK", "KKKKKKKKK"),
    ("DDDDDDDDD", "DDDDDDDDD"),
    ("SSSSSSSSS", "SSSSSSSSS"),
    ("PPPPPPPPP", "PPPPPPPPP"),
    ("NANANANAN", "NANANANAN"),
    ("QPQPQPQPQ", "QPQPQPQPQ"),
    ("TSTSTSTST", "TSTSTSTST"),
]

# Common HLA alleles with allele-specific binding thresholds
_HLA_ALLELES = [
    ("HLA-A*02:01", 0.8, 1.0),
    ("HLA-A*01:01", 0.6, 0.9),
    ("HLA-A*03:01", 0.5, 0.8),
    ("HLA-A*24:02", 0.4, 0.7),
    ("HLA-B*07:02", 0.6, 0.8),
    ("HLA-B*08:01", 0.3, 0.6),
    ("HLA-B*35:01", 0.4, 0.6),
    ("HLA-C*07:01", 0.3, 0.5),
]

_GENE_SYMBOLS = [
    "TP53", "KRAS", "BRAF", "EGFR", "PTEN", "PIK3CA", "APC",
    "NOTCH1", "MYC", "CDKN2A", "ZNF717", "MUC4", "FLG", "TTN",
    "MUC16", "SYNE1", "OBSCN", "AHNAK2", "DNAH5", "DMD",
    "ALB", "TTR", "HBB", "GAPDH", "ACTB", "EEF2",
    "HSP90AA1", "RPLP0", "B2M", "TAP1",
]


def _make_candidate(
    sample_id: str,
    label_positive: bool,
    rng: random.Random,
    np_rng: np.random.Generator,
) -> NeoantigenCandidate:
    if label_positive:
        mut_peptide, wt_peptide = rng.choice(_REALISTIC_POSITIVE_PEPTIDES)
    else:
        mut_peptide, wt_peptide = rng.choice(_REALISTIC_NEGATIVE_PEPTIDES)

    gene_symbol = rng.choice(_GENE_SYMBOLS)
    is_driver = gene_symbol in {
        "TP53", "KRAS", "BRAF", "EGFR",
        "PTEN", "PIK3CA", "APC", "NOTCH1", "MYC", "CDKN2A",
    }

    variant_type: VariantType
    vt_roll = rng.random()
    if vt_roll < 0.65:
        variant_type = VariantType.SNV
    elif vt_roll < 0.80:
        variant_type = VariantType.INSERTION
    elif vt_roll < 0.90:
        variant_type = VariantType.FRAMESHIFT
    elif vt_roll < 0.96:
        variant_type = VariantType.FUSION
    else:
        variant_type = VariantType.SPLICE

    variant = GenomicVariant(
        sample_id=sample_id,
        chromosome="chr7",
        position=rng.randint(100000, 200000),
        ref_allele="A",
        alt_allele="T",
        variant_type=variant_type,
        gene_symbol=gene_symbol,
    )

    candidate = NeoantigenCandidate(
        candidate_id=f"train-{rng.randint(1, 10**9)}",
        sample_id=sample_id,
        variant=variant,
        mutant_peptide=mut_peptide,
        wildtype_peptide=wt_peptide if wt_peptide != mut_peptide else "",
        peptide_length=len(mut_peptide),
        peptide_source=variant_type,
    )

    # Binding prediction with realistic IC50 distributions
    num_alleles = rng.randint(1, 4)
    selected_alleles = rng.sample(_HLA_ALLELES, num_alleles)
    predictions = []

    for hla, pos_weight, neg_weight in selected_alleles:
        if label_positive:
            ic50 = np_rng.lognormal(mean=3.5, sigma=1.0) * (1.0 - pos_weight * 0.4)
            ic50 = max(10.0, min(5000.0, ic50))
        else:
            ic50 = np_rng.lognormal(mean=6.5, sigma=1.5) * (1.0 + neg_weight * 0.3)
            ic50 = max(100.0, min(50000.0, ic50))

        if ic50 < 50:
            binding_level = BindingLevel.STRONG_BINDER
            percent_rank = max(0.01, np_rng.lognormal(mean=-1.0, sigma=1.0))
        elif ic50 < 500:
            binding_level = BindingLevel.WEAK_BINDER
            percent_rank = np_rng.lognormal(mean=0.5, sigma=0.8)
        elif ic50 < 5000:
            binding_level = BindingLevel.WEAK_BINDER
            percent_rank = np_rng.lognormal(mean=1.5, sigma=0.6)
        else:
            binding_level = BindingLevel.NON_BINDER
            percent_rank = np_rng.lognormal(mean=3.0, sigma=0.5)

        percent_rank = max(0.01, min(50.0, percent_rank))

        presentation_score = (
            max(0.0, min(0.99, 1.0 - np_rng.lognormal(mean=0.0, sigma=0.5) * ic50 / 50000.0))
            if label_positive
            else max(0.0, min(0.3, 0.1 + rng.random() * 0.2))
        )

        half_life = predict_stability_half_life(mut_peptide, ic50)
        stability_rank = percent_rank * (0.8 + rng.random() * 0.4)

        predictions.append(BindingPrediction(
            peptide=mut_peptide,
            hla_allele=hla,
            predictor="mhcflurry",
            ic50=ic50,
            percent_rank=percent_rank,
            binding_level=binding_level,
            presentation_score=presentation_score,
            half_life=half_life,
            stability_rank=stability_rank,
        ))

    candidate.binding_predictions = predictions
    candidate.best_binding = min(
        predictions,
        key=lambda p: (p.percent_rank if p.percent_rank is not None else 999.0,
                       p.ic50 if p.ic50 is not None else 999999.0),
    )

    # Expression with log-normal distribution
    if label_positive:
        gene_tpm = np_rng.lognormal(mean=2.5, sigma=1.0)
        mutant_expr = np_rng.lognormal(mean=1.5, sigma=0.8)
        ccf = min(1.0, max(0.3, np_rng.beta(a=8, b=3)))
        clonal = ccf > 0.7 or rng.random() < 0.6
    else:
        gene_tpm = np_rng.lognormal(mean=0.5, sigma=1.5)
        mutant_expr = max(0.0, np_rng.lognormal(mean=-0.5, sigma=1.2) - 0.5)
        ccf = max(0.01, min(0.6, np_rng.beta(a=2, b=8)))
        clonal = ccf > 0.5 and rng.random() < 0.2

    vaf = ccf * (0.4 + rng.random() * 0.2)
    rna_vaf = vaf * (0.7 + rng.random() * 0.3)

    candidate.expression = ExpressionData(
        gene_symbol=gene_symbol,
        gene_tpm=round(gene_tpm, 2),
        transcript_tpm=round(max(0, gene_tpm * (0.8 + rng.random() * 0.4)), 2),
        vaf=round(vaf, 4),
        rna_vaf=round(rna_vaf, 4),
        mutant_expression=round(mutant_expr, 4),
    )

    candidate.clonal_status = ClonalStatus.CLONAL if clonal else ClonalStatus.SUBCLONAL
    candidate.cancer_cell_fraction = round(ccf, 4)

    # Antigen processing using PSSM models
    cleavage_score = predict_proteasome_cleavage(mut_peptide)
    tap_score = predict_tap_transport(mut_peptide)
    b2m_mut = rng.random() < 0.015  # ~1.5% B2M mutation rate

    candidate.processing_features = AntigenProcessingFeatures(
        proteasome_cleavage_score=round(cleavage_score, 4),
        tap_transport_score=round(tap_score, 4),
        b2m_mutated=b2m_mut,
    )

    # Self-similarity and agretopicity (computed by FeatureEngineer)
    candidate.self_similarity = None
    candidate.agretopicity = None

    # MS evidence
    ms_detected = rng.random() < (0.35 if label_positive else 0.03)
    candidate.ms_evidence = MSEvidence(
        peptide=mut_peptide,
        detected=ms_detected,
        pepquery_score=round(np_rng.uniform(0.75, 0.99), 4) if ms_detected else None,
        spectral_count=rng.randint(1, 8) if ms_detected else 0,
    )

    # HLA LOH (rare event)
    candidate.presenting_allele_loh = rng.random() < 0.03
    candidate.is_driver_gene = is_driver

    return candidate


def generate_synthetic_training_set(
    n_samples: int = 400,
    positive_fraction: float = 0.2,
    random_seed: int = 42,
) -> tuple[pd.DataFrame, np.ndarray]:
    """Build a labeled feature matrix with biologically realistic distributions.

    The synthetic data generator models:
      - Log-normal IC50 distributions (realistic for binding affinities)
      - Beta-distributed cancer cell fractions (clonal vs subclonal)
      - Log-normal gene expression (realistic for RNA-seq TPM)
      - HLA allele-specific binding preferences
      - Feature correlations (strong binding + high expression = more likely MS detected)
      - Realistic variant type frequencies (SNV > indel > frameshift > fusion > splice)

    Returns:
        Tuple of (feature_matrix_with_metadata, binary_labels).
    """
    rng = random.Random(random_seed)
    np_rng = np.random.default_rng(random_seed)
    engineer = FeatureEngineer()
    candidates: list[NeoantigenCandidate] = []

    n_positive = int(n_samples * positive_fraction)
    labels: list[int] = [1] * n_positive + [0] * (n_samples - n_positive)
    rng.shuffle(labels)

    for index, label in enumerate(labels):
        sample_id = f"S{index % 20:03d}"
        candidates.append(_make_candidate(sample_id, bool(label), rng, np_rng))

    candidates = engineer.enrich_candidates(candidates)
    matrix = engineer.build_feature_matrix(candidates)
    matrix.insert(0, "sample_id", [candidate.sample_id for candidate in candidates])
    matrix.insert(0, "candidate_id", [candidate.candidate_id for candidate in candidates])
    matrix["label"] = labels
    return matrix, np.array(labels, dtype=int)


def load_training_table(path: str) -> tuple[pd.DataFrame, np.ndarray]:
    table = pd.read_csv(path, sep="\t")
    if "label" not in table.columns:
        raise ValueError("Training table must include a 'label' column")
    labels = table["label"].astype(int).to_numpy()
    features = table.drop(columns=["label"])
    for column in ("candidate_id", "sample_id"):
        if column in features.columns:
            features = features.drop(columns=[column])
    missing = [column for column in FEATURE_COLUMNS if column not in features.columns]
    if missing:
        raise ValueError(f"Training table missing feature columns: {missing}")
    return features[FEATURE_COLUMNS], labels


def train_demo_model(
    output_path: str,
    n_samples: int = 400,
    random_seed: int = 42,
    training_table: Optional[str] = None,
) -> dict[str, float]:
    from sklearn.model_selection import train_test_split

    from neoflux.ml_ranker import NeoFluxRanker

    if training_table:
        features, labels = load_training_table(training_table)
    else:
        table, labels = generate_synthetic_training_set(
            n_samples=n_samples,
            random_seed=random_seed,
        )
        features = table[FEATURE_COLUMNS]

    x_train, x_val, y_train, y_val = train_test_split(
        features,
        labels,
        test_size=0.25,
        random_state=random_seed,
        stratify=labels,
    )

    ranker = NeoFluxRanker()
    metrics = ranker.train(x_train, y_train, x_val, y_val)
    ranker.save(output_path)
    return metrics


def train_publication_model(
    training_table: str,
    output_path: str,
    metrics_path: Optional[str] = None,
    test_size: float = 0.2,
    random_seed: int = 42,
) -> dict[str, float]:
    """Train the publication ensemble on a labeled table and evaluate it.

    Performs a stratified held-out split, trains the XGBoost/LightGBM/logistic
    ensemble, and reports held-out AUROC for the ensemble alongside a
    binding-affinity-only baseline. The two ROC curves are compared on the
    identical test set with a DeLong test (Sun & Xu 2014).

    Args:
        training_table: TSV with FEATURE_COLUMNS + a 'label' column.
        output_path: where to pickle the trained ensemble.
        metrics_path: optional JSON path for the evaluation metrics.
        test_size: held-out fraction (0.2 -> ~20% test).
        random_seed: split/seed reproducibility.

    Returns:
        Dict of evaluation metrics (counts, AUROCs, DeLong z/p-value).
    """
    import json

    from sklearn.impute import SimpleImputer
    from sklearn.linear_model import LogisticRegression
    from sklearn.metrics import roc_auc_score
    from sklearn.model_selection import train_test_split

    from neoflux.evaluation import delong_roc_test
    from neoflux.ml_ranker import NeoFluxRanker

    features, labels = load_training_table(training_table)

    x_train, x_test, y_train, y_test = train_test_split(
        features,
        labels,
        test_size=test_size,
        random_state=random_seed,
        stratify=labels,
    )

    ranker = NeoFluxRanker()
    ranker.train(x_train, y_train, x_test, y_test)
    ranker.save(output_path)

    # Full-ensemble held-out scores.
    ensemble_scores = ranker.predict(x_test)
    ensemble_auroc = float(roc_auc_score(y_test, ensemble_scores))

    # Binding-affinity-only baseline: a logistic model that sees nothing but
    # MHC binding affinity + percentile rank. This isolates the predictive
    # value of binding alone, the standard first-pass neoantigen filter.
    binding_cols = ["best_ic50", "best_percent_rank"]
    baseline_imputer = SimpleImputer(strategy="median")
    xb_train = baseline_imputer.fit_transform(x_train[binding_cols])
    xb_test = baseline_imputer.transform(x_test[binding_cols])
    baseline = LogisticRegression(max_iter=1000, class_weight="balanced")
    baseline.fit(xb_train, y_train)
    baseline_scores = baseline.predict_proba(xb_test)[:, 1]
    baseline_auroc = float(roc_auc_score(y_test, baseline_scores))

    delong = delong_roc_test(
        np.asarray(y_test), ensemble_scores, baseline_scores
    )

    metrics = {
        "n_total": int(len(labels)),
        "n_train": int(len(y_train)),
        "n_test": int(len(y_test)),
        "n_positive": int(labels.sum()),
        "n_negative": int((1 - labels).sum()),
        "ensemble_auroc": round(ensemble_auroc, 4),
        "baseline_auroc": round(baseline_auroc, 4),
        "delong_z": round(delong["z"], 4),
        "delong_p_value": delong["p_value"],
    }

    if metrics_path:
        with open(metrics_path, "w") as handle:
            json.dump(metrics, handle, indent=2)

    return metrics
