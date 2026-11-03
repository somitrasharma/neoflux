import numpy as np
import pandas as pd
import pytest

from neoflux.features import FEATURE_COLUMNS
from neoflux.ml_ranker import NeoFluxRanker, heuristic_score_candidates
from neoflux.models import (
    BindingLevel,
    BindingPrediction,
    GenomicVariant,
    NeoantigenCandidate,
    VariantType,
)
from neoflux.training_data import generate_synthetic_training_set, train_demo_model


def _candidate(score_hint: float) -> NeoantigenCandidate:
    variant = GenomicVariant(
        sample_id="S1",
        chromosome="chr1",
        position=100,
        ref_allele="A",
        alt_allele="T",
        variant_type=VariantType.SNV,
    )
    candidate = NeoantigenCandidate(
        candidate_id=f"cand-{score_hint}",
        sample_id="S1",
        variant=variant,
        mutant_peptide="GLYVGLG" if score_hint > 0.5 else "AAAAAAAAA",
        peptide_length=9,
        peptide_source=VariantType.SNV,
    )
    candidate.best_binding = BindingPrediction(
        peptide=candidate.mutant_peptide,
        hla_allele="HLA-A*02:01",
        predictor="mhcflurry",
        ic50=100.0 if score_hint > 0.5 else 3000.0,
        percent_rank=0.5 if score_hint > 0.5 else 20.0,
        binding_level=BindingLevel.STRONG_BINDER if score_hint > 0.5 else BindingLevel.NON_BINDER,
        presentation_score=score_hint,
    )
    return candidate


def test_heuristic_ranker_orders_higher_binding_first():
    low = _candidate(0.1)
    high = _candidate(0.9)
    ranked = heuristic_score_candidates([low, high])
    assert ranked[0].candidate_id == high.candidate_id
    assert ranked[0].neoflux_rank == 1


def test_train_and_predict_demo_model(tmp_path):
    pytest.importorskip("xgboost")
    pytest.importorskip("lightgbm")
    pytest.importorskip("sklearn")

    model_path = tmp_path / "demo_model.pkl"
    metrics = train_demo_model(str(model_path), n_samples=120, random_seed=7)
    assert model_path.exists()
    assert "ensemble_auroc" in metrics or "xgboost_auroc" in metrics

    ranker = NeoFluxRanker.load(model_path)
    table, labels = generate_synthetic_training_set(n_samples=40, random_seed=7)
    features = table[FEATURE_COLUMNS]
    scores = ranker.predict(features)
    assert len(scores) == len(labels)
    assert np.all((scores >= 0.0) & (scores <= 1.0))

    candidates = [_candidate(0.9), _candidate(0.1)]
    scored = ranker.predict_candidates(candidates)
    assert scored[0].neoflux_score is not None
    assert scored[0].neoflux_rank == 1


def test_feature_importance_after_training(tmp_path):
    pytest.importorskip("xgboost")
    pytest.importorskip("lightgbm")
    pytest.importorskip("sklearn")

    model_path = tmp_path / "demo_model.pkl"
    train_demo_model(str(model_path), n_samples=80, random_seed=3)
    ranker = NeoFluxRanker.load(model_path)
    importance = ranker.get_feature_importance()
    assert isinstance(importance, pd.DataFrame)
    assert "feature" in importance.columns
    assert len(importance) == len(FEATURE_COLUMNS)
