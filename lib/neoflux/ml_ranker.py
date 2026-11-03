"""
ML-based neoantigen ranking model for NeoFlux.

Implements a gradient-boosted ensemble (XGBoost + LightGBM + logistic
regression baseline) for predicting neoantigen immunogenicity. The model
outputs a calibrated probability score in [0, 1] representing the
likelihood that a candidate is truly immunogenic (presented AND
recognized by T cells).

Training data sources:
    - TESLA consortium validated neoantigens (ground truth)
    - IEDB T-cell epitope assays (positive/negative)
    - HLA Ligand Atlas immunopeptidomics (presentation ground truth)
"""

from __future__ import annotations

import logging
import pickle
from pathlib import Path
from typing import Optional

import numpy as np
import pandas as pd

from neoflux.features import FEATURE_COLUMNS, FeatureEngineer
from neoflux.models import NeoantigenCandidate, VariantType

logger = logging.getLogger("neoflux")


def heuristic_score_candidates(
    candidates: list[NeoantigenCandidate],
) -> list[NeoantigenCandidate]:
    """
    Rank candidates with an adaptive literature-weighted multi-factor model.

    Weights are derived from published neoantigen prioritization studies:
      - Binding: 30% (TESLA and Wells et al. 2020; Yadav et al. 2014)
      - Expression: 20% (McGranahan et al. 2016; Rizvi et al. 2015)
      - Clonality: 15% (McGranahan et al. 2016; Gejman et al. 2022)
      - Processing: 10% (proteasome cleavage and TAP transport contributions;
        Fritsch et al. 2014; Peters et al. 2003)
      - Sequence: 10% (foreignness and self-similarity; Chowell et al. 2015)
      - MS evidence: 15% (empirical presentation validation; Yadav et al. 2014)

    When data for a category is entirely missing (e.g. no expression data),
    its weight is redistributed proportionally among available categories.
    This prevents missing data from artificially capping the score.

    Immune escape penalties (multiplicative):
      - HLA LOH: 0.5x (McGranahan et al. 2017; Rosenthal et al. 2019)
      - B2M mutation: 0.5x (Zaretsky et al. 2016)
      - Both present: 0.3x
    """
    engineer = FeatureEngineer()
    candidates = engineer.enrich_candidates(candidates)

    BASE_WEIGHTS = {
        "binding": 0.30,
        "expression": 0.20,
        "clonality": 0.15,
        "processing": 0.10,
        "sequence": 0.10,
        "ms": 0.15,
    }

    for candidate in candidates:
        # -- Binding score (0-1) --
        binding_score = 0.2
        if candidate.best_binding:
            if candidate.best_binding.percent_rank is not None:
                binding_score = max(0.0, min(1.0, 1.0 - candidate.best_binding.percent_rank / 5.0))
            elif candidate.best_binding.ic50 is not None:
                binding_score = max(0.0, min(1.0, 1.0 - candidate.best_binding.ic50 / 5000.0))
            if candidate.best_binding.presentation_score is not None:
                pres = candidate.best_binding.presentation_score
                binding_score = 0.6 * binding_score + 0.4 * pres

        # -- Expression score (0-1) --
        expression_score = 0.2
        if candidate.expression and candidate.expression.mutant_expression is not None:
            me = candidate.expression.mutant_expression
            expression_score = max(0.0, min(1.0, np.log1p(me) / np.log1p(50.0)))
        elif candidate.expression and candidate.expression.gene_tpm is not None:
            gt = candidate.expression.gene_tpm
            expression_score = max(0.0, min(1.0, np.log1p(gt) / np.log1p(100.0)))

        # -- Clonality score (0-1) --
        clonality_score = 0.5
        if candidate.cancer_cell_fraction is not None:
            clonality_score = max(0.0, min(1.0, candidate.cancer_cell_fraction))

        # -- Processing score (0-1) --
        processing_score = 0.3
        if candidate.processing_features:
            ps = candidate.processing_features.proteasome_cleavage_score
            ts = candidate.processing_features.tap_transport_score
            if ps is not None and ts is not None:
                processing_score = 0.5 * ps + 0.5 * ts
            elif ps is not None:
                processing_score = ps
            elif ts is not None:
                processing_score = ts

        # -- Sequence score (0-1) --
        sequence_score = 0.5
        if candidate.self_similarity is not None:
            sequence_score = max(0.0, min(1.0, 1.0 - candidate.self_similarity))
        if candidate.agretopicity is not None and candidate.agretopicity > 1.0:
            sequence_score = min(1.0, sequence_score * (1.0 + np.log(candidate.agretopicity) * 0.3))

        # -- MS evidence score (0-1) --
        ms_score = 0.0
        if candidate.ms_evidence and candidate.ms_evidence.detected:
            ms_score = 0.8
            if candidate.ms_evidence.pepquery_score is not None:
                ms_score = 0.5 + 0.5 * candidate.ms_evidence.pepquery_score
            sc = candidate.ms_evidence.spectral_count
            if sc is not None and sc > 1:
                ms_score = min(1.0, ms_score + 0.1 * min(candidate.ms_evidence.spectral_count, 5))

        # -- Adaptive re-weighting: detect available categories --
        scores = {
            "binding": binding_score,
            "expression": expression_score,
            "clonality": clonality_score,
            "processing": processing_score,
            "sequence": sequence_score,
            "ms": ms_score,
        }
        available = {}
        if candidate.best_binding and candidate.best_binding.ic50 is not None:
            available["binding"] = True
        has_expr = (
            candidate.expression is not None
            and (
                candidate.expression.mutant_expression is not None
                or candidate.expression.gene_tpm is not None
            )
        )
        if has_expr:
            available["expression"] = True
        if candidate.cancer_cell_fraction is not None:
            available["clonality"] = True
        has_proc = (
            candidate.processing_features is not None
            and candidate.processing_features.proteasome_cleavage_score is not None
        )
        if has_proc:
            available["processing"] = True
        if candidate.self_similarity is not None:
            available["sequence"] = True
        if candidate.ms_evidence and candidate.ms_evidence.detected:
            available["ms"] = True

        if available:
            total_weight = sum(BASE_WEIGHTS[k] for k in available)
            raw_score = sum(BASE_WEIGHTS[k] / total_weight * scores[k] for k in available)
        else:
            raw_score = sum(BASE_WEIGHTS[k] * scores[k] for k in scores)

        # Immune escape penalties (multiplicative)
        immune_penalty = 1.0
        has_loh = candidate.presenting_allele_loh
        has_b2m = (candidate.processing_features is not None
                   and candidate.processing_features.b2m_mutated)

        if has_loh and has_b2m:
            immune_penalty = 0.3
        elif has_loh:
            immune_penalty = 0.5
        elif has_b2m:
            immune_penalty = 0.5

        # Foreignness boost for non-canonical peptide sources
        # Literature evidence:
        #   - Yang et al. 2019 (Nat Med): fusion neoantigens stimulate T cell responses
        #   - Laumont et al. 2018 (Sci Transl Med): ~90% of TSAs from noncoding regions
        #   - Genome Biology 2023: SV-derived neoantigens have lower self-similarity
        #   - Frameshift INDEL neoepitopes: ~21% immunogenic vs ~3% for SNVs
        foreignness_boost = 1.0
        if not candidate.wildtype_peptide:
            foreignness_boost = 1.15  # fusion, nCORF — entirely novel, no central tolerance
        elif candidate.peptide_source in (
            VariantType.FRAMESHIFT, VariantType.FUSION, VariantType.SPLICE
        ):
            foreignness_boost = 1.11  # partially novel sequences

        candidate.neoflux_score = float(
            max(0.0, min(1.0, raw_score * immune_penalty * foreignness_boost))
        )

    candidates.sort(key=lambda cand: cand.neoflux_score or 0.0, reverse=True)
    for index, candidate in enumerate(candidates, start=1):
        candidate.neoflux_rank = index
    return candidates

# Default hyperparameters for ensemble members
DEFAULT_XGBOOST_PARAMS = {
    "objective": "binary:logistic",
    "eval_metric": "aucpr",
    "max_depth": 6,
    "learning_rate": 0.05,
    "n_estimators": 500,
    "subsample": 0.8,
    "colsample_bytree": 0.8,
    "min_child_weight": 5,
    "scale_pos_weight": 10,
    "tree_method": "hist",
    "random_state": 42,
}

DEFAULT_LGBM_PARAMS = {
    "objective": "binary",
    "metric": "average_precision",
    "max_depth": 6,
    "learning_rate": 0.05,
    "n_estimators": 500,
    "subsample": 0.8,
    "colsample_bytree": 0.8,
    "min_child_samples": 10,
    "is_unbalance": True,
    "random_state": 42,
    "verbose": -1,
}


class NeoFluxRanker:
    """
    Ensemble ML model for neoantigen prioritization.

    Architecture:
        1. XGBoost gradient-boosted trees
        2. LightGBM gradient-boosted trees
        3. Logistic regression baseline

    The final score is a weighted average of the three models, calibrated
    via Platt scaling (sigmoid) to produce a probability in [0, 1].

    XGBoost and LightGBM are chosen because they:
        - Handle missing features natively (critical for optional MS data)
        - Provide feature importance and SHAP compatibility
        - Are computationally efficient for the feature set size
        - Handle class imbalance well

    The logistic regression serves as an interpretable baseline and
    regularizes the ensemble.
    """

    def __init__(
        self,
        xgb_params: Optional[dict] = None,
        lgbm_params: Optional[dict] = None,
        ensemble_weights: tuple[float, float, float] = (0.4, 0.4, 0.2),
    ):
        """
        Args:
            xgb_params: XGBoost hyperparameters (defaults used if None).
            lgbm_params: LightGBM hyperparameters (defaults used if None).
            ensemble_weights: Weights for (XGBoost, LightGBM, LogReg).
        """
        self.xgb_params = xgb_params or DEFAULT_XGBOOST_PARAMS
        self.lgbm_params = lgbm_params or DEFAULT_LGBM_PARAMS
        self.ensemble_weights = ensemble_weights
        self.feature_engineer = FeatureEngineer()

        # Models (initialized during training)
        self.xgb_model = None
        self.lgbm_model = None
        self.logreg_model = None

        # Calibration (Platt scaling parameters)
        self.calibration_a: float = 1.0
        self.calibration_b: float = 0.0

        self.is_trained = False

    def train(
        self,
        X: pd.DataFrame,
        y: np.ndarray,
        X_val: Optional[pd.DataFrame] = None,
        y_val: Optional[np.ndarray] = None,
    ) -> dict[str, float]:
        """
        Train the ensemble model on labeled neoantigen data.

        Args:
            X: Feature matrix (n_samples, n_features). Columns must match
               FEATURE_COLUMNS order.
            y: Binary labels (1 = immunogenic, 0 = not immunogenic).
            X_val: Optional validation feature matrix.
            y_val: Optional validation labels.

        Returns:
            Dict of training metrics (AUROC, AUPRC per model and ensemble).

        NOTE: Requires xgboost, lightgbm, and scikit-learn to be installed.
        These are listed as optional dependencies in pyproject.toml.
        """
        try:
            import lightgbm as lgb
            import xgboost as xgb
            from sklearn.impute import SimpleImputer
            from sklearn.linear_model import LogisticRegression
            from sklearn.metrics import average_precision_score, roc_auc_score
        except ImportError as e:
            raise ImportError(
                f"ML training requires xgboost, lightgbm, and scikit-learn: {e}"
            )

        logger.info(
            "Training NeoFlux ranker on %d samples (%d positive, %d negative)",
            len(y), int(y.sum()), int((1 - y).sum()),
        )

        metrics = {}

        # 1. XGBoost
        self.xgb_model = xgb.XGBClassifier(**self.xgb_params)
        eval_set = [(X_val, y_val)] if X_val is not None else None
        self.xgb_model.fit(
            X, y,
            eval_set=eval_set,
            verbose=False,
        )

        # 2. LightGBM
        self.lgbm_model = lgb.LGBMClassifier(**self.lgbm_params)
        self.lgbm_model.fit(
            X, y,
            eval_set=[(X_val, y_val)] if X_val is not None else None,
        )

        # 3. Logistic regression (impute missing for LR)
        imputer = SimpleImputer(strategy="median")
        X_imputed = imputer.fit_transform(X)
        self.logreg_model = LogisticRegression(
            max_iter=1000,
            class_weight="balanced",
            random_state=42,
        )
        self.logreg_model.fit(X_imputed, y)
        self._imputer = imputer

        # Validation metrics
        if X_val is not None and y_val is not None:
            preds_xgb = self.xgb_model.predict_proba(X_val)[:, 1]
            preds_lgbm = self.lgbm_model.predict_proba(X_val)[:, 1]
            X_val_imp = imputer.transform(X_val)
            preds_lr = self.logreg_model.predict_proba(X_val_imp)[:, 1]

            # Ensemble prediction
            preds_ensemble = (
                self.ensemble_weights[0] * preds_xgb
                + self.ensemble_weights[1] * preds_lgbm
                + self.ensemble_weights[2] * preds_lr
            )

            for name, preds in [
                ("xgboost", preds_xgb),
                ("lightgbm", preds_lgbm),
                ("logreg", preds_lr),
                ("ensemble", preds_ensemble),
            ]:
                metrics[f"{name}_auroc"] = roc_auc_score(y_val, preds)
                metrics[f"{name}_auprc"] = average_precision_score(y_val, preds)
                logger.info(
                    "  %s: AUROC=%.4f AUPRC=%.4f",
                    name, metrics[f"{name}_auroc"], metrics[f"{name}_auprc"],
                )

            # Platt scaling calibration on validation set
            self._calibrate(preds_ensemble, y_val)

        self.is_trained = True
        return metrics

    def predict(self, X: pd.DataFrame) -> np.ndarray:
        """
        Predict neoantigen immunogenicity scores.

        Args:
            X: Feature matrix (n_candidates, n_features).

        Returns:
            Array of calibrated scores in [0, 1].
        """
        if not self.is_trained:
            raise RuntimeError("Model has not been trained. Call train() first.")

        preds_xgb = self.xgb_model.predict_proba(X)[:, 1]
        preds_lgbm = self.lgbm_model.predict_proba(X)[:, 1]
        X_imp = self._imputer.transform(X)
        preds_lr = self.logreg_model.predict_proba(X_imp)[:, 1]

        raw_ensemble = (
            self.ensemble_weights[0] * preds_xgb
            + self.ensemble_weights[1] * preds_lgbm
            + self.ensemble_weights[2] * preds_lr
        )

        # Apply Platt calibration
        calibrated = 1.0 / (
            1.0 + np.exp(self.calibration_a * raw_ensemble + self.calibration_b)
        )
        return calibrated

    def predict_candidates(
        self, candidates: list[NeoantigenCandidate]
    ) -> list[NeoantigenCandidate]:
        """
        Score a list of NeoantigenCandidate objects and attach scores + ranks.

        This is the primary inference entry point. It:
        1. Enriches candidates with derived features
        2. Builds the feature matrix
        3. Runs the ensemble model
        4. Attaches scores and ranks to each candidate
        5. Computes SHAP explanations

        Args:
            candidates: List of NeoantigenCandidate objects.

        Returns:
            The same list, sorted by score (descending), with neoflux_score,
            neoflux_rank, and shap_top_features populated.
        """
        # Enrich derived features
        candidates = self.feature_engineer.enrich_candidates(candidates)

        # Build feature matrix
        X = self.feature_engineer.build_feature_matrix(candidates)

        # Predict
        scores = self.predict(X)

        # Immune-escape gating (applied regardless of scorer): if the presenting
        # allele is lost (HLA LOH) or B2M is disrupted, the peptide cannot reach
        # the cell surface, so the immunogenicity score is penalized downstream
        # of the model. Mirrors the multipliers used by the heuristic scorer.
        for i, candidate in enumerate(candidates):
            has_loh = candidate.presenting_allele_loh
            has_b2m = (
                candidate.processing_features is not None
                and candidate.processing_features.b2m_mutated
            )
            if has_loh and has_b2m:
                scores[i] = scores[i] * 0.3
            elif has_loh or has_b2m:
                scores[i] = scores[i] * 0.5

        # Compute SHAP values
        shap_explanations = self._compute_shap(X)

        # Attach scores and rank
        ranked_indices = np.argsort(-scores)
        for rank_pos, idx in enumerate(ranked_indices):
            candidates[idx].neoflux_score = float(scores[idx])
            candidates[idx].neoflux_rank = rank_pos + 1
            if shap_explanations is not None:
                candidates[idx].shap_top_features = shap_explanations[idx]

        # Sort by score descending
        candidates.sort(key=lambda c: c.neoflux_score or 0, reverse=True)

        logger.info(
            "Scored %d candidates. Top score: %.4f, Median: %.4f",
            len(candidates),
            scores.max() if len(scores) > 0 else 0,
            float(np.median(scores)) if len(scores) > 0 else 0,
        )
        return candidates

    def _compute_shap(
        self, X: pd.DataFrame, top_k: int = 3
    ) -> Optional[list[dict[str, float]]]:
        """
        Compute SHAP feature importance for each candidate.

        Returns top_k most influential features per candidate as a dict
        of {feature_name: shap_value}.

        NOTE: Uses TreeSHAP for XGBoost (fast, exact for trees). The
        SHAP values from XGBoost are used as the primary explanation
        since it receives the highest ensemble weight.
        """
        try:
            import shap
        except ImportError:
            logger.warning(
                "shap package not installed; skipping SHAP explanations"
            )
            return None

        explainer = shap.TreeExplainer(self.xgb_model)
        shap_values = explainer.shap_values(X)

        explanations = []
        feature_names = list(X.columns)

        for i in range(len(X)):
            abs_shap = np.abs(shap_values[i])
            top_indices = np.argsort(-abs_shap)[:top_k]
            top_features = {
                feature_names[j]: float(shap_values[i][j])
                for j in top_indices
            }
            explanations.append(top_features)

        return explanations

    def _calibrate(
        self, raw_scores: np.ndarray, y_true: np.ndarray
    ) -> None:
        """Fit Platt scaling (logistic calibration) on validation predictions."""
        try:
            from sklearn.linear_model import LogisticRegression as LR
        except ImportError:
            return

        calibrator = LR()
        calibrator.fit(raw_scores.reshape(-1, 1), y_true)
        self.calibration_a = float(-calibrator.coef_[0][0])
        self.calibration_b = float(-calibrator.intercept_[0])
        logger.info(
            "Platt calibration: a=%.4f, b=%.4f",
            self.calibration_a, self.calibration_b,
        )

    def save(self, path: str | Path) -> None:
        """
        Save the trained model to disk.

        Saves all ensemble members, calibration parameters, and metadata
        as a pickle file.

        Args:
            path: Output file path (typically .pkl).
        """
        path = Path(path)
        state = {
            "xgb_model": self.xgb_model,
            "lgbm_model": self.lgbm_model,
            "logreg_model": self.logreg_model,
            "imputer": self._imputer if hasattr(self, "_imputer") else None,
            "calibration_a": self.calibration_a,
            "calibration_b": self.calibration_b,
            "ensemble_weights": self.ensemble_weights,
            "feature_columns": FEATURE_COLUMNS,
            "is_trained": self.is_trained,
        }
        with open(path, "wb") as f:
            pickle.dump(state, f)
        logger.info("Saved NeoFlux ranker to %s", path)

    @classmethod
    def load(cls, path: str | Path) -> NeoFluxRanker:
        """
        Load a trained model from disk.

        Args:
            path: Path to the saved .pkl file.

        Returns:
            A NeoFluxRanker instance ready for inference.
        """
        path = Path(path)
        with open(path, "rb") as f:
            state = pickle.load(f)

        ranker = cls(ensemble_weights=state["ensemble_weights"])
        ranker.xgb_model = state["xgb_model"]
        ranker.lgbm_model = state["lgbm_model"]
        ranker.logreg_model = state["logreg_model"]
        ranker._imputer = state["imputer"]
        ranker.calibration_a = state["calibration_a"]
        ranker.calibration_b = state["calibration_b"]
        ranker.is_trained = state["is_trained"]

        logger.info("Loaded NeoFlux ranker from %s", path)
        return ranker

    def get_feature_importance(self) -> pd.DataFrame:
        """
        Get feature importance from the XGBoost model (gain-based).

        Returns:
            DataFrame with columns [feature, importance], sorted descending.
        """
        if self.xgb_model is None:
            raise RuntimeError("Model not trained.")

        importance = self.xgb_model.feature_importances_
        df = pd.DataFrame({
            "feature": FEATURE_COLUMNS,
            "importance": importance,
        }).sort_values("importance", ascending=False)
        return df
