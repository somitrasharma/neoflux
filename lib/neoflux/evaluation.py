"""
Model evaluation utilities for the NeoFlux ranker.

Provides held-out performance metrics and a DeLong test for comparing two
correlated ROC curves (e.g. the full ensemble vs. a binding-affinity-only
baseline evaluated on the same test set). The DeLong implementation follows
the fast algorithm of Sun & Xu (2014), IEEE Signal Processing Letters.
"""

from __future__ import annotations

import numpy as np
from scipy import stats


def _compute_midrank(x: np.ndarray) -> np.ndarray:
    """Mid-rank transform used by the fast DeLong covariance estimator."""
    order = np.argsort(x)
    x_sorted = x[order]
    n = len(x)
    cumulative = np.zeros(n, dtype=float)
    i = 0
    while i < n:
        j = i
        while j < n and x_sorted[j] == x_sorted[i]:
            j += 1
        cumulative[i:j] = 0.5 * (i + j - 1) + 1
        i = j
    out = np.empty(n, dtype=float)
    out[order] = cumulative
    return out


def _fast_delong(predictions_sorted: np.ndarray, n_positive: int):
    """Fast DeLong AUC + covariance for a stack of prediction vectors.

    Args:
        predictions_sorted: array (k_models, n_samples) with all positive
            samples placed before all negative samples along axis 1.
        n_positive: number of positive samples (the leading columns).

    Returns:
        (aucs, covariance) where aucs is length-k and covariance is (k, k).
    """
    m = n_positive
    n = predictions_sorted.shape[1] - m
    positive = predictions_sorted[:, :m]
    negative = predictions_sorted[:, m:]
    k = predictions_sorted.shape[0]

    tx = np.empty([k, m], dtype=float)
    ty = np.empty([k, n], dtype=float)
    tz = np.empty([k, m + n], dtype=float)
    for r in range(k):
        tx[r, :] = _compute_midrank(positive[r, :])
        ty[r, :] = _compute_midrank(negative[r, :])
        tz[r, :] = _compute_midrank(predictions_sorted[r, :])

    aucs = tz[:, :m].sum(axis=1) / m / n - (m + 1.0) / 2.0 / n
    v01 = (tz[:, :m] - tx[:, :]) / n
    v10 = 1.0 - (tz[:, m:] - ty[:, :]) / m
    sx = np.cov(v01)
    sy = np.cov(v10)
    delongcov = sx / m + sy / n
    return aucs, np.atleast_2d(delongcov)


def delong_roc_test(
    y_true: np.ndarray,
    scores_a: np.ndarray,
    scores_b: np.ndarray,
) -> dict[str, float]:
    """Compare two correlated ROC AUCs on the same labels via the DeLong test.

    Args:
        y_true: binary labels (1 = positive).
        scores_a: prediction scores for model A (e.g. the ensemble).
        scores_b: prediction scores for model B (e.g. binding-only baseline).

    Returns:
        Dict with auc_a, auc_b, z statistic, and two-sided p-value.
    """
    y_true = np.asarray(y_true, dtype=int)
    order = np.argsort(-y_true, kind="mergesort")  # positives (label 1) first
    n_positive = int((y_true == 1).sum())

    predictions_sorted = np.vstack(
        (np.asarray(scores_a, dtype=float)[order],
         np.asarray(scores_b, dtype=float)[order])
    )
    aucs, cov = _fast_delong(predictions_sorted, n_positive)

    var = cov[0, 0] + cov[1, 1] - 2.0 * cov[0, 1]
    if var <= 0:
        var = np.finfo(float).tiny
    z = (aucs[0] - aucs[1]) / np.sqrt(var)
    p = 2.0 * stats.norm.sf(abs(z))
    return {
        "auc_a": float(aucs[0]),
        "auc_b": float(aucs[1]),
        "z": float(z),
        "p_value": float(p),
    }
