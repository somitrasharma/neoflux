"""
Cohort-level analysis and aggregation for NeoFlux.

Computes summary statistics and comparative metrics across a cohort
of patients, suitable for SITC-format reporting and benchmarking.
"""

from __future__ import annotations

import logging
from collections import Counter

import numpy as np
import pandas as pd

from neoflux.models import (
    CohortSummary,
    PatientSample,
    VariantType,
)

logger = logging.getLogger("neoflux")


class CohortAnalyzer:
    """
    Aggregates neoantigen results across a cohort of patients.

    Produces summary statistics and data structures suitable for
    generating SITC-format figures and tables.
    """

    def __init__(self, cohort_name: str = "NeoFlux_Cohort"):
        self.cohort_name = cohort_name
        self.patients: list[PatientSample] = []

    def add_patient(self, patient: PatientSample) -> None:
        """Add a patient's results to the cohort."""
        self.patients.append(patient)

    def compute_summary(self) -> CohortSummary:
        """
        Compute cohort-level summary statistics.

        Returns:
            CohortSummary object with aggregate metrics.
        """
        all_candidates = []
        variant_type_counter: Counter = Counter()
        ms_detected_top = 0
        top_total = 0
        hla_loh_patients = 0

        for patient in self.patients:
            all_candidates.extend(patient.candidates)

            for cand in patient.candidates:
                variant_type_counter[cand.peptide_source.value] += 1

            # MS detection rate among top 20 candidates per patient
            top_candidates = sorted(
                patient.candidates,
                key=lambda c: c.neoflux_score or 0,
                reverse=True,
            )[:20]
            for cand in top_candidates:
                top_total += 1
                if cand.ms_evidence and cand.ms_evidence.detected:
                    ms_detected_top += 1

            # HLA LOH prevalence
            if any(a.is_loh for a in patient.hla_alleles):
                hla_loh_patients += 1

        candidates_per_patient = [
            len(p.candidates) for p in self.patients
        ]

        summary = CohortSummary(
            cohort_name=self.cohort_name,
            n_patients=len(self.patients),
            samples=[p.sample_id for p in self.patients],
            total_candidates=len(all_candidates),
            median_candidates_per_patient=(
                float(np.median(candidates_per_patient))
                if candidates_per_patient else None
            ),
            variant_type_counts=dict(variant_type_counter),
            ms_detection_rate=(
                ms_detected_top / max(top_total, 1) if top_total > 0 else None
            ),
            hla_loh_prevalence=(
                hla_loh_patients / max(len(self.patients), 1)
            ),
        )

        logger.info(
            "Cohort summary: %d patients, %d total candidates, "
            "median %.0f/patient",
            summary.n_patients,
            summary.total_candidates,
            summary.median_candidates_per_patient or 0,
        )
        return summary

    def build_candidate_dataframe(self) -> pd.DataFrame:
        """
        Build a cohort-wide DataFrame of all candidates with output columns.

        This DataFrame is the primary input for cohort-level plotting
        and statistical analysis.

        Returns:
            DataFrame with one row per candidate across all patients.
        """
        rows = []
        for patient in self.patients:
            for cand in patient.candidates:
                rows.append(cand.to_output_row())

        df = pd.DataFrame(rows)
        logger.info(
            "Built cohort DataFrame: %d candidates x %d columns",
            len(df), len(df.columns),
        )
        return df

    def compute_neoantigen_burden(self) -> pd.DataFrame:
        """
        Compute neoantigen burden per patient, stratified by variant type.

        Returns:
            DataFrame with columns: sample_id, total, SNV, indel,
            frameshift, fusion, splice, non_canonical_orf
        """
        rows = []
        for patient in self.patients:
            counts = Counter(
                cand.peptide_source.value for cand in patient.candidates
            )
            row = {"sample_id": patient.sample_id, "total": len(patient.candidates)}
            for vt in VariantType:
                row[vt.value] = counts.get(vt.value, 0)
            rows.append(row)

        return pd.DataFrame(rows)

    def compute_score_vs_ms_correlation(self) -> pd.DataFrame:
        """
        Compute correlation between NeoFlux scores and MS detection.

        For SITC figure: do higher-scoring candidates show higher
        MS detection rates?

        Returns:
            DataFrame with columns: score_bin, n_candidates, n_detected,
            detection_rate
        """
        all_candidates = []
        for patient in self.patients:
            all_candidates.extend(patient.candidates)

        if not all_candidates:
            return pd.DataFrame()

        scores = [c.neoflux_score or 0 for c in all_candidates]
        detected = [
            1 if (c.ms_evidence and c.ms_evidence.detected) else 0
            for c in all_candidates
        ]

        # Bin scores into deciles
        bins = np.linspace(0, 1, 11)
        bin_labels = [f"{bins[i]:.1f}-{bins[i+1]:.1f}" for i in range(len(bins)-1)]

        score_bins = np.digitize(scores, bins) - 1
        score_bins = np.clip(score_bins, 0, len(bin_labels) - 1)

        rows = []
        for i, label in enumerate(bin_labels):
            mask = score_bins == i
            n_total = int(mask.sum())
            n_det = int(np.array(detected)[mask].sum()) if n_total > 0 else 0
            rows.append({
                "score_bin": label,
                "n_candidates": n_total,
                "n_detected": n_det,
                "detection_rate": n_det / max(n_total, 1),
            })

        return pd.DataFrame(rows)

    def compute_hla_loh_impact(self) -> pd.DataFrame:
        """
        Analyze the impact of HLA LOH on neoantigen presentation.

        For SITC figure: how many top candidates are lost due to LOH
        of their presenting allele?

        Returns:
            DataFrame with columns: sample_id, total_candidates,
            loh_affected, fraction_lost
        """
        rows = []
        for patient in self.patients:
            total = len(patient.candidates)
            loh_affected = sum(
                1 for c in patient.candidates if c.presenting_allele_loh
            )
            rows.append({
                "sample_id": patient.sample_id,
                "total_candidates": total,
                "loh_affected": loh_affected,
                "fraction_lost": loh_affected / max(total, 1),
            })

        return pd.DataFrame(rows)
