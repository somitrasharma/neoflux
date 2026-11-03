"""
Immunopeptidomics evidence scoring for NeoFlux.

Handles integration of LC-MS/MS data into the neoantigen prioritization
pipeline. Parses FragPipe/MSFragger search results and PepQuery2 targeted
validation output to produce per-peptide evidence scores.
"""

from __future__ import annotations

import csv
import logging
from pathlib import Path

from neoflux.models import MSEvidence, NeoantigenCandidate

logger = logging.getLogger("neoflux")

# Default FDR threshold for accepting MS identifications
DEFAULT_FDR_THRESHOLD = 0.01

# PepQuery2 p-value threshold for confident identification
DEFAULT_PEPQUERY_PVALUE_THRESHOLD = 0.05


class MSEvidenceScorer:
    """
    Scores neoantigen candidates based on immunopeptidomics evidence.

    Workflow:
    1. Parse FragPipe/MSFragger search results (PSM-level)
    2. Parse PepQuery2 targeted validation results
    3. For each candidate peptide, determine if it was detected
    4. Assign confidence scores based on PSM quality and PepQuery2 metrics

    This module is designed to be optional: if no MS data is provided,
    all MSEvidence fields are left as None/default, and the ML model
    handles missing features natively.
    """

    def __init__(
        self,
        fdr_threshold: float = DEFAULT_FDR_THRESHOLD,
        pepquery_pvalue_threshold: float = DEFAULT_PEPQUERY_PVALUE_THRESHOLD,
    ):
        self.fdr_threshold = fdr_threshold
        self.pepquery_pvalue_threshold = pepquery_pvalue_threshold

        # Internal peptide lookup tables
        self._fragpipe_peptides: dict[str, dict] = {}
        self._pepquery_results: dict[str, dict] = {}

    def load_fragpipe_results(self, path: str | Path) -> int:
        """
        Parse FragPipe/MSFragger peptide-level search results.

        Expects a TSV file with at minimum:
            Peptide, SpectralCount, Probability, (optional) BestScore

        This is the output of the FragPipe combined_peptide.tsv or
        equivalent MSFragger output.

        Args:
            path: Path to FragPipe peptide results TSV.

        Returns:
            Number of peptides loaded.

        NOTE: The exact column names depend on the FragPipe version.
        This parser is flexible and logs warnings for missing columns.
        Adjust column mappings as needed for your FragPipe configuration.
        """
        path = Path(path)
        if not path.exists():
            logger.warning("FragPipe results not found: %s", path)
            return 0

        count = 0
        with open(path, newline="") as f:
            reader = csv.DictReader(f, delimiter="\t")
            for row in reader:
                peptide = row.get("Peptide", row.get("peptide", "")).strip()
                # Remove flanking residues if present (e.g., K.PEPTIDE.R)
                if "." in peptide:
                    parts = peptide.split(".")
                    if len(parts) == 3:
                        peptide = parts[1]

                # Remove modification annotations (e.g., M[147])
                clean_peptide = ""
                skip = False
                for ch in peptide:
                    if ch == "[":
                        skip = True
                    elif ch == "]":
                        skip = False
                    elif not skip and ch.isalpha():
                        clean_peptide += ch.upper()
                peptide = clean_peptide

                if not peptide:
                    continue

                spectral_count = int(
                    row.get("SpectralCount", row.get("Spectral Count", 0))
                )
                probability = float(
                    row.get("Probability", row.get("MaxPepProb", 1.0))
                )
                best_score = float(
                    row.get("BestScore", row.get("Best Score", 0.0))
                )

                self._fragpipe_peptides[peptide] = {
                    "spectral_count": spectral_count,
                    "probability": probability,
                    "best_score": best_score,
                }
                count += 1

        logger.info(
            "Loaded %d peptides from FragPipe results: %s", count, path
        )
        return count

    def load_pepquery_results(self, path: str | Path) -> int:
        """
        Parse PepQuery2 targeted validation results.

        Expects a TSV file with columns:
            peptide, pvalue, score, n_db_hits, confident

        This is the standard PepQuery2 psm_rank.txt output.

        Args:
            path: Path to PepQuery2 results TSV.

        Returns:
            Number of results loaded.

        NOTE: PepQuery2 is invoked as a separate Nextflow process
        (modules/immunopeptidomics/pepquery2.nf). This method parses
        its output for downstream scoring.
        """
        path = Path(path)
        if not path.exists():
            logger.warning("PepQuery2 results not found: %s", path)
            return 0

        count = 0
        with open(path, newline="") as f:
            reader = csv.DictReader(f, delimiter="\t")
            for row in reader:
                peptide = row.get("peptide", "").strip().upper()
                if not peptide:
                    continue

                pvalue = float(row.get("pvalue", 1.0))
                score = float(row.get("score", 0.0))
                confident = row.get("confident", "No").strip().lower() == "yes"

                self._pepquery_results[peptide] = {
                    "pvalue": pvalue,
                    "score": score,
                    "confident": confident,
                }
                count += 1

        logger.info(
            "Loaded %d PepQuery2 results from: %s", count, path
        )
        return count

    def score_peptide(self, peptide: str) -> MSEvidence:
        """
        Score a single peptide based on loaded MS evidence.

        Combines FragPipe search results with PepQuery2 validation
        to produce a final MSEvidence object.

        Args:
            peptide: Amino acid sequence to look up.

        Returns:
            MSEvidence object with detection status and scores.
        """
        peptide = peptide.upper().strip()
        evidence = MSEvidence(peptide=peptide)

        # Check FragPipe results
        fp_data = self._fragpipe_peptides.get(peptide)
        if fp_data:
            evidence.spectral_count = fp_data["spectral_count"]
            evidence.best_search_engine_score = fp_data["best_score"]

        # Check PepQuery2 results
        pq_data = self._pepquery_results.get(peptide)
        if pq_data:
            evidence.pepquery_pvalue = pq_data["pvalue"]
            evidence.pepquery_score = pq_data["score"]

            # Mark as detected if PepQuery2 confirms with confidence
            if (
                pq_data["confident"]
                and pq_data["pvalue"] < self.pepquery_pvalue_threshold
            ):
                evidence.detected = True

        # If no PepQuery2 but FragPipe found it with good spectral evidence
        if not evidence.detected and fp_data:
            if fp_data["spectral_count"] >= 2 and fp_data["probability"] > 0.9:
                evidence.detected = True

        return evidence

    def score_candidates(
        self, candidates: list[NeoantigenCandidate]
    ) -> list[NeoantigenCandidate]:
        """
        Score all candidates and attach MS evidence.

        Args:
            candidates: List of NeoantigenCandidate objects (modified in place).

        Returns:
            The same list with ms_evidence populated.
        """
        detected_count = 0
        for cand in candidates:
            cand.ms_evidence = self.score_peptide(cand.mutant_peptide)
            if cand.ms_evidence.detected:
                detected_count += 1

        logger.info(
            "MS evidence scoring: %d/%d candidates detected (%.1f%%)",
            detected_count, len(candidates),
            100.0 * detected_count / max(len(candidates), 1),
        )
        return candidates

    @property
    def has_data(self) -> bool:
        """Whether any MS data has been loaded."""
        return bool(self._fragpipe_peptides or self._pepquery_results)
