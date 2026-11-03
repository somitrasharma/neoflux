"""
Peptide sequence generation from somatic variants.

Generates mutant and wildtype peptide sequences for all supported variant
types: SNVs, indels/frameshifts, gene fusions, splice variants, and
non-canonical ORFs. Each generator produces peptide windows of configurable
length (default 8-11 for MHC-I, 15 for MHC-II).
"""

from __future__ import annotations

import logging
from dataclasses import dataclass
from typing import Optional

from neoflux.models import (
    GenomicVariant,
    NeoantigenCandidate,
    VariantType,
)
from neoflux.utils import generate_candidate_id, validate_peptide

logger = logging.getLogger("neoflux")

# Default peptide length range for MHC-I (8-11mers) and MHC-II (15mers)
DEFAULT_MHC_I_LENGTHS = (8, 9, 10, 11)
DEFAULT_MHC_II_LENGTHS = (15,)


@dataclass
class PeptideWindow:
    """A peptide window extracted from a protein sequence around a mutation."""
    peptide: str
    wildtype_peptide: Optional[str]
    start_pos: int
    end_pos: int
    mutation_position_in_peptide: int
    source: VariantType


class PeptideGenerator:
    """
    Generates candidate neoantigen peptides from somatic variants.

    This class implements the core logic for extracting peptide windows
    around somatic mutations in protein sequences. It handles all variant
    types defined in the NeoFlux pipeline.
    """

    def __init__(
        self,
        mhc_i_lengths: tuple[int, ...] = DEFAULT_MHC_I_LENGTHS,
        mhc_ii_lengths: tuple[int, ...] = DEFAULT_MHC_II_LENGTHS,
        flanking_context: int = 30,
    ):
        """
        Args:
            mhc_i_lengths: Peptide lengths to generate for MHC class I.
            mhc_ii_lengths: Peptide lengths to generate for MHC class II.
            flanking_context: Number of flanking amino acids to include
                around the mutation site when extracting protein context.
        """
        self.mhc_i_lengths = mhc_i_lengths
        self.mhc_ii_lengths = mhc_ii_lengths
        self.flanking_context = flanking_context
        self.all_lengths = sorted(set(mhc_i_lengths + mhc_ii_lengths))

    def generate_from_snv(
        self,
        variant: GenomicVariant,
        mutant_protein: str,
        wildtype_protein: str,
        mutation_pos: int,
    ) -> list[PeptideWindow]:
        """
        Generate peptide windows around a single amino acid substitution (SNV).

        For each peptide length k, generates all k-mer windows that contain
        the mutated residue, paired with the corresponding wildtype k-mer.

        Args:
            variant: The source GenomicVariant.
            mutant_protein: Full mutant protein sequence.
            wildtype_protein: Full wildtype protein sequence.
            mutation_pos: 0-based position of the mutated residue in the protein.

        Returns:
            List of PeptideWindow objects.

        NOTE: In the full pipeline, mutant_protein and wildtype_protein are
        derived from VEP annotations + reference FASTA via pyensembl or
        pyfaidx. This function receives the pre-computed sequences.
        """
        windows = []

        for k in self.all_lengths:
            # Sliding window: all k-mers containing the mutation site
            start_min = max(0, mutation_pos - k + 1)
            start_max = min(mutation_pos, len(mutant_protein) - k)

            for start in range(start_min, start_max + 1):
                end = start + k
                mt_pep = mutant_protein[start:end]
                wt_pep = wildtype_protein[start:end]

                if not validate_peptide(mt_pep):
                    continue
                if mt_pep == wt_pep:
                    continue  # No actual amino acid change in this window

                mut_pos_in_pep = mutation_pos - start

                windows.append(PeptideWindow(
                    peptide=mt_pep,
                    wildtype_peptide=wt_pep,
                    start_pos=start,
                    end_pos=end,
                    mutation_position_in_peptide=mut_pos_in_pep,
                    source=VariantType.SNV,
                ))

        logger.debug(
            "Generated %d peptide windows for SNV at %s pos %d",
            len(windows), variant.gene_symbol, mutation_pos,
        )
        return windows

    def generate_from_indel(
        self,
        variant: GenomicVariant,
        mutant_protein: str,
        wildtype_protein: str,
        indel_start: int,
    ) -> list[PeptideWindow]:
        """
        Generate peptide windows for insertions, deletions, and frameshifts.

        For frameshifts, the mutant protein sequence diverges from wildtype
        at the indel site and continues with a novel sequence until the next
        stop codon. All k-mers spanning or downstream of the indel are
        generated.

        Args:
            variant: The source GenomicVariant.
            mutant_protein: Mutant protein (may be truncated at new stop).
            wildtype_protein: Wildtype protein.
            indel_start: 0-based position where the indel begins.

        Returns:
            List of PeptideWindow objects.
        """
        windows = []
        is_frameshift = variant.variant_type == VariantType.FRAMESHIFT

        for k in self.all_lengths:
            if is_frameshift:
                # For frameshifts, generate all k-mers from indel_start onward
                max_start = min(
                    len(mutant_protein) - k,
                    indel_start + self.flanking_context,
                )
                start_range = range(max(0, indel_start - k + 1), max_start + 1)
            else:
                # For in-frame indels, use the same windowing as SNVs
                start_min = max(0, indel_start - k + 1)
                start_max = min(indel_start, len(mutant_protein) - k)
                start_range = range(start_min, start_max + 1)

            for start in start_range:
                end = start + k
                mt_pep = mutant_protein[start:end]

                # For frameshifts, wildtype may not have a corresponding region
                wt_pep = None
                if start + k <= len(wildtype_protein):
                    wt_pep = wildtype_protein[start:end]

                if not validate_peptide(mt_pep):
                    continue
                if wt_pep and mt_pep == wt_pep:
                    continue

                windows.append(PeptideWindow(
                    peptide=mt_pep,
                    wildtype_peptide=wt_pep,
                    start_pos=start,
                    end_pos=end,
                    mutation_position_in_peptide=max(0, indel_start - start),
                    source=variant.variant_type,
                ))

        logger.debug(
            "Generated %d peptide windows for %s at %s pos %d",
            len(windows), variant.variant_type.value,
            variant.gene_symbol, indel_start,
        )
        return windows

    def generate_from_fusion(
        self,
        variant: GenomicVariant,
        fusion_protein: str,
        junction_pos: int,
    ) -> list[PeptideWindow]:
        """
        Generate peptide windows spanning a gene fusion junction.

        The junction peptides are the novel sequences not present in either
        parent protein. Only k-mers that span the actual junction point
        are generated (these are the neoantigen-producing windows).

        Args:
            variant: The source GenomicVariant (with fusion_partner fields).
            fusion_protein: Predicted fusion protein sequence.
            junction_pos: 0-based position of the junction in fusion_protein.

        Returns:
            List of PeptideWindow objects (no wildtype, since fusions are novel).
        """
        windows = []

        for k in self.all_lengths:
            # Only generate k-mers that actually span the junction
            start_min = max(0, junction_pos - k + 1)
            start_max = min(junction_pos, len(fusion_protein) - k)

            for start in range(start_min, start_max + 1):
                end = start + k
                pep = fusion_protein[start:end]

                if not validate_peptide(pep):
                    continue

                windows.append(PeptideWindow(
                    peptide=pep,
                    wildtype_peptide=None,  # No WT for fusions
                    start_pos=start,
                    end_pos=end,
                    mutation_position_in_peptide=junction_pos - start,
                    source=VariantType.FUSION,
                ))

        logger.debug(
            "Generated %d junction peptides for fusion %s--%s",
            len(windows),
            variant.fusion_partner_5p,
            variant.fusion_partner_3p,
        )
        return windows

    def generate_from_splice_variant(
        self,
        variant: GenomicVariant,
        aberrant_protein: str,
        reference_protein: str,
        divergence_pos: int,
    ) -> list[PeptideWindow]:
        """
        Generate peptide windows from aberrant splicing events.

        Aberrant splicing (exon skipping, intron retention, alternative
        splice site usage) can produce novel protein sequences at the
        splice junction. This method generates k-mers around the point
        where the aberrant and reference proteins diverge.

        Args:
            variant: The source GenomicVariant (with splice_event_type).
            aberrant_protein: Protein from the aberrantly spliced transcript.
            reference_protein: Protein from the canonical transcript.
            divergence_pos: 0-based position where sequences diverge.

        Returns:
            List of PeptideWindow objects.

        NOTE: In the full pipeline, aberrant_protein is derived by
        translating the rMATS/SUPPA2-predicted transcript using the
        modified transcript FASTA.
        """
        windows = []

        for k in self.all_lengths:
            start_min = max(0, divergence_pos - k + 1)
            start_max = min(
                divergence_pos + self.flanking_context,
                len(aberrant_protein) - k,
            )

            for start in range(start_min, start_max + 1):
                end = start + k
                mt_pep = aberrant_protein[start:end]

                wt_pep = None
                if start + k <= len(reference_protein):
                    wt_pep = reference_protein[start:end]

                if not validate_peptide(mt_pep):
                    continue
                if wt_pep and mt_pep == wt_pep:
                    continue

                windows.append(PeptideWindow(
                    peptide=mt_pep,
                    wildtype_peptide=wt_pep,
                    start_pos=start,
                    end_pos=end,
                    mutation_position_in_peptide=max(0, divergence_pos - start),
                    source=VariantType.SPLICE,
                ))

        logger.debug(
            "Generated %d peptide windows for splice variant in %s (%s)",
            len(windows), variant.gene_symbol, variant.splice_event_type,
        )
        return windows

    def generate_from_noncanonical_orf(
        self,
        sample_id: str,
        orf_id: str,
        orf_protein: str,
        orf_type: str = "novel",
    ) -> list[PeptideWindow]:
        """
        Generate peptide windows from non-canonical open reading frames.

        Non-canonical ORFs include:
        - Upstream ORFs (uORFs) in 5'UTRs
        - Downstream ORFs (dORFs) in 3'UTRs
        - ORFs in lncRNAs
        - ORFs from endogenous retroelements
        - Intergenic ORFs

        These are entirely novel sequences with no wildtype counterpart
        in the canonical proteome.

        Args:
            sample_id: Patient/sample identifier.
            orf_id: Identifier for the non-canonical ORF.
            orf_protein: Translated protein sequence from the ncORF.
            orf_type: Type of non-canonical ORF (for annotation).

        Returns:
            List of PeptideWindow objects.

        NOTE: In the full pipeline, ncORFs are identified by scanning
        RNA-seq data with tools like Ribo-TISH, ORFquant, or a custom
        3-frame translation + ORF finder. The input here is the
        pre-translated sequence.
        """
        windows = []

        for k in self.all_lengths:
            for start in range(len(orf_protein) - k + 1):
                end = start + k
                pep = orf_protein[start:end]

                if not validate_peptide(pep):
                    continue

                windows.append(PeptideWindow(
                    peptide=pep,
                    wildtype_peptide=None,
                    start_pos=start,
                    end_pos=end,
                    mutation_position_in_peptide=0,
                    source=VariantType.NON_CANONICAL_ORF,
                ))

        logger.debug(
            "Generated %d peptide windows from ncORF %s (%s)",
            len(windows), orf_id, orf_type,
        )
        return windows

    def windows_to_candidates(
        self,
        windows: list[PeptideWindow],
        variant: GenomicVariant,
        sample_id: str,
        hla_allele: str = "",
    ) -> list[NeoantigenCandidate]:
        """
        Convert PeptideWindow objects into NeoantigenCandidate models.

        This is a convenience method for downstream pipeline stages. The
        returned candidates have binding, expression, and ML fields unset
        (to be filled by later modules).

        Args:
            windows: List of PeptideWindow from any generator method.
            variant: Source GenomicVariant.
            sample_id: Patient/sample identifier.
            hla_allele: HLA allele string (for candidate ID generation).

        Returns:
            List of NeoantigenCandidate models.
        """
        candidates = []
        for win in windows:
            cid = generate_candidate_id(
                sample_id=sample_id,
                variant_id=variant.variant_id,
                peptide=win.peptide,
                hla_allele=hla_allele,
            )
            candidate = NeoantigenCandidate(
                candidate_id=cid,
                sample_id=sample_id,
                variant=variant,
                mutant_peptide=win.peptide,
                wildtype_peptide=win.wildtype_peptide,
                peptide_length=len(win.peptide),
                peptide_source=win.source,
            )
            candidates.append(candidate)
        return candidates
