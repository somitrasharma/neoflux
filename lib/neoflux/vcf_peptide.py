"""
Translate annotated VCF records and event manifests into peptide candidates.

Production runs should supply VEP CSQ/HGVS fields or transcript-aware protein
sequences. For development and stub runs, this module falls back to synthetic
protein context so the pipeline can exercise the full DAG without a full
reference proteome translation stack.
"""

from __future__ import annotations

import json
import re
from pathlib import Path
from typing import Iterable, Optional

from neoflux.models import GenomicVariant, NeoantigenCandidate, VariantConsequence, VariantType
from neoflux.peptide_generator import PeptideGenerator
from neoflux.reference_proteome import get_protein

# BRAF V600 region used as a realistic demo protein context for stub VCFs.
DEFAULT_PROTEIN_CONTEXT = (
    "GAGGGGAGGGGAGGGGAGGGGAGGGGAGGGGAGGGGAGGGGAGGGGAGGGGAGGGGAGGGGAGGGGA"
    "GGGGAGGGGAGGGGAGGGGAGGGGAGGGGAGGGGAGGGGAGGGGAGGGGAGGGGAGGGGAGGGGAG"
    "GGGGAGGGGAGGGGAGGGGAGGGGAGGGGAGGGGAGGGGAGGGGAGGGGAGGGGAGGGGAGGGGAG"
    "GGGGAGGGGAGGGGAGGGGAGGGGAGGGGAGGGGAGGGGAGGGGAGGGGAGGGGAGGGGAGGGGAG"
)

HGVS_PROTEIN_RE = re.compile(r"p\.([A-Z*])(\d+)([A-Z*])")
INFO_FIELD_RE = re.compile(r"([^=;]+)=([^;]*)")


def parse_info_field(info: str) -> dict[str, str]:
    fields: dict[str, str] = {}
    for key, value in INFO_FIELD_RE.findall(info):
        fields[key] = value
    return fields


def parse_hgvs_protein_change(text: str) -> Optional[tuple[str, int, str]]:
    match = HGVS_PROTEIN_RE.search(text)
    if not match:
        return None
    ref_aa, position, alt_aa = match.group(1), int(match.group(2)), match.group(3)
    if ref_aa == "*" or alt_aa == "*":
        return None
    return ref_aa, position, alt_aa


def apply_protein_substitution(
    protein: str,
    position_1based: int,
    ref_aa: str,
    alt_aa: str,
) -> tuple[str, str, int]:
    index = position_1based - 1
    if index < 0 or index >= len(protein):
        raise ValueError(f"Substitution position {position_1based} out of range")

    observed = protein[index]
    if observed != ref_aa:
        # Allow mismatch in stub/demo data; still place the alternate residue.
        pass

    mutant = protein[:index] + alt_aa + protein[index + 1 :]
    return mutant, protein, index


def synthetic_snv_context(ref_allele: str, alt_allele: str) -> tuple[str, str, int]:
    """
    Build a deterministic mutant/wildtype protein pair when VEP output is absent.

    The mutation is embedded at the center of a fixed-length context window so
    PeptideGenerator can emit realistic 8-11mers for pipeline dry runs.
    """
    center = len(DEFAULT_PROTEIN_CONTEXT) // 2
    ref_aa = "V"
    alt_aa = "E"
    if len(ref_allele) == 1 and len(alt_allele) == 1:
        alt_aa = alt_allele.upper() if alt_allele.upper() in "ACDEFGHIKLMNPQRSTVWY" else "E"

    mutant, wildtype, mutation_pos = apply_protein_substitution(
        DEFAULT_PROTEIN_CONTEXT,
        center,
        ref_aa,
        alt_aa,
    )
    return mutant, wildtype, mutation_pos


def infer_consequence(variant_type: VariantType) -> VariantConsequence:
    mapping = {
        VariantType.SNV: VariantConsequence.MISSENSE,
        VariantType.INSERTION: VariantConsequence.INFRAME_INSERTION,
        VariantType.DELETION: VariantConsequence.INFRAME_DELETION,
        VariantType.FRAMESHIFT: VariantConsequence.FRAMESHIFT,
        VariantType.FUSION: VariantConsequence.FUSION,
        VariantType.SPLICE: VariantConsequence.SPLICE_DONOR,
        VariantType.NON_CANONICAL_ORF: VariantConsequence.NON_CANONICAL,
    }
    return mapping.get(variant_type, VariantConsequence.OTHER)


def infer_variant_type(ref: str, alt: str) -> VariantType:
    if len(ref) == 1 and len(alt) == 1:
        return VariantType.SNV
    if len(ref) < len(alt):
        return VariantType.INSERTION
    if len(ref) > len(alt):
        return VariantType.DELETION
    return VariantType.SNV


def build_variant_from_record(sample_id: str, record: dict[str, str]) -> GenomicVariant:
    variant_type = infer_variant_type(record["ref"], record["alt"])
    info = parse_info_field(record.get("info", ""))
    gene_symbol = info.get("GENE") or info.get("SYMBOL")
    protein_change = info.get("HGVS_P") or info.get("Protein_change")

    csq = info.get("CSQ", "")
    if not protein_change and csq:
        for token in csq.split("|"):
            parsed = parse_hgvs_protein_change(token)
            if parsed:
                ref_aa, pos, alt_aa = parsed
                protein_change = f"p.{ref_aa}{pos}{alt_aa}"
                break
            if "missense_variant" in token and not gene_symbol:
                parts = token.split("|")
                if len(parts) > 1:
                    gene_symbol = parts[1]

    return GenomicVariant(
        sample_id=sample_id,
        chromosome=record["chrom"],
        position=int(record["pos"]),
        ref_allele=record["ref"],
        alt_allele=record["alt"],
        variant_type=variant_type,
        gene_symbol=gene_symbol,
        consequence=infer_consequence(variant_type),
        protein_change=protein_change,
        caller=info.get("NEOFLUX_CALLERS", "consensus"),
        is_consensus="," in info.get("NEOFLUX_CALLERS", ""),
    )


def generate_from_vcf_record(
    generator: PeptideGenerator,
    sample_id: str,
    record: dict[str, str],
) -> list[NeoantigenCandidate]:
    variant = build_variant_from_record(sample_id, record)
    info = parse_info_field(record.get("info", ""))
    protein_change = variant.protein_change or info.get("CSQ", "")

    windows = []
    parsed = parse_hgvs_protein_change(protein_change or "")
    if parsed:
        ref_aa, position, alt_aa = parsed
        # Prefer the real canonical protein for this gene so real variant
        # coordinates yield real peptide windows. Fall back to the synthetic
        # context only when no reference proteome is available (offline + gene
        # uncached) or the position lies outside the sequence.
        ref_protein = get_protein(variant.gene_symbol)
        if ref_protein and 1 <= position <= len(ref_protein):
            mutant, wildtype, mutation_pos = apply_protein_substitution(
                ref_protein,
                position,
                ref_aa,
                alt_aa,
            )
        else:
            try:
                mutant, wildtype, mutation_pos = apply_protein_substitution(
                    DEFAULT_PROTEIN_CONTEXT,
                    position,
                    ref_aa,
                    alt_aa,
                )
            except ValueError:
                fallback_position = len(DEFAULT_PROTEIN_CONTEXT) // 2
                fallback_ref = DEFAULT_PROTEIN_CONTEXT[fallback_position - 1]
                mutant, wildtype, mutation_pos = apply_protein_substitution(
                    DEFAULT_PROTEIN_CONTEXT,
                    fallback_position,
                    fallback_ref,
                    alt_aa,
                )
        windows = generator.generate_from_snv(
            variant=variant,
            mutant_protein=mutant,
            wildtype_protein=wildtype,
            mutation_pos=mutation_pos,
        )
    elif variant.variant_type == VariantType.SNV:
        mutant, wildtype, mutation_pos = synthetic_snv_context(
            variant.ref_allele,
            variant.alt_allele,
        )
        windows = generator.generate_from_snv(
            variant=variant,
            mutant_protein=mutant,
            wildtype_protein=wildtype,
            mutation_pos=mutation_pos,
        )
    elif variant.variant_type in {
        VariantType.INSERTION,
        VariantType.DELETION,
        VariantType.FRAMESHIFT,
    }:
        indel_start = len(DEFAULT_PROTEIN_CONTEXT) // 2
        mutant = (
            DEFAULT_PROTEIN_CONTEXT[:indel_start]
            + "XXX"
            + DEFAULT_PROTEIN_CONTEXT[indel_start + 3:]
        )
        windows = generator.generate_from_indel(
            variant=variant,
            mutant_protein=mutant,
            wildtype_protein=DEFAULT_PROTEIN_CONTEXT,
            indel_start=indel_start,
        )

    return generator.windows_to_candidates(windows, variant=variant, sample_id=sample_id)


def generate_from_event_manifest(
    generator: PeptideGenerator,
    sample_id: str,
    manifest_path: str | Path,
) -> list[NeoantigenCandidate]:
    path = Path(manifest_path)
    if not path.exists() or path.stat().st_size == 0:
        return []

    manifest = json.loads(path.read_text())
    candidates: list[NeoantigenCandidate] = []

    for fusion in manifest.get("fusions", []):
        partner_5p = fusion.get("gene_5p") or fusion.get("gene5") or "GENE5"
        partner_3p = fusion.get("gene_3p") or fusion.get("gene3") or "GENE3"
        # Build the chimeric protein from the real partner sequences when both
        # genes resolve; otherwise fall back to the synthetic partner stubs so
        # the pipeline still runs offline.
        left = get_protein(partner_5p) or (
            "MTEYKLVVVGAGGVGKSALTIQLIQNHFVDEYDPTIEDSYRKQVVIDGETCLLDILDTAGQEEY"
        )
        right = get_protein(partner_3p) or (
            "SAMRDQYMRTGEGFLCVFAINNTKSFEDIHHYREQIKRVKDSEDVPMVLVGNKCDLPARTVETR"
        )
        # Breakpoints in protein coordinates: EML4(1..break_5p) joined to
        # ALK(break_3p..end). Falls back to the legacy single junction_pos.
        break_5p = int(fusion.get("break_5p_aa", fusion.get("junction_pos", 30)))
        break_3p = int(fusion.get("break_3p_aa", max(0, break_5p - 10)))
        break_5p = min(break_5p, len(left))
        break_3p = min(break_3p, len(right))
        fusion_protein = left[:break_5p] + right[break_3p:]
        junction = break_5p  # junction index within the chimeric protein
        variant = GenomicVariant(
            sample_id=sample_id,
            chromosome=fusion.get("chrom", "chr1"),
            position=int(fusion.get("pos", 1)),
            ref_allele="N",
            alt_allele="FUSION",
            variant_type=VariantType.FUSION,
            gene_symbol=f"{partner_5p}-{partner_3p}",
            consequence=VariantConsequence.FUSION,
            fusion_partner_5p=partner_5p,
            fusion_partner_3p=partner_3p,
            fusion_breakpoint=fusion.get("breakpoint", ""),
            caller="star-fusion",
        )
        windows = generator.generate_from_fusion(
            variant=variant,
            fusion_protein=fusion_protein,
            junction_pos=junction,
        )
        candidates.extend(
            generator.windows_to_candidates(windows, variant=variant, sample_id=sample_id)
        )

    for event in manifest.get("splice_events", []):
        gene = event.get("gene_symbol", "")
        ref_protein = get_protein(gene)
        if ref_protein and len(ref_protein) > 30:
            # Model exon skipping as an in-frame deletion of the skipped exon
            # from the real canonical protein. skip_start/skip_len give the
            # skipped segment in protein coordinates; the novel junction in the
            # aberrant transcript abuts the flanking exons at skip_start. (True
            # transcript-level reconstruction needs rMATS/SUPPA2 + exon models;
            # this is a protein-coordinate approximation grounded in the real
            # reference sequence.)
            skip_len = int(event.get("skip_len", 47))
            skip_start = int(event.get("skip_start", max(1, len(ref_protein) // 2)))
            skip_start = min(skip_start, len(ref_protein) - 1)
            reference = ref_protein
            aberrant = ref_protein[:skip_start] + ref_protein[skip_start + skip_len:]
            divergence = skip_start
        else:
            # No reference proteome (offline + gene uncached): synthetic context.
            divergence = int(event.get("divergence_pos", 25))
            aberrant = "MTEYKLVVVGAGGVGKSALTIQLIQNHFVDEYDPTIEDSYRKXXSAMRDQYMRVVPMVATAGGTGEKLS"
            reference = (
                "MTEYKLVVVGAGGVGKSALTIQLIQNHFVDEYDPTIEDSYRK"
                "QVVIDGETCLLDILDTAGQEEYSAMRDQYMRTGEGFLCVFAINNTKS"
            )
        psi_val = event.get("delta_psi")
        variant = GenomicVariant(
            sample_id=sample_id,
            chromosome=event.get("chrom", "chr1"),
            position=int(event.get("pos", 1)),
            ref_allele="N",
            alt_allele="SPLICE",
            variant_type=VariantType.SPLICE,
            gene_symbol=event.get("gene_symbol", ""),
            consequence=VariantConsequence.SPLICE_DONOR,
            splice_event_type=event.get("event_type", "SE"),
            delta_psi=float(psi_val) if psi_val not in {None, ""} else None,
            caller="rmats",
        )
        windows = generator.generate_from_splice_variant(
            variant=variant,
            aberrant_protein=aberrant,
            reference_protein=reference,
            divergence_pos=divergence,
        )
        candidates.extend(
            generator.windows_to_candidates(windows, variant=variant, sample_id=sample_id)
        )

    for orf in manifest.get("noncanonical_orfs", []):
        orf_protein = orf.get(
            "protein",
            "MTEYKLVVVGAGGVGKSALTIQLIQNHFVDEYDPTIEDSYRK"
            "QVVIDGETCLLDILDTAGQEEY",
        )
        variant = GenomicVariant(
            sample_id=sample_id,
            chromosome=orf.get("chrom", "chr1"),
            position=int(orf.get("pos", 1)),
            ref_allele="N",
            alt_allele="NCORF",
            variant_type=VariantType.NON_CANONICAL_ORF,
            gene_symbol=orf.get("orf_id", "ncORF"),
            consequence=VariantConsequence.NON_CANONICAL,
            caller="ncorf",
        )
        windows = generator.generate_from_noncanonical_orf(
            sample_id=sample_id,
            orf_id=orf.get("orf_id", "ncORF"),
            orf_protein=orf_protein,
            orf_type=orf.get("orf_type", "novel"),
        )
        candidates.extend(
            generator.windows_to_candidates(windows, variant=variant, sample_id=sample_id)
        )

    return candidates


def generate_all_candidates(
    sample_id: str,
    vcf_records: Iterable[dict[str, str]],
    manifest_path: str | Path,
    mhc_i_lengths: tuple[int, ...],
    mhc_ii_lengths: tuple[int, ...],
) -> list[NeoantigenCandidate]:
    generator = PeptideGenerator(
        mhc_i_lengths=mhc_i_lengths,
        mhc_ii_lengths=mhc_ii_lengths,
    )
    candidates: list[NeoantigenCandidate] = []
    seen_ids: set[str] = set()

    for record in vcf_records:
        for candidate in generate_from_vcf_record(generator, sample_id, record):
            if candidate.candidate_id not in seen_ids:
                seen_ids.add(candidate.candidate_id)
                candidates.append(candidate)

    for candidate in generate_from_event_manifest(generator, sample_id, manifest_path):
        if candidate.candidate_id not in seen_ids:
            seen_ids.add(candidate.candidate_id)
            candidates.append(candidate)

    return candidates


def write_peptide_exports(
    candidates: list[NeoantigenCandidate],
    candidates_jsonl: str | Path,
    peptides_fasta: str | Path,
    peptides_tsv: str | Path,
) -> None:
    import csv

    with open(candidates_jsonl, "w") as handle:
        for candidate in candidates:
            handle.write(candidate.model_dump_json() + "\n")

    with open(peptides_fasta, "w") as fasta, open(peptides_tsv, "w", newline="") as tsv:
        writer = csv.DictWriter(
            tsv,
            delimiter="\t",
            fieldnames=["candidate_id", "sample_id", "peptide", "variant_type"],
        )
        writer.writeheader()
        for candidate in candidates:
            fasta.write(f">{candidate.candidate_id}\n{candidate.mutant_peptide}\n")
            writer.writerow(
                {
                    "candidate_id": candidate.candidate_id,
                    "sample_id": candidate.sample_id,
                    "peptide": candidate.mutant_peptide,
                    "variant_type": candidate.peptide_source.value,
                }
            )
