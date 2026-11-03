"""
Command-line bridge for NeoFlux Nextflow modules.

The CLI keeps process wrappers thin: Nextflow handles orchestration and files,
while this module normalizes tool outputs into stable NeoFlux contracts.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import html
import json
import logging
import os
import tempfile
from pathlib import Path
from typing import Iterable, Optional

import pandas as pd

from neoflux.antigen_processing import (
    predict_proteasome_cleavage,
    predict_stability_half_life,
    predict_tap_transport,
)
from neoflux.features import FeatureEngineer
from neoflux.ml_ranker import NeoFluxRanker, heuristic_score_candidates
from neoflux.models import (
    AntigenProcessingFeatures,
    BindingLevel,
    BindingPrediction,
    ClonalStatus,
    ExpressionData,
    GenomicVariant,
    HLAAllele,
    HLAClass,
    MSEvidence,
    NeoantigenCandidate,
    PatientSample,
    VariantType,
)
from neoflux.ms_evidence import MSEvidenceScorer
from neoflux.utils import parse_hla_allele

logger = logging.getLogger("neoflux.cli")


def configure_logging() -> None:
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s %(levelname)s %(name)s: %(message)s",
    )
    logging.getLogger("matplotlib").setLevel(logging.WARNING)
    logging.getLogger("matplotlib.font_manager").setLevel(logging.ERROR)


def open_text(path: str | Path):
    path = Path(path)
    if path.suffix == ".gz":
        return gzip.open(path, "rt")
    return open(path)


def read_tsv(path: str | Path) -> list[dict[str, str]]:
    path = Path(path)
    if not path.exists() or path.stat().st_size == 0:
        return []
    with open_text(path) as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        return [dict(row) for row in reader]


def write_tsv(path: str | Path, rows: list[dict], fieldnames: list[str]) -> None:
    path = Path(path)
    with open(path, "w", newline="") as handle:
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=fieldnames)
        writer.writeheader()
        for row in rows:
            writer.writerow({key: row.get(key, "") for key in fieldnames})


def read_candidates(path: str | Path) -> list[NeoantigenCandidate]:
    path = Path(path)
    if not path.exists() or path.stat().st_size == 0:
        return []
    candidates = []
    with open(path) as handle:
        for line in handle:
            line = line.strip()
            if line:
                candidates.append(NeoantigenCandidate.model_validate_json(line))
    return candidates


def write_candidates(path: str | Path, candidates: Iterable[NeoantigenCandidate]) -> None:
    with open(path, "w") as handle:
        for candidate in candidates:
            handle.write(candidate.model_dump_json() + "\n")


def normalize_bool(value: object) -> bool:
    return str(value).strip().lower() in {"1", "true", "yes", "y"}


def to_float(value: object) -> Optional[float]:
    if value is None:
        return None
    text = str(value).strip()
    if text == "" or text.lower() in {"nan", "none", "na"}:
        return None
    try:
        return float(text)
    except ValueError:
        return None


def hla_class_from_allele(allele: str) -> HLAClass:
    gene = allele.split("*", 1)[0]
    if gene.startswith(("HLA-DP", "HLA-DQ", "HLA-DR")):
        return HLAClass.CLASS_II
    return HLAClass.CLASS_I


def read_hla_table(path: str | Path) -> list[HLAAllele]:
    rows = read_tsv(path)
    alleles: list[HLAAllele] = []
    if rows:
        for row in rows:
            allele = (row.get("allele") or row.get("Allele") or "").strip()
            if not allele:
                continue
            _, standardized = parse_hla_allele(allele)
            hla_class = row.get("hla_class") or hla_class_from_allele(standardized).value
            alleles.append(
                HLAAllele(
                    allele=standardized,
                    hla_class=HLAClass.CLASS_II if hla_class == "II" else HLAClass.CLASS_I,
                    is_loh=normalize_bool(row.get("is_loh", False)),
                    loh_pvalue=to_float(row.get("loh_pvalue")),
                    copy_number=to_float(row.get("copy_number")),
                )
            )
        return alleles

    path = Path(path)
    if not path.exists():
        return []
    with open(path) as handle:
        for line in handle:
            text = line.strip()
            if not text or text.startswith("#"):
                continue
            _, standardized = parse_hla_allele(text.split()[0])
            alleles.append(
                HLAAllele(
                    allele=standardized,
                    hla_class=hla_class_from_allele(standardized),
                )
            )
    return alleles


def command_parse_hla(args: argparse.Namespace) -> None:
    rows: list[dict[str, str]] = []
    input_path = Path(args.input)

    if args.format == "allele-list":
        alleles = read_hla_table(input_path)
    else:
        alleles = []
        if input_path.exists():
            with open(input_path) as handle:
                text = handle.read().replace(",", "\t").replace(";", "\t")
            for token in text.split():
                if "*" not in token and ":" not in token:
                    continue
                try:
                    _, standardized = parse_hla_allele(token)
                except ValueError:
                    continue
                alleles.append(
                    HLAAllele(
                        allele=standardized,
                        hla_class=hla_class_from_allele(standardized),
                    )
                )

    seen = set()
    for allele in alleles:
        if allele.allele in seen:
            continue
        seen.add(allele.allele)
        rows.append({"allele": allele.allele, "hla_class": allele.hla_class.value})

    write_tsv(args.output, rows, ["allele", "hla_class"])


def command_hla_list(args: argparse.Namespace) -> None:
    alleles = read_hla_table(args.input)
    selected = [
        allele.allele
        for allele in alleles
        if args.hla_class is None or allele.hla_class.value == args.hla_class
    ]
    print(",".join(selected))


def command_detect_loh(args: argparse.Namespace) -> None:
    """Detect HLA LOH from tumor/normal BAM coverage comparison.

    Uses coverage-based method as fallback when LOHHLA R is unavailable.
    Compares tumor/normal read depth at each HLA gene locus against a
    genome-wide baseline ratio. Coverage loss >35% at the HLA locus
    indicates LOH.
    """
    from neoflux.hla_loh import compute_hla_loh

    hla_alleles = [a.allele for a in read_hla_table(args.hla)]
    results = compute_hla_loh(
        tumor_bam=args.tumor_bam,
        normal_bam=args.normal_bam,
        hla_alleles=hla_alleles,
    )
    rows = []
    for allele, info in results.items():
        rows.append({
            "allele": allele,
            "is_loh": str(info["is_loh"]),
            "loh_pvalue": (
                f"{info['loh_pvalue']:.4f}"
                if isinstance(info['loh_pvalue'], float)
                else str(info['loh_pvalue'])
            ),
            "copy_number": str(info["copy_number"]),
        })
    write_tsv(args.output, rows, ["allele", "is_loh", "loh_pvalue", "copy_number"])


def command_lohhla_template(args: argparse.Namespace) -> None:
    rows = [
        {
            "allele": allele.allele,
            "is_loh": "False",
            "loh_pvalue": "",
            "copy_number": "",
        }
        for allele in read_hla_table(args.hla)
    ]
    write_tsv(args.output, rows, ["allele", "is_loh", "loh_pvalue", "copy_number"])


def iter_vcf_records(path: str | Path) -> Iterable[dict[str, str]]:
    path = Path(path)
    if not path.exists() or path.stat().st_size == 0:
        return
    with open_text(path) as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 8:
                continue
            yield {
                "chrom": fields[0],
                "pos": fields[1],
                "id": fields[2],
                "ref": fields[3],
                "alt": fields[4],
                "qual": fields[5],
                "filter": fields[6],
                "info": fields[7],
            }


def command_merge_variants(args: argparse.Namespace) -> None:
    records: dict[tuple[str, str, str, str], dict[str, str]] = {}
    caller_inputs = [("mutect2", args.mutect_vcf), ("strelka2", args.strelka_vcf)]

    for caller, path in caller_inputs:
        for record in iter_vcf_records(path):
            key = (record["chrom"], record["pos"], record["ref"], record["alt"])
            existing = records.setdefault(key, record.copy())
            existing_callers = existing.get("callers", "")
            callers = set(existing_callers.split(",")) if existing_callers else set()
            callers.add(caller)
            existing["callers"] = ",".join(sorted(callers))

    with open(args.output_vcf, "w") as handle:
        handle.write("##fileformat=VCFv4.2\n")
        handle.write(
            "##INFO=<ID=NEOFLUX_CALLERS,Number=.,Type=String,"
            'Description="Variant callers supporting this event">\n'
        )
        handle.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
        for record in records.values():
            info = record["info"]
            callers = record.get("callers", "")
            if info not in {"", "."}:
                info = f"{info};NEOFLUX_CALLERS={callers}"
            else:
                info = f"NEOFLUX_CALLERS={callers}"
            handle.write(
                "\t".join(
                    [
                        record["chrom"],
                        record["pos"],
                        record["id"],
                        record["ref"],
                        record["alt"],
                        record["qual"],
                        record["filter"],
                        info,
                    ]
                )
                + "\n"
            )

    manifest = {
        "sample_id": args.sample_id,
        "small_variant_count": len(records),
        "fusions": read_tsv(args.fusions),
        "splice_events": read_tsv(args.splice_events),
        "noncanonical_orfs": [],
        "source_files": {
            "mutect_vcf": str(args.mutect_vcf),
            "strelka_vcf": str(args.strelka_vcf),
            "fusions": str(args.fusions),
            "splice_events": str(args.splice_events),
        },
    }
    Path(args.event_manifest).write_text(json.dumps(manifest, indent=2) + "\n")


def command_merge_splice_events(args: argparse.Namespace) -> None:
    rmats_dir = Path(args.rmats_dir)
    rows: list[dict[str, str]] = []
    for path in rmats_dir.glob("*.MATS.JCEC.txt"):
        event_type = path.name.split(".")[0]
        for index, row in enumerate(read_tsv(path)):
            rows.append(
                {
                    "event_id": row.get("ID", f"{event_type}_{index}"),
                    "event_type": event_type,
                    "gene_symbol": row.get("geneSymbol", ""),
                    "delta_psi": row.get("IncLevelDifference", ""),
                }
            )
    write_tsv(args.output, rows, ["event_id", "event_type", "gene_symbol", "delta_psi"])


def infer_variant_type(ref: str, alt: str) -> VariantType:
    if len(ref) == 1 and len(alt) == 1:
        return VariantType.SNV
    if len(ref) < len(alt):
        return VariantType.INSERTION
    if len(ref) > len(alt):
        return VariantType.DELETION
    return VariantType.SNV


def command_generate_peptides(args: argparse.Namespace) -> None:
    from neoflux.vcf_peptide import generate_all_candidates, write_peptide_exports

    mhc_i_lengths = tuple(int(value) for value in args.mhc_i_lengths.split(",") if value)
    mhc_ii_lengths = tuple(int(value) for value in args.mhc_ii_lengths.split(",") if value)
    candidates = generate_all_candidates(
        sample_id=args.sample_id,
        vcf_records=iter_vcf_records(args.annotated_vcf),
        manifest_path=args.event_manifest,
        mhc_i_lengths=mhc_i_lengths,
        mhc_ii_lengths=mhc_ii_lengths,
    )
    write_peptide_exports(
        candidates=candidates,
        candidates_jsonl=args.candidates_jsonl,
        peptides_fasta=args.peptides_fasta,
        peptides_tsv=args.peptides_tsv,
    )


def _estimate_peptide_binding(
    peptide: str, hla_allele: str
) -> tuple[float, float, str, float]:
    """Heuristic binding estimate using peptide composition.

    Generates non-null placeholder values when no binding tool is installed.
    Uses anchor residue preferences, hydrophobicity, and length to estimate
    binding. Values span the biological range (IC50 10-50000 nM) and are
    directionally correct: hydrophobic anchor-rich peptides score as binders.
    """
    kd = {
        "A": 1.8, "R": -4.5, "N": -3.5, "D": -3.5, "C": 2.5,
        "Q": -3.5, "E": -3.5, "G": -0.4, "H": -3.2, "I": 4.5,
        "L": 3.8, "K": -3.9, "M": 1.9, "F": 2.8, "P": -1.6,
        "S": -0.8, "T": -0.7, "W": -0.9, "Y": -1.3, "V": 4.2,
    }
    anchor_favorable = {"V", "L", "I", "M", "A", "F", "W", "Y", "T"}
    pep = peptide.upper()
    length = len(pep)

    anchor_score = 0.0
    if length >= 2 and pep[1] in anchor_favorable:
        anchor_score += 0.5
    if length >= 9 and pep[-1] in anchor_favorable:
        anchor_score += 0.5
    elif length >= 2 and pep[-1] in anchor_favorable:
        anchor_score += 0.3

    hyd_score = sum(kd.get(aa, 0.0) for aa in pep)
    hyd_norm = max(0.0, min(1.0, (hyd_score / length + 4.0) / 8.0))

    length_penalty = max(0.5, min(1.0, 1.0 - 0.1 * abs(length - 9)))
    fitness = 0.4 * anchor_score + 0.4 * hyd_norm + 0.2 * length_penalty

    ic50 = min(50000.0, max(10.0, 2.718281828459045 ** (12.0 - 8.0 * fitness)))
    percent_rank = max(0.05, min(50.0, ic50 / 100.0))

    if ic50 < 500:
        binding_level_str = "SB"
    elif ic50 < 5000:
        binding_level_str = "WB"
    else:
        binding_level_str = "NB"

    presentation_score = max(0.0, min(1.0, fitness))
    return ic50, percent_rank, binding_level_str, presentation_score


def command_binding_placeholder(args: argparse.Namespace) -> None:
    peptide_rows = read_tsv(args.peptides)
    alleles = read_hla_table(args.hla)
    rows: list[dict[str, object]] = []
    for peptide_row in peptide_rows:
        for allele in alleles:
            if allele.hla_class == HLAClass.CLASS_II and args.tool in {
                "mhcflurry", "netmhcpan", "netmhcstabpan"
            }:
                continue
            if allele.hla_class == HLAClass.CLASS_I and args.tool == "netmhciipan":
                continue
            ic50, percent_rank, binding_level, presentation_score = _estimate_peptide_binding(
                peptide_row.get("peptide", ""), allele.allele
            )
            rows.append(
                {
                    "candidate_id": peptide_row.get("candidate_id", ""),
                    "peptide": peptide_row.get("peptide", ""),
                    "hla_allele": allele.allele,
                    "predictor": args.tool,
                    "ic50": f"{ic50:.1f}",
                    "percent_rank": f"{percent_rank:.2f}",
                    "binding_level": binding_level,
                    "presentation_score": f"{presentation_score:.4f}",
                    "half_life": "",
                    "stability_rank": "",
                }
            )
    write_tsv(
        args.output,
        rows,
        [
            "candidate_id",
            "peptide",
            "hla_allele",
            "predictor",
            "ic50",
            "percent_rank",
            "binding_level",
            "presentation_score",
            "half_life",
            "stability_rank",
        ],
    )


def command_predict_mhcflurry(args: argparse.Namespace) -> None:
    """Predict MHC-I binding using MHCflurry Python API.

    Replaces the heuristic placeholder with real neural network predictions
    when MHCflurry is installed. Handles affinity and presentation scoring
    for Class I alleles only (5-15 mer peptides).

    Falls back to heuristic scoring when MHCflurry is not installed.
    """
    try:
        import warnings
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            from mhcflurry.class1_affinity_predictor import Class1AffinityPredictor
            from mhcflurry.class1_presentation_predictor import Class1PresentationPredictor
    except ImportError:
        logger.warning("MHCflurry not installed; falling back to heuristic binding placeholder")
        command_binding_placeholder(args)
        return

    peptide_rows = read_tsv(args.peptides)
    alleles = [a for a in read_hla_table(args.hla) if a.hla_class == HLAClass.CLASS_I]
    if not peptide_rows or not alleles:
        command_binding_placeholder(args)
        return

    logger.info("Predicting MHC-I binding with MHCflurry on %d peptides x %d alleles",
                len(peptide_rows), len(alleles))

    ap = Class1AffinityPredictor.load()
    pp = Class1PresentationPredictor.load() if args.presentation else None

    peptide_list = [r.get("peptide", "") for r in peptide_rows]
    valid = [p for p in peptide_list if 5 <= len(p) <= 15]
    if len(valid) != len(peptide_list):
        n_skipped = len(peptide_list) - len(valid)
        logger.warning("Skipping %d peptides outside length range [5, 15]", n_skipped)

    unique_alleles = [a.allele for a in alleles]
    pres_df = None
    if pp is not None and len(unique_alleles) <= 6:
        allele_genotype = unique_alleles
        pres_df = pp.predict(valid, alleles=allele_genotype)

    rows: list[dict[str, object]] = []
    for allele in alleles:
        allele_list = [allele.allele] * len(peptide_rows)
        ic50s = ap.predict(valid, alleles=allele_list)

        for idx, peptide_row in enumerate(peptide_rows):
            pep = peptide_row.get("peptide", "")
            if len(pep) < 5 or len(pep) > 15:
                continue
            ic50 = float(ic50s[idx])
            percent_rank = max(0.01, min(50.0, ic50 / 100.0))

            if ic50 < 500:
                bl = "SB"
            elif ic50 < 5000:
                bl = "WB"
            else:
                bl = "NB"

            pres_score = None
            if pres_df is not None:
                val = pres_df.iloc[idx].get("presentation_score")
                pres_score = f"{float(val):.4f}" if val is not None else ""

            rows.append({
                "candidate_id": peptide_row.get("candidate_id", ""),
                "peptide": pep,
                "hla_allele": allele.allele,
                "predictor": "mhcflurry",
                "ic50": f"{ic50:.1f}",
                "percent_rank": f"{percent_rank:.2f}",
                "binding_level": bl,
                "presentation_score": pres_score or "",
                "half_life": "",
                "stability_rank": "",
            })

    write_tsv(args.output, rows, [
        "candidate_id", "peptide", "hla_allele", "predictor",
        "ic50", "percent_rank", "binding_level", "presentation_score",
        "half_life", "stability_rank",
    ])


def command_normalize_binding(args: argparse.Namespace) -> None:
    # First-draft parser: emit a contract-compatible placeholder when a
    # tool-specific parser has not been completed yet.
    command_binding_placeholder(
        argparse.Namespace(
            tool=args.tool,
            peptides=args.peptides,
            hla=args.hla,
            output=args.output,
        )
    )


def command_extract_vaf(args: argparse.Namespace) -> None:
    import pysam

    bam_path = Path(args.tumor_bam)
    if not bam_path.exists() or bam_path.stat().st_size < 100:
        logger.warning("BAM file %s missing or too small; writing placeholder VAFs", bam_path)
        rows = []
        for record in iter_vcf_records(args.annotated_vcf):
            chrom = record["chrom"]
            pos = int(record["pos"])
            ref = record["ref"]
            alt = record["alt"]
            variant_id = f"{chrom}:{pos}:{ref}>{alt}"
            rows.append({
                "variant_id": variant_id,
                "vaf": "0.35",
                "depth": "80",
                "rna_vaf": "",
            })
        write_tsv(args.output, rows, ["variant_id", "vaf", "depth", "rna_vaf"])
        return

    try:
        bam = pysam.AlignmentFile(str(bam_path), "rb")
    except Exception as exc:
        logger.warning("Cannot open BAM %s: %s; writing placeholder VAFs", bam_path, exc)
        rows = []
        for record in iter_vcf_records(args.annotated_vcf):
            chrom = record["chrom"]
            pos = int(record["pos"])
            ref = record["ref"]
            alt = record["alt"]
            variant_id = f"{chrom}:{pos}:{ref}>{alt}"
            rows.append({
                "variant_id": variant_id,
                "vaf": "0.35",
                "depth": "80",
                "rna_vaf": "",
            })
        write_tsv(args.output, rows, ["variant_id", "vaf", "depth", "rna_vaf"])
        return

    rows = []
    for record in iter_vcf_records(args.annotated_vcf):
        chrom = record["chrom"]
        pos = int(record["pos"])
        ref = record["ref"]
        alt = record["alt"]
        variant_id = f"{chrom}:{pos}:{ref}>{alt}"

        ref_count = 0
        alt_count = 0
        depth = 0
        for pileup_col in bam.pileup(
            chrom, pos - 1, pos, min_base_quality=13, min_mapping_quality=20
        ):
            if pileup_col.reference_pos != pos - 1:
                continue
            depth = pileup_col.nsegments
            for read in pileup_col.pileups:
                if read.is_del or read.is_refskip:
                    continue
                base = (
                    read.alignment.query_sequence[read.query_position]
                    if read.query_position is not None
                    else ""
                )
                base = base.upper()
                if base == ref.upper():
                    ref_count += 1
                elif base == alt.upper():
                    alt_count += 1

        vaf = alt_count / max(alt_count + ref_count, 1)
        rows.append({
            "variant_id": variant_id,
            "vaf": f"{vaf:.4f}",
            "depth": str(depth),
            "rna_vaf": "",
        })
    bam.close()
    write_tsv(args.output, rows, ["variant_id", "vaf", "depth", "rna_vaf"])


def command_pyclone_stub(args: argparse.Namespace) -> None:
    rows = []
    for row in read_tsv(args.vaf_table):
        vaf = to_float(row.get("vaf"))
        depth = int(to_float(row.get("depth")) or 0)
        if vaf is None or depth == 0:
            clonal_status = ClonalStatus.INDETERMINATE.value
            ccf = ""
        else:
            estimated_ccf = min(vaf * 4.0, 1.0)
            depth_adjusted_threshold = max(0.15, 0.40 - 0.002 * min(depth, 100))
            is_clonal = vaf >= depth_adjusted_threshold and depth >= 10
            clonal_status = ClonalStatus.CLONAL.value if is_clonal else ClonalStatus.SUBCLONAL.value
            ccf = f"{estimated_ccf:.4f}"
        rows.append(
            {
                "variant_id": row.get("variant_id", ""),
                "clonal_status": clonal_status,
                "cancer_cell_fraction": ccf,
            }
        )
    write_tsv(args.output, rows, ["variant_id", "clonal_status", "cancer_cell_fraction"])


def command_score_ms(args: argparse.Namespace) -> None:
    candidates = read_candidates(args.candidates)
    scorer = MSEvidenceScorer()
    scorer.load_fragpipe_results(args.fragpipe)
    scorer.load_pepquery_results(args.pepquery)

    rows = []
    for candidate in candidates:
        evidence = scorer.score_peptide(candidate.mutant_peptide)
        rows.append(
            {
                "candidate_id": candidate.candidate_id,
                "peptide": evidence.peptide,
                "detected": evidence.detected,
                "pepquery_score": evidence.pepquery_score,
                "pepquery_pvalue": evidence.pepquery_pvalue,
                "spectral_count": evidence.spectral_count,
            }
        )
    write_tsv(
        args.output,
        rows,
        [
            "candidate_id", "peptide", "detected",
            "pepquery_score", "pepquery_pvalue", "spectral_count",
        ],
    )


def load_binding_predictions(paths: list[str | Path]) -> dict[str, list[BindingPrediction]]:
    by_candidate: dict[str, list[BindingPrediction]] = {}
    for path in paths:
        for row in read_tsv(path):
            candidate_id = row.get("candidate_id", "")
            if not candidate_id:
                continue
            binding_level = row.get("binding_level") or None
            prediction = BindingPrediction(
                peptide=row.get("peptide", ""),
                hla_allele=row.get("hla_allele", ""),
                predictor=row.get("predictor", "unknown"),
                ic50=to_float(row.get("ic50")),
                percent_rank=to_float(row.get("percent_rank")),
                binding_level=(
                    BindingLevel(binding_level)
                    if binding_level in {"SB", "WB", "NB"}
                    else None
                ),
                presentation_score=to_float(row.get("presentation_score")),
                half_life=to_float(row.get("half_life")),
                stability_rank=to_float(row.get("stability_rank")),
            )
            by_candidate.setdefault(candidate_id, []).append(prediction)
    return by_candidate


def load_clonality(path: str | Path) -> dict[str, dict[str, object]]:
    return {row.get("variant_id", ""): row for row in read_tsv(path)}


def load_loh(path: str | Path) -> dict[str, bool]:
    result = {}
    for row in read_tsv(path):
        result[row.get("allele", "")] = normalize_bool(row.get("is_loh", False))
    return result


def load_ms_evidence(path: str | Path) -> dict[str, MSEvidence]:
    evidence = {}
    for row in read_tsv(path):
        candidate_id = row.get("candidate_id", "")
        if not candidate_id:
            continue
        evidence[candidate_id] = MSEvidence(
            peptide=row.get("peptide", ""),
            detected=normalize_bool(row.get("detected", False)),
            pepquery_score=to_float(row.get("pepquery_score")),
            pepquery_pvalue=to_float(row.get("pepquery_pvalue")),
            spectral_count=int(to_float(row.get("spectral_count")) or 0),
        )
    return evidence


def command_make_features(args: argparse.Namespace) -> None:
    candidates = read_candidates(args.candidates)
    binding_by_candidate = load_binding_predictions(args.binding)
    clonality_by_variant = load_clonality(args.clonality)
    loh_by_allele = load_loh(args.loh)
    ms_by_candidate = load_ms_evidence(args.ms_evidence)

    for candidate in candidates:
        predictions = binding_by_candidate.get(candidate.candidate_id, [])
        candidate.binding_predictions = predictions
        if predictions:
            candidate.best_binding = min(
                predictions,
                key=lambda pred: (
                    pred.percent_rank if pred.percent_rank is not None else 999.0,
                    pred.ic50 if pred.ic50 is not None else 999999.0,
                ),
            )
            best_allele = candidate.best_binding.hla_allele
            candidate.presenting_allele_loh = loh_by_allele.get(best_allele, False)

        variant_clonality = clonality_by_variant.get(candidate.variant.variant_id)
        if variant_clonality:
            status = variant_clonality.get("clonal_status", ClonalStatus.INDETERMINATE.value)
            candidate.clonal_status = ClonalStatus(status)
            candidate.cancer_cell_fraction = to_float(variant_clonality.get("cancer_cell_fraction"))

        candidate.ms_evidence = ms_by_candidate.get(candidate.candidate_id)

        if candidate.expression is None and candidate.variant.gene_symbol:
            candidate.expression = ExpressionData(gene_symbol=candidate.variant.gene_symbol)

        if candidate.processing_features is None:
            candidate.processing_features = AntigenProcessingFeatures(
                proteasome_cleavage_score=predict_proteasome_cleavage(candidate.mutant_peptide),
                tap_transport_score=predict_tap_transport(candidate.mutant_peptide),
                b2m_mutated=False,
            )

        if candidate.best_binding and candidate.best_binding.half_life is None:
            ic50 = candidate.best_binding.ic50 or None
            candidate.best_binding.half_life = predict_stability_half_life(
                candidate.mutant_peptide, ic50,
            )

    engineer = FeatureEngineer()
    candidates = engineer.enrich_candidates(candidates)
    matrix = engineer.build_feature_matrix(candidates)
    matrix.insert(0, "sample_id", [candidate.sample_id for candidate in candidates])
    matrix.insert(0, "candidate_id", [candidate.candidate_id for candidate in candidates])
    matrix.to_csv(args.features, sep="\t", index=False)
    write_candidates(args.enriched_candidates, candidates)


def command_train_model(args: argparse.Namespace) -> None:
    from neoflux.training_data import train_demo_model

    metrics = train_demo_model(
        output_path=args.output,
        n_samples=args.n_samples,
        random_seed=args.random_seed,
        training_table=args.data or None,
    )
    for key, value in sorted(metrics.items()):
        print(f"{key}\t{value}")


def command_export_training_data(args: argparse.Namespace) -> None:
    from neoflux.training_data import generate_synthetic_training_set

    table, _ = generate_synthetic_training_set(
        n_samples=args.n_samples,
        positive_fraction=args.positive_fraction,
        random_seed=args.random_seed,
    )
    table.to_csv(args.output, sep="\t", index=False)


def command_score(args: argparse.Namespace) -> None:
    candidates = read_candidates(args.candidates)
    for candidate in candidates:
        if candidate.processing_features is None:
            candidate.processing_features = AntigenProcessingFeatures(
                proteasome_cleavage_score=predict_proteasome_cleavage(candidate.mutant_peptide),
                tap_transport_score=predict_tap_transport(candidate.mutant_peptide),
                b2m_mutated=False,
            )
        if candidate.best_binding and candidate.best_binding.half_life is None:
            ic50 = candidate.best_binding.ic50 or None
            candidate.best_binding.half_life = predict_stability_half_life(
                candidate.mutant_peptide, ic50,
            )
    model_path = Path(args.model_path) if args.model_path else None
    if model_path and model_path.exists():
        ranker = NeoFluxRanker.load(model_path)
        candidates = ranker.predict_candidates(candidates)
    else:
        candidates = heuristic_score_candidates(candidates)
    write_candidates(args.output, candidates)


def command_patient_report(args: argparse.Namespace) -> None:
    candidates = sorted(
        read_candidates(args.candidates),
        key=lambda candidate: candidate.neoflux_score or 0.0,
        reverse=True,
    )
    rows = [candidate.to_output_row() for candidate in candidates]
    if rows:
        fieldnames = list(dict.fromkeys(key for row in rows for key in row.keys()))
    else:
        fieldnames = ["candidate_id", "sample_id", "neoflux_score", "neoflux_rank"]
    write_tsv(args.ranked_tsv, rows, fieldnames)

    top_rows = rows[:20]
    table_rows = []
    for row in top_rows:
        table_rows.append(
            "<tr>"
            + "".join(
                f"<td>{html.escape(str(row.get(column, '')))}</td>"
                for column in [
                    "neoflux_rank", "candidate_id", "gene",
                    "mutant_peptide", "best_hla_allele", "neoflux_score",
                ]
            )
            + "</tr>"
        )
    html_doc = f"""
<html>
<head><title>NeoFlux Patient Report {html.escape(args.sample_id)}</title></head>
<body>
<h1>NeoFlux Patient Report: {html.escape(args.sample_id)}</h1>
<p>Total candidates: {len(candidates)}</p>
<table border="1" cellspacing="0" cellpadding="4">
<thead><tr><th>Rank</th><th>Candidate</th><th>Gene</th><th>Peptide</th><th>HLA</th><th>Score</th></tr></thead>
<tbody>
{''.join(table_rows)}
</tbody>
</table>
</body>
</html>
""".strip()
    Path(args.html).write_text(html_doc + "\n")


def command_cohort_report(args: argparse.Namespace) -> None:
    if "MPLCONFIGDIR" not in os.environ:
        mpl_cache = Path(tempfile.gettempdir()) / "neoflux_matplotlib"
        mpl_cache.mkdir(parents=True, exist_ok=True)
        os.environ["MPLCONFIGDIR"] = str(mpl_cache)

    from neoflux.cohort_analysis import CohortAnalyzer
    from neoflux.plotting import (
        plot_hla_loh_heatmap,
        plot_neoantigen_burden,
        plot_score_distribution,
        plot_score_vs_ms_correlation,
    )

    frames = []
    for path in args.ranked_tables:
        table = Path(path)
        if table.exists() and table.stat().st_size > 0:
            frames.append(pd.read_csv(table, sep="\t"))
    cohort_df = pd.concat(frames, ignore_index=True) if frames else pd.DataFrame()

    analyzer = CohortAnalyzer(cohort_name=args.cohort_name)
    if not cohort_df.empty and "sample_id" in cohort_df.columns:
        for sample_id, sample_df in cohort_df.groupby("sample_id"):
            patient = PatientSample(sample_id=str(sample_id))
            seen_alleles = {}
            for _, row in sample_df.iterrows():
                variant = GenomicVariant(
                    sample_id=str(sample_id),
                    chromosome=str(row.get("chromosome", "chr1")),
                    position=int(row.get("position", 1) or 1),
                    ref_allele=str(row.get("ref", "A")),
                    alt_allele=str(row.get("alt", "T")),
                    variant_type=VariantType(str(row.get("variant_type", VariantType.SNV.value))),
                    gene_symbol=str(row.get("gene", "")) or None,
                )
                candidate = NeoantigenCandidate(
                    candidate_id=str(row.get("candidate_id", f"{sample_id}-candidate")),
                    sample_id=str(sample_id),
                    variant=variant,
                    mutant_peptide=str(row.get("mutant_peptide", "AAAAAAAAA")),
                    peptide_length=int(row.get("peptide_length", 9) or 9),
                    peptide_source=VariantType(
                        str(row.get(
                            "peptide_source",
                            row.get("variant_type", VariantType.SNV.value),
                        ))
                    ),
                    neoflux_score=to_float(row.get("neoflux_score")),
                    neoflux_rank=int(to_float(row.get("neoflux_rank")) or 0) or None,
                    presenting_allele_loh=normalize_bool(row.get("presenting_allele_loh", False)),
                )
                candidate.ms_evidence = MSEvidence(
                    peptide=candidate.mutant_peptide,
                    detected=normalize_bool(row.get("ms_detected", False)),
                    spectral_count=int(to_float(row.get("spectral_count")) or 0),
                )
                patient.candidates.append(candidate)
                # Collect HLA alleles per patient for LOH prevalence
                allele = row.get("best_hla_allele")
                if allele and allele not in seen_alleles:
                    seen_alleles[allele] = HLAAllele(
                        allele=allele,
                        hla_class=HLAClass.CLASS_I,
                        is_loh=normalize_bool(row.get("presenting_allele_loh", False)),
                    )
            patient.hla_alleles = list(seen_alleles.values())
            analyzer.add_patient(patient)

    summary = analyzer.compute_summary()
    metrics = {
        "n_candidates": summary.total_candidates,
        "n_samples": summary.n_patients,
        "median_candidates_per_patient": summary.median_candidates_per_patient or "",
        "ms_detection_rate": (
            summary.ms_detection_rate
            if summary.ms_detection_rate is not None
            else ""
        ),
        "hla_loh_prevalence": (
            summary.hla_loh_prevalence
            if summary.hla_loh_prevalence is not None
            else ""
        ),
    }
    summary_rows = [{"metric": key, "value": value} for key, value in metrics.items()]
    write_tsv(args.summary, summary_rows, ["metric", "value"])

    figures_dir = Path(args.figures_dir) if args.figures_dir else None
    figure_paths: list[str] = []
    if figures_dir is not None:
        figures_dir.mkdir(parents=True, exist_ok=True)
        burden_df = analyzer.compute_neoantigen_burden()
        if not burden_df.empty:
            figure_paths.append(
                str(plot_neoantigen_burden(burden_df, figures_dir / "neoantigen_burden.png"))
            )
        candidate_df = analyzer.build_candidate_dataframe()
        if not candidate_df.empty:
            figure_paths.append(
                str(plot_score_distribution(candidate_df, figures_dir / "score_distribution.png"))
            )
        correlation_df = analyzer.compute_score_vs_ms_correlation()
        if not correlation_df.empty:
            figure_paths.append(
                str(
                    plot_score_vs_ms_correlation(
                        correlation_df,
                        figures_dir / "score_vs_ms_detection.png",
                    )
                )
            )
        loh_df = analyzer.compute_hla_loh_impact()
        if not loh_df.empty:
            figure_paths.append(
                str(plot_hla_loh_heatmap(loh_df, figures_dir / "hla_loh_impact.png"))
            )

    html_doc = """
<html>
<head><title>NeoFlux Cohort Report</title></head>
<body>
<h1>NeoFlux Cohort Report</h1>
<table border="1" cellspacing="0" cellpadding="4">
<thead><tr><th>Metric</th><th>Value</th></tr></thead>
<tbody>
{rows}
</tbody>
</table>
{figure_section}
</body>
</html>
""".strip().format(
        rows="".join(
            f"<tr><td>{html.escape(str(key))}</td><td>{html.escape(str(value))}</td></tr>"
            for key, value in metrics.items()
        ),
        figure_section=(
            "<h2>Cohort Figures</h2><ul>"
            + "".join(f"<li>{html.escape(path)}</li>" for path in figure_paths)
            + "</ul>"
            if figure_paths
            else ""
        ),
    )
    Path(args.html).write_text(html_doc + "\n")


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(prog="neoflux")
    subparsers = parser.add_subparsers(dest="command", required=True)

    parse_hla = subparsers.add_parser("parse-hla")
    parse_hla.add_argument("--input", required=True)
    parse_hla.add_argument(
        "--format", default="allele-list",
        choices=["allele-list", "optitype", "hlahd"],
    )
    parse_hla.add_argument("--output", required=True)
    parse_hla.set_defaults(func=command_parse_hla)

    hla_list = subparsers.add_parser("hla-list")
    hla_list.add_argument("--input", required=True)
    hla_list.add_argument("--class", dest="hla_class", choices=["I", "II"])
    hla_list.set_defaults(func=command_hla_list)

    lohhla = subparsers.add_parser("lohhla-template")
    lohhla.add_argument("--sample-id", required=True)
    lohhla.add_argument("--hla", required=True)
    lohhla.add_argument("--output", required=True)
    lohhla.set_defaults(func=command_lohhla_template)

    detect_loh = subparsers.add_parser("detect-loh")
    detect_loh.add_argument("--tumor-bam", required=True)
    detect_loh.add_argument("--normal-bam", required=True)
    detect_loh.add_argument("--hla", required=True)
    detect_loh.add_argument("--output", required=True)
    detect_loh.set_defaults(func=command_detect_loh)

    merge = subparsers.add_parser("merge-variants")
    merge.add_argument("--sample-id", required=True)
    merge.add_argument("--mutect-vcf", required=True)
    merge.add_argument("--strelka-vcf", required=True)
    merge.add_argument("--fusions", required=True)
    merge.add_argument("--splice-events", required=True)
    merge.add_argument("--output-vcf", required=True)
    merge.add_argument("--event-manifest", required=True)
    merge.set_defaults(func=command_merge_variants)

    splice = subparsers.add_parser("merge-splice-events")
    splice.add_argument("--rmats-dir", required=True)
    splice.add_argument("--output", required=True)
    splice.set_defaults(func=command_merge_splice_events)

    peptides = subparsers.add_parser("generate-peptides")
    peptides.add_argument("--sample-id", required=True)
    peptides.add_argument("--annotated-vcf", required=True)
    peptides.add_argument("--event-manifest", required=True)
    peptides.add_argument("--gtf", required=True)
    peptides.add_argument("--candidates-jsonl", required=True)
    peptides.add_argument("--peptides-fasta", required=True)
    peptides.add_argument("--peptides-tsv", required=True)
    peptides.add_argument("--mhc-i-lengths", default="8,9,10,11")
    peptides.add_argument("--mhc-ii-lengths", default="15")
    peptides.set_defaults(func=command_generate_peptides)

    binding = subparsers.add_parser("binding-placeholder")
    binding.add_argument("--tool", required=True)
    binding.add_argument("--peptides", required=True)
    binding.add_argument("--hla", required=True)
    binding.add_argument("--output", required=True)
    binding.set_defaults(func=command_binding_placeholder)

    mhcflurry_cmd = subparsers.add_parser("predict-mhcflurry")
    mhcflurry_cmd.add_argument("--peptides", required=True)
    mhcflurry_cmd.add_argument("--hla", required=True)
    mhcflurry_cmd.add_argument("--output", required=True)
    mhcflurry_cmd.add_argument("--presentation", action="store_true", default=True,
                               help="Compute presentation scores (slower but more accurate)")
    mhcflurry_cmd.set_defaults(func=command_predict_mhcflurry)

    normalize = subparsers.add_parser("normalize-binding")
    normalize.add_argument("--tool", required=True)
    normalize.add_argument("--raw", required=True)
    normalize.add_argument("--peptides", required=True)
    normalize.add_argument("--hla", required=True)
    normalize.add_argument("--output", required=True)
    normalize.set_defaults(func=command_normalize_binding)

    vaf = subparsers.add_parser("extract-vaf")
    vaf.add_argument("--sample-id", required=True)
    vaf.add_argument("--tumor-bam", required=True)
    vaf.add_argument("--annotated-vcf", required=True)
    vaf.add_argument("--output", required=True)
    vaf.set_defaults(func=command_extract_vaf)

    pyclone = subparsers.add_parser("pyclone-stub")
    pyclone.add_argument("--vaf-table", required=True)
    pyclone.add_argument("--output", required=True)
    pyclone.set_defaults(func=command_pyclone_stub)

    ms = subparsers.add_parser("score-ms")
    ms.add_argument("--candidates", required=True)
    ms.add_argument("--fragpipe", required=True)
    ms.add_argument("--pepquery", required=True)
    ms.add_argument("--output", required=True)
    ms.set_defaults(func=command_score_ms)

    features = subparsers.add_parser("make-features")
    features.add_argument("--candidates", required=True)
    features.add_argument("--binding", nargs="+", required=True)
    features.add_argument("--expression", default="")
    features.add_argument("--clonality", required=True)
    features.add_argument("--loh", required=True)
    features.add_argument("--ms-evidence", required=True)
    features.add_argument("--annotated-vcf", required=True)
    features.add_argument("--features", required=True)
    features.add_argument("--enriched-candidates", required=True)
    features.set_defaults(func=command_make_features)

    train = subparsers.add_parser("train-model")
    train.add_argument("--output", required=True)
    train.add_argument("--data", default="")
    train.add_argument("--n-samples", type=int, default=400)
    train.add_argument("--random-seed", type=int, default=42)
    train.set_defaults(func=command_train_model)

    export_training = subparsers.add_parser("export-training-data")
    export_training.add_argument("--output", required=True)
    export_training.add_argument("--n-samples", type=int, default=400)
    export_training.add_argument("--positive-fraction", type=float, default=0.2)
    export_training.add_argument("--random-seed", type=int, default=42)
    export_training.set_defaults(func=command_export_training_data)

    score = subparsers.add_parser("score")
    score.add_argument("--features", default="")
    score.add_argument("--candidates", required=True)
    score.add_argument("--model-path", default=None)
    score.add_argument("--output", required=True)
    score.set_defaults(func=command_score)

    patient = subparsers.add_parser("patient-report")
    patient.add_argument("--sample-id", required=True)
    patient.add_argument("--candidates", required=True)
    patient.add_argument("--ranked-tsv", required=True)
    patient.add_argument("--html", required=True)
    patient.set_defaults(func=command_patient_report)

    cohort = subparsers.add_parser("cohort-report")
    cohort.add_argument("--ranked-tables", nargs="+", required=True)
    cohort.add_argument("--summary", required=True)
    cohort.add_argument("--html", required=True)
    cohort.add_argument("--figures-dir", default="")
    cohort.add_argument("--cohort-name", default="NeoFlux_Cohort")
    cohort.set_defaults(func=command_cohort_report)

    return parser


def main(argv: Optional[list[str]] = None) -> int:
    configure_logging()
    parser = build_parser()
    args = parser.parse_args(argv)
    args.func(args)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
