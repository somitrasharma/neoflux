#!/usr/bin/env python3
"""Demonstration runner for NeoFlux SITC 2026 abstract.

Generates realistic synthetic data for 3-5 "patients" using known
cancer mutation patterns, runs the full Python pipeline, and produces
output tables + figures ready for the abstract.

Usage:
    python scripts/run_demo.py --output-dir demo_results

This script does NOT need Nextflow or Docker. It calls the Python CLI
directly. Real data can replace the synthetic VCFs later by running:
    python scripts/run_demo.py --vcf-dir my_vcfs --output-dir results
"""

import argparse
import csv
import json
import subprocess
import sys
import tempfile
from pathlib import Path


ROOT = Path(__file__).resolve().parent.parent
CLI = [sys.executable, "-m", "neoflux.cli"]

# Realistic "patients" with known driver mutations and HLA types
PATIENTS = [
    {
        "id": "MEL001",
        "hla_alleles": ["HLA-A*02:01", "HLA-A*01:01", "HLA-B*07:02", "HLA-B*08:01", "HLA-C*07:01"],
        "mutations": [
            # gene, chrom, pos, ref, alt, HGVSp (V600E is a real melanoma driver)
            ("BRAF", "7", "140453136", "A", "T", "p.V600E"),
            ("NRAS", "1", "115256529", "T", "C", "p.Q61R"),
            ("TP53", "17", "7674253", "C", "T", "p.R248W"),
            ("CDKN2A", "9", "21971188", "G", "A", "p.R80*"),
            ("PTEN", "10", "89624250", "A", "G", "p.D326N"),
            ("EGFR", "7", "55249063", "C", "T", "p.G719S"),
        ],
    },
    {
        "id": "MEL002",
        "hla_alleles": ["HLA-A*03:01", "HLA-A*24:02", "HLA-B*35:01", "HLA-B*44:02", "HLA-C*04:01"],
        "mutations": [
            ("BRAF", "7", "140453136", "A", "T", "p.V600E"),
            ("NRAS", "1", "115258747", "G", "C", "p.G12V"),
            ("TP53", "17", "7673802", "C", "T", "p.R175H"),
            ("APC", "5", "112175221", "T", "A", "p.E1317Q"),
        ],
    },
    {
        "id": "MEL003",
        "hla_alleles": ["HLA-A*02:01", "HLA-A*11:01", "HLA-B*15:01", "HLA-B*40:01", "HLA-C*03:04"],
        "mutations": [
            ("NRAS", "1", "115256530", "A", "T", "p.Q61L"),
            ("TP53", "17", "7674220", "G", "A", "p.R273C"),
            ("PTEN", "10", "89692915", "C", "T", "p.R130Q"),
            ("KRAS", "12", "25398281", "G", "A", "p.G12D"),
            ("BRAF", "7", "140453136", "A", "T", "p.V600E"),
            ("MUC16", "19", "8959491", "G", "A", "p.E3405K"),
        ],
    },
    {
        "id": "MEL004",
        "hla_alleles": ["HLA-A*01:01", "HLA-A*32:01", "HLA-B*08:01", "HLA-B*51:01", "HLA-C*06:02"],
        "mutations": [
            ("BRAF", "7", "140453136", "A", "T", "p.V600E"),
            ("TP53", "17", "7675088", "C", "T", "p.Y220C"),
            ("KRAS", "12", "25398284", "C", "T", "p.G12V"),
            ("ARID1A", "1", "27022898", "A", "G", "p.Q1868R"),
            ("SMAD4", "18", "48575510", "T", "C", "p.L440P"),
        ],
    },
    {
        "id": "MEL005",
        "hla_alleles": ["HLA-A*02:01", "HLA-A*68:01", "HLA-B*15:01", "HLA-B*53:01", "HLA-C*04:01"],
        "mutations": [
            ("NRAS", "1", "115258744", "G", "T", "p.G12C"),
            ("TP53", "17", "7674245", "T", "C", "p.V272M"),
            ("CTNNB1", "3", "41266101", "C", "T", "p.S37F"),
            ("EGFR", "7", "55249071", "G", "A", "p.L858R"),
            ("IDH1", "2", "209113113", "G", "A", "p.R132C"),
            ("KEAP1", "19", "10597482", "C", "T", "p.G333C"),
            ("ARID2", "12", "46124961", "T", "C", "p.I1235V"),
        ],
    },
]


def build_parser():
    p = argparse.ArgumentParser(description="Run NeoFlux demo for SITC abstract")
    p.add_argument("--output-dir", default="demo_results", help="Output directory")
    p.add_argument("--n-patients", type=int, default=5, help="Number of patients (max 5)")
    p.add_argument("--skip-mhcflurry", action="store_true",
                   help="Skip MHCflurry (use heuristic placeholder)")
    return p


def make_vcf(patient: dict, output_path: Path) -> None:
    lines = ["##fileformat=VCFv4.2", "##source=NeoFlux_demo"]
    header = ["#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO"]
    lines.append("\t".join(header))

    for gene, chrom, pos, ref, alt, hgvsp in patient["mutations"]:
        info = f"GENE={gene};HGVS_P={hgvsp}"
        lines.append(f"{chrom}\t{pos}\t.\t{ref}\t{alt}\t.\tPASS\t{info}")

    output_path.write_text("\n".join(lines) + "\n")


def make_hla(patient: dict, output_path: Path) -> None:
    with open(output_path, "w", newline="") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(["allele", "hla_class"])
        for allele in patient["hla_alleles"]:
            w.writerow([allele, "I"])


def make_event_manifest(patient: dict, output_path: Path) -> None:
    """Event manifest with non-canonical sources for relevant patients.
    
    Adds realistic fusion and splice events to specific patients so the
    pipeline demonstrates non-canonical candidate generation."""
    pid = patient["id"]
    manifest = {"fusions": [], "splice_events": [], "noncanonical_orfs": []}
    
    if pid == "MEL003":
        # EML4-ALK fusion (canonical lung cancer fusion, E13;A20)
        # Canonical EML4-ALK variant 3 (E13;A20): EML4(1-495) fused to the ALK
        # kinase domain (ALK residue ~1058 to C-terminus). Breakpoints given in
        # protein coordinates against canonical EML4/ALK (UniProt Q9HC35/Q9UM73).
        manifest["fusions"].append({
            "gene_5p": "EML4", "gene_3p": "ALK",
            "break_5p_aa": 495, "break_3p_aa": 1058,
            "chrom": "chr2", "pos": 42522617,
            "breakpoint": "chr2:42522617",
        })
        # Also a non-canonical ORF from a LINE-1 retroelement
        manifest["noncanonical_orfs"].append({
            "orf_id": "L1HS_ORF1",
            "orf_type": "retroelement",
            "protein": "MSTLNRQKSYPGPAGPAGPAGPAGPAGPAGPAGPAGPAGPAGPAGPAGPAGPAGPAGPAGPAGPG",
        })
    elif pid == "MEL005":
        # MET exon 14 skipping (splice variant, ~3% of NSCLC)
        # MET exon 14 encodes the juxtamembrane domain (~residues 963-1009 of
        # canonical MET, P08581); skipping deletes it in-frame, creating a novel
        # junction. skip_start/skip_len place the skipped exon in protein coords.
        manifest["splice_events"].append({
            "gene_symbol": "MET", "event_type": "SE",
            "skip_start": 962, "skip_len": 47,
            "chrom": "chr7", "pos": 116411751,
            "delta_psi": 0.92,
        })
    
    output_path.write_text(json.dumps(manifest) + "\n")


def _cli(*args, cwd=None) -> subprocess.CompletedProcess:
    cmd = CLI + [str(a) for a in args]
    result = subprocess.run(cmd, capture_output=True, text=True, cwd=cwd or ROOT)
    if result.returncode != 0:
        print(f"CLI failed: {' '.join(cmd)}")
        print(f"STDERR: {result.stderr[:500]}")
        result.check_returncode()
    return result


def main():
    args = build_parser().parse_args()
    out = Path(args.output_dir)
    out.mkdir(parents=True, exist_ok=True)

    tmp = Path(tempfile.mkdtemp(prefix="neoflux_demo_"))
    gtf = ROOT / "assets" / "ref.gtf"
    n_patients = min(args.n_patients, len(PATIENTS))

    ranked_tables = []
    summary_rows = []

    for patient in PATIENTS[:n_patients]:
        pid = patient["id"]
        print(f"\n{'='*60}")
        print(f"Processing {pid}...")
        print(f"{'='*60}")

        pdir = out / pid
        pdir.mkdir(exist_ok=True)

        vcf = tmp / f"{pid}.vcf"
        hla = tmp / f"{pid}_hla.tsv"
        events = tmp / f"{pid}_events.json"
        make_vcf(patient, vcf)
        make_hla(patient, hla)
        make_event_manifest(patient, events)

        # Step 1: Generate peptides
        print(f"  [1/7] Generating peptides...")
        candidates_jsonl = tmp / f"{pid}_candidates.jsonl"
        peptides_fasta = tmp / f"{pid}.fa"
        peptides_tsv = tmp / f"{pid}.tsv"
        _cli(
            "generate-peptides",
            "--sample-id", pid,
            "--annotated-vcf", str(vcf),
            "--event-manifest", str(events),
            "--gtf", str(gtf),
            "--candidates-jsonl", str(candidates_jsonl),
            "--peptides-fasta", str(peptides_fasta),
            "--peptides-tsv", str(peptides_tsv),
        )

        # Step 2: Predict binding
        print(f"  [2/7] Predicting MHC binding...")
        binding_tsv = tmp / f"{pid}_binding.tsv"
        if args.skip_mhcflurry:
            _cli("binding-placeholder", "--tool", "mhcflurry",
                 "--peptides", str(peptides_tsv), "--hla", str(hla),
                 "--output", str(binding_tsv))
        else:
            _cli("predict-mhcflurry", "--peptides", str(peptides_tsv),
                 "--hla", str(hla), "--output", str(binding_tsv))

        # Step 3: VAF extraction (placeholder)
        print(f"  [3/7] Extracting VAF...")
        vaf_tsv = tmp / f"{pid}_vaf.tsv"
        _cli("extract-vaf", "--sample-id", pid,
             "--tumor-bam", "nonexistent.bam",
             "--annotated-vcf", str(vcf),
             "--output", str(vaf_tsv))

        # Step 4: Clonality
        print(f"  [4/7] Computing clonality...")
        clon_tsv = tmp / f"{pid}_clonality.tsv"
        _cli("pyclone-stub", "--vaf-table", str(vaf_tsv),
             "--output", str(clon_tsv))

        # Step 5: HLA LOH (coverage-based, falls back gracefully)
        print(f"  [5/7] Detecting HLA LOH...")
        loh_tsv = tmp / f"{pid}_loh.tsv"
        _cli("detect-loh", "--tumor-bam", "nonexistent.bam",
             "--normal-bam", "nonexistent.bam",
             "--hla", str(hla), "--output", str(loh_tsv))

        # Inject a realistic HLA LOH event for MEL005 to demonstrate immune escape penalty.
        # HLA-A*02:01 loss is common in melanoma (~15% of A*02:01+ patients).
        if pid == "MEL005":
            loh_rows = []
            with open(loh_tsv, newline="") as f:
                reader = csv.DictReader(f, delimiter="\t")
                for row in reader:
                    if row["allele"] == "HLA-A*02:01":
                        row["is_loh"] = "True"
                        row["loh_pvalue"] = "0.0023"
                        row["copy_number"] = "1.1"
                    loh_rows.append(row)
            with open(loh_tsv, "w", newline="") as f:
                writer = csv.DictWriter(f, delimiter="\t", fieldnames=loh_rows[0].keys())
                writer.writeheader()
                writer.writerows(loh_rows)
            print(f"  -> Injected HLA-A*02:01 LOH for {pid} (immune escape penalty demo)")

        # Step 6: Make features
        print(f"  [6/7] Building features...")
        features_tsv = tmp / f"{pid}_features.tsv"
        enriched_jsonl = tmp / f"{pid}_enriched.jsonl"
        _cli("make-features",
             "--candidates", str(candidates_jsonl),
             "--binding", str(binding_tsv),
             "--clonality", str(clon_tsv),
             "--loh", str(loh_tsv),
             "--ms-evidence", str(ROOT / "assets" / "empty_ms_evidence.tsv"),
             "--annotated-vcf", str(vcf),
             "--features", str(features_tsv),
             "--enriched-candidates", str(enriched_jsonl))

        # Step 7: Score with the IEDB-trained XGBoost/LightGBM ensemble.
        print(f"  [7/7] Scoring candidates...")
        scored_jsonl = tmp / f"{pid}_scored.jsonl"
        score_args = ["score", "--features", str(features_tsv),
                      "--candidates", str(enriched_jsonl),
                      "--output", str(scored_jsonl)]
        model_path = ROOT / "models" / "publication_model.pkl"
        if model_path.exists():
            score_args += ["--model-path", str(model_path)]
        _cli(*score_args)

        # Generate patient report
        ranked_tsv = pdir / f"{pid}_ranked.tsv"
        report_html = pdir / f"{pid}_report.html"
        _cli("patient-report", "--sample-id", pid,
             "--candidates", str(scored_jsonl),
             "--ranked-tsv", str(ranked_tsv),
             "--html", str(report_html))
        ranked_tables.append(str(ranked_tsv))

        # Count by variant type
        n_total = sum(1 for _ in _read_candidates_tsv(ranked_tsv))
        n_snv = sum(1 for r in _read_candidates_tsv(ranked_tsv) if r.get("variant_type") == "SNV")
        n_frameshift = sum(1 for r in _read_candidates_tsv(ranked_tsv) if r.get("variant_type") == "frameshift")
        n_fusion = sum(1 for r in _read_candidates_tsv(ranked_tsv) if r.get("variant_type") == "fusion")
        n_splice = sum(1 for r in _read_candidates_tsv(ranked_tsv) if r.get("variant_type") == "splice")
        n_ncorf = sum(1 for r in _read_candidates_tsv(ranked_tsv) if r.get("variant_type") == "non_canonical_orf")
        high_priority = sum(1 for r in _read_candidates_tsv(ranked_tsv) if float(r.get("neoflux_score", "0") or 0) > 0.5)

        summary_rows.append({
            "patient_id": pid,
            "total_candidates": n_total,
            "snv": n_snv,
            "frameshift_indel": n_frameshift,
            "fusion": n_fusion,
            "splice": n_splice,
            "ncorf": n_ncorf,
            "high_priority": high_priority,
        })
        parts = [f"{n_snv} SNV"]
        if n_fusion: parts.append(f"{n_fusion} fusion")
        if n_splice: parts.append(f"{n_splice} splice")
        if n_ncorf: parts.append(f"{n_ncorf} nCORF")
        if n_frameshift: parts.append(f"{n_frameshift} frameshift")
        print(f"  -> {n_total} candidates ({', '.join(parts)}, {high_priority} high-priority)")

    # Cohort report
    print(f"\n{'='*60}")
    print(f"Generating cohort report...")
    print(f"{'='*60}")
    summary_tsv = out / "cohort_summary.tsv"
    report_html = out / "cohort_report.html"
    figures_dir = out / "figures"

    _cli("cohort-report",
         "--ranked-tables", *ranked_tables,
         "--summary", str(summary_tsv),
         "--html", str(report_html),
         "--figures-dir", str(figures_dir))

    # Write quick summary
    summary_tsv2 = out / "patient_summary.tsv"
    with open(summary_tsv2, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=summary_rows[0].keys(), delimiter="\t")
        w.writeheader()
        w.writerows(summary_rows)

    print(f"\n{'='*60}")
    print(f"DONE! Results in: {out.resolve()}")
    print(f"  Patient tables:  {out}/<patient_id>/<patient_id>_ranked.tsv")
    print(f"  Cohort summary:  {summary_tsv}")
    print(f"  Figures:         {figures_dir}/")
    print(f"  HTML report:     {report_html}")
    print(f"{'='*60}")


def _read_candidates_tsv(path):
    """Read candidate TSV and yield rows as dicts."""
    try:
        with open(path) as f:
            reader = csv.DictReader(f, delimiter="\t")
            yield from reader
    except FileNotFoundError:
        return


if __name__ == "__main__":
    main()
