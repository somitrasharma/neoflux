import json
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent


def _cli(*args):
    result = subprocess.run(
        [sys.executable, "-m", "neoflux.cli"] + [str(a) for a in args],
        cwd=ROOT,
        capture_output=True,
        text=True,
    )
    assert result.returncode == 0, f"CLI failed: {result.stderr}"
    return result


def test_parse_hla_from_allele_list(tmp_path):
    tsv = tmp_path / "alleles.tsv"
    tsv.write_text("allele\nHLA-A*02:01\nHLA-B*07:02\n")
    out = tmp_path / "parsed.tsv"
    _cli("parse-hla", "--input", str(tsv), "--format", "allele-list", "--output", str(out))
    lines = out.read_text().splitlines()
    assert len(lines) == 3
    assert "HLA-A*02:01" in lines[1]


def test_generate_peptides_produces_all_variant_types(tmp_path):
    vcf = tmp_path / "demo.vcf"
    vcf.write_text(
        "##fileformat=VCFv4.2\n"
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
        "chr7\t140453136\t.\tA\tT\t.\tPASS\tGENE=BRAF;HGVS_P=p.V600E\n"
    )
    events = tmp_path / "events.json"
    events.write_text(
        json.dumps({
            "fusions": [
                {"gene_5p": "EML4", "gene_3p": "ALK",
                 "junction_pos": 35, "chrom": "chr2", "pos": 29446394}
            ],
            "splice_events": [
                {"event_type": "SE", "gene_symbol": "TP53",
                 "delta_psi": "0.4", "divergence_pos": 25}
            ],
            "noncanonical_orfs": [
                {"orf_id": "nc1",
                 "protein": "MTEYKLVVVGAGGVGKSALTIQLIQNHFVDEYDPTIEDSYRKQVVIDGETCLLDILDTAGQEEY"}
            ],
        })
    )
    jsonl = tmp_path / "candidates.jsonl"
    fasta = tmp_path / "peptides.fa"
    tsv = tmp_path / "peptides.tsv"
    _cli(
        "generate-peptides",
        "--sample-id", "S1",
        "--annotated-vcf", str(vcf),
        "--event-manifest", str(events),
        "--gtf", str(ROOT / "assets" / "ref.gtf"),
        "--candidates-jsonl", str(jsonl),
        "--peptides-fasta", str(fasta),
        "--peptides-tsv", str(tsv),
    )
    assert jsonl.stat().st_size > 0
    assert fasta.stat().st_size > 0
    assert tsv.stat().st_size > 0
    candidates = [json.loads(line) for line in jsonl.read_text().splitlines() if line.strip()]
    sources = {c["peptide_source"] for c in candidates}
    assert "SNV" in sources
    assert "fusion" in sources
    assert "splice" in sources
    assert "non_canonical_orf" in sources


def test_end_to_end_heuristic_ranking_produces_scores_and_ranks(tmp_path):
    vcf = tmp_path / "demo.vcf"
    vcf.write_text(
        "##fileformat=VCFv4.2\n"
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
        "chr7\t140453136\t.\tA\tT\t.\tPASS\tGENE=BRAF;HGVS_P=p.V600E\n"
    )
    events = tmp_path / "events.json"
    events.write_text(
        json.dumps({
            "fusions": [],
            "splice_events": [],
            "noncanonical_orfs": [],
        })
    )
    hla = ROOT / "assets" / "hla_alleles.tsv"
    c_jsonl = tmp_path / "candidates.jsonl"
    fasta = tmp_path / "peptides.fa"
    peptides = tmp_path / "peptides.tsv"
    binding = tmp_path / "binding.tsv"
    vaf = tmp_path / "vaf.tsv"
    clon = tmp_path / "clonality.tsv"
    loh = tmp_path / "loh.tsv"
    features = tmp_path / "features.tsv"
    enriched = tmp_path / "enriched.jsonl"
    scored = tmp_path / "scored.jsonl"
    ranked = tmp_path / "ranked.tsv"

    _cli(
        "generate-peptides", "--sample-id", "S1",
        "--annotated-vcf", str(vcf),
        "--event-manifest", str(events),
        "--gtf", str(ROOT / "assets" / "ref.gtf"),
        "--candidates-jsonl", str(c_jsonl),
        "--peptides-fasta", str(fasta),
        "--peptides-tsv", str(peptides),
    )
    _cli(
        "binding-placeholder", "--tool", "mhcflurry",
        "--peptides", str(peptides), "--hla", str(hla),
        "--output", str(binding),
    )
    _cli(
        "extract-vaf", "--sample-id", "S1",
        "--tumor-bam", "assets/tumor.bam",
        "--annotated-vcf", str(vcf), "--output", str(vaf),
    )
    _cli(
        "pyclone-stub", "--vaf-table", str(vaf), "--output", str(clon),
    )
    _cli(
        "lohhla-template", "--sample-id", "S1",
        "--hla", str(hla), "--output", str(loh),
    )
    _cli(
        "make-features", "--candidates", str(c_jsonl),
        "--binding", str(binding), "--clonality", str(clon),
        "--loh", str(loh),
        "--ms-evidence", str(ROOT / "assets" / "empty_ms_evidence.tsv"),
        "--annotated-vcf", str(vcf), "--features", str(features),
        "--enriched-candidates", str(enriched),
    )
    _cli(
        "score", "--features", str(features),
        "--candidates", str(enriched), "--output", str(scored),
    )
    _cli(
        "patient-report", "--sample-id", "S1",
        "--candidates", str(scored), "--ranked-tsv", str(ranked),
        "--html", str(tmp_path / "report.html"),
    )

    lines = ranked.read_text().splitlines()
    assert len(lines) > 1
    assert "neoflux_score" in lines[0]
    assert "neoflux_rank" in lines[0]


def test_cohort_report_generates_summary_and_figures(tmp_path):
    ranked = tmp_path / "ranked.tsv"
    ranked.write_text("candidate_id\tsample_id\tneoflux_score\tneoflux_rank\nNF-001\tS1\t0.85\t1\nNF-002\tS1\t0.42\t2\n")
    _cli(
        "cohort-report",
        "--ranked-tables", str(ranked),
        "--summary", str(tmp_path / "summary.tsv"),
        "--html", str(tmp_path / "report.html"),
        "--figures-dir", str(tmp_path / "figures"),
    )
    assert (tmp_path / "summary.tsv").exists()
    assert (tmp_path / "figures").is_dir()
