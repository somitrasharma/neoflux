from neoflux.models import GenomicVariant, VariantType
from neoflux.peptide_generator import PeptideGenerator
from neoflux.vcf_peptide import (
    generate_all_candidates,
    generate_from_vcf_record,
    write_peptide_exports,
)


def test_snv_generator_emits_mutant_only_windows():
    generator = PeptideGenerator(mhc_i_lengths=(9,), mhc_ii_lengths=())
    variant = GenomicVariant(
        sample_id="S1",
        chromosome="chr7",
        position=1,
        ref_allele="A",
        alt_allele="T",
        variant_type=VariantType.SNV,
        gene_symbol="BRAF",
    )
    mutant = "AAAAAVVVVAAAAA"
    wildtype = "AAAAAAVVVVAAAAA"
    windows = generator.generate_from_snv(
        variant=variant,
        mutant_protein=mutant,
        wildtype_protein=wildtype,
        mutation_pos=5,
    )

    assert windows
    assert all(window.peptide != window.wildtype_peptide for window in windows)
    assert all(len(window.peptide) == 9 for window in windows)


def test_fusion_generator_spans_junction():
    generator = PeptideGenerator(mhc_i_lengths=(9,), mhc_ii_lengths=())
    variant = GenomicVariant(
        sample_id="S1",
        chromosome="chr2",
        position=1,
        ref_allele="N",
        alt_allele="FUSION",
        variant_type=VariantType.FUSION,
        fusion_partner_5p="EML4",
        fusion_partner_3p="ALK",
    )
    fusion_protein = "ACDEFGHIKLMNPQRSTVWYACDEFGHIK"
    junction_pos = 12
    windows = generator.generate_from_fusion(
        variant=variant,
        fusion_protein=fusion_protein,
        junction_pos=junction_pos,
    )

    assert windows
    assert all(window.wildtype_peptide is None for window in windows)


def test_vcf_and_manifest_generation_writes_exports(tmp_path):
    vcf_path = tmp_path / "demo.vcf"
    vcf_path.write_text(
        "##fileformat=VCFv4.2\n"
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
        "chr7\t140453136\t.\tA\tT\t.\tPASS\tGENE=BRAF;HGVS_P=p.V600E\n"
    )
    manifest_path = tmp_path / "events.json"
    manifest_path.write_text(
        '{"fusions":[{"gene_5p":"EML4","gene_3p":"ALK","junction_pos":35}],'
        '"splice_events":[{"event_type":"SE","gene_symbol":"TP53","delta_psi":"0.4"}],'
        '"noncanonical_orfs":[{"orf_id":"nc1","protein":"MTEYKLVVVGAGGVGKSALTIQLIQNHFVDEYDPTIEDSYRKQVVIDGETCLLDILDTAGQEEY"}]}'
    )

    def _records():
        for line in vcf_path.read_text().splitlines():
            if line.startswith("#") or not line.strip():
                continue
            fields = line.split("\t")
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

    candidates = generate_all_candidates(
        sample_id="S1",
        vcf_records=_records(),
        manifest_path=manifest_path,
        mhc_i_lengths=(8, 9, 10, 11),
        mhc_ii_lengths=(15,),
    )
    assert candidates
    assert any(candidate.peptide_source == VariantType.SNV for candidate in candidates)
    assert any(candidate.peptide_source == VariantType.FUSION for candidate in candidates)
    assert any(candidate.peptide_source == VariantType.SPLICE for candidate in candidates)
    assert any(
        candidate.peptide_source == VariantType.NON_CANONICAL_ORF
        for candidate in candidates
    )

    jsonl = tmp_path / "candidates.jsonl"
    fasta = tmp_path / "peptides.fa"
    tsv = tmp_path / "peptides.tsv"
    write_peptide_exports(candidates, jsonl, fasta, tsv)
    assert jsonl.stat().st_size > 0
    assert fasta.stat().st_size > 0
    assert "candidate_id" in tsv.read_text()


def test_generate_from_vcf_record_uses_hgvs_annotation():
    generator = PeptideGenerator(mhc_i_lengths=(9,), mhc_ii_lengths=())
    record = {
        "chrom": "chr7",
        "pos": "140453136",
        "ref": "A",
        "alt": "T",
        "info": "GENE=BRAF;HGVS_P=p.V600E",
    }
    candidates = generate_from_vcf_record(generator, "S1", record)
    assert candidates
