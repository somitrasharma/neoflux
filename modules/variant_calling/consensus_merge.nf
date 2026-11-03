// Merge SNV/indel calls plus unconventional event manifests into one contract.

process CONSENSUS_MERGE {
    tag "${meta.id}"

    input:
    tuple val(meta), path(mutect_vcf), path(strelka_vcf), path(fusions), path(splice_events)

    output:
    tuple val(meta), path("${meta.id}.consensus.vcf"), path("${meta.id}.events.json"), emit: merged

    script:
    """
    python -m neoflux.cli merge-variants \\
        --sample-id ${meta.id} \\
        --mutect-vcf ${mutect_vcf} \\
        --strelka-vcf ${strelka_vcf} \\
        --fusions ${fusions} \\
        --splice-events ${splice_events} \\
        --output-vcf ${meta.id}.consensus.vcf \\
        --event-manifest ${meta.id}.events.json
    """

    stub:
    """
    printf "##fileformat=VCFv4.2\\n#CHROM\\tPOS\\tID\\tREF\\tALT\\tQUAL\\tFILTER\\tINFO\\nchr7\\t140453136\\t.\\tA\\tT\\t.\\tPASS\\tGENE=BRAF;HGVS_P=p.V600E;NEOFLUX_CALLERS=mutect2,strelka2\\n" > ${meta.id}.consensus.vcf
    cat > ${meta.id}.events.json <<'EOF'
    {
      "sample_id": "${meta.id}",
      "fusions": [
        {
          "gene_5p": "EML4",
          "gene_3p": "ALK",
          "junction_pos": 35,
          "chrom": "chr2",
          "pos": 29446394
        }
      ],
      "splice_events": [
        {
          "event_id": "SE_001",
          "event_type": "SE",
          "gene_symbol": "TP53",
          "delta_psi": "0.42",
          "divergence_pos": 25,
          "chrom": "chr17",
          "pos": 7675088
        }
      ],
      "noncanonical_orfs": [
        {
          "orf_id": "ncORF_001",
          "orf_type": "uORF",
          "protein": "MTEYKLVVVGAGGVGKSALTIQLIQNHFVDEYDPTIEDSYRKQVVIDGETCLLDILDTAGQEEY"
        }
      ]
    }
    EOF
    """
}
