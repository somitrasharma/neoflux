// Transcript-aware neoantigen peptide generation.
// Python module hook for SNV, indel, frameshift, fusion, splice, and ncORF
// candidates. Production sequence extraction should use VEP protein context,
// STAR-Fusion predicted proteins, rMATS transcript structures, and an ncORF
// caller such as ORFquant or Ribo-TISH.

process PEPTIDE_GENERATION {
    tag "${meta.id}"

    input:
    tuple val(meta), path(annotated_vcf), path(event_manifest)
    path gtf

    output:
    tuple val(meta), path("${meta.id}.candidates.jsonl"), path("${meta.id}.peptides.fa"), path("${meta.id}.peptides.tsv"), emit: peptides

    script:
    """
    python -m neoflux.cli generate-peptides \\
        --sample-id ${meta.id} \\
        --annotated-vcf ${annotated_vcf} \\
        --event-manifest ${event_manifest} \\
        --gtf ${gtf} \\
        --candidates-jsonl ${meta.id}.candidates.jsonl \\
        --peptides-fasta ${meta.id}.peptides.fa \\
        --peptides-tsv ${meta.id}.peptides.tsv \\
        --mhc-i-lengths ${params.peptide_lengths_i} \\
        --mhc-ii-lengths ${params.peptide_lengths_ii}
    """

    stub:
    """
    touch ${meta.id}.candidates.jsonl
    printf "candidate_id\\tsample_id\\tpeptide\\tvariant_type\\n" > ${meta.id}.peptides.tsv
    touch ${meta.id}.peptides.fa
    """
}
