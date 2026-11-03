// Optional pMHC stability prediction with NetMHCstabpan.

process NETMHCSTABPAN {
    tag "${meta.id}"

    input:
    tuple val(meta), path(candidates_jsonl), path(peptides_fasta), path(peptides_tsv), path(hla_tsv)

    output:
    tuple val(meta), path("${meta.id}.netmhcstabpan.tsv"), emit: predictions

    script:
    """
    netMHCstabpan \\
        -f ${peptides_fasta} \\
        -a \$(python -m neoflux.cli hla-list --input ${hla_tsv} --class I) \\
        > ${meta.id}.netmhcstabpan.raw.txt

    python -m neoflux.cli normalize-binding \\
        --tool netmhcstabpan \\
        --raw ${meta.id}.netmhcstabpan.raw.txt \\
        --peptides ${peptides_tsv} \\
        --hla ${hla_tsv} \\
        --output ${meta.id}.netmhcstabpan.tsv
    """

    stub:
    """
    printf "candidate_id\\tpeptide\\thla_allele\\tpredictor\\thalf_life\\tstability_rank\\n" > ${meta.id}.netmhcstabpan.tsv
    """
}
