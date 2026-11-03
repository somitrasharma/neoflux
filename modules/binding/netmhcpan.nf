// Optional class I binding prediction with NetMHCpan.
// NetMHCpan requires an academic license and a local/custom container.

process NETMHCPAN {
    tag "${meta.id}"

    input:
    tuple val(meta), path(candidates_jsonl), path(peptides_fasta), path(peptides_tsv), path(hla_tsv)

    output:
    tuple val(meta), path("${meta.id}.netmhcpan.tsv"), emit: predictions

    script:
    """
    netMHCpan \\
        -BA \\
        -f ${peptides_fasta} \\
        -a \$(python -m neoflux.cli hla-list --input ${hla_tsv} --class I) \\
        -xls \\
        -xlsfile ${meta.id}.netmhcpan.raw.xls \\
        > ${meta.id}.netmhcpan.log

    python -m neoflux.cli normalize-binding \\
        --tool netmhcpan \\
        --raw ${meta.id}.netmhcpan.raw.xls \\
        --peptides ${peptides_tsv} \\
        --hla ${hla_tsv} \\
        --output ${meta.id}.netmhcpan.tsv
    """

    stub:
    """
    printf "candidate_id\\tpeptide\\thla_allele\\tpredictor\\tic50\\tpercent_rank\\tbinding_level\\n" > ${meta.id}.netmhcpan.tsv
    """
}
