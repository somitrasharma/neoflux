// Optional class II binding prediction with NetMHCIIpan.
// NetMHCIIpan requires an academic license and a local/custom container.

process NETMHCIIPAN {
    tag "${meta.id}"

    input:
    tuple val(meta), path(candidates_jsonl), path(peptides_fasta), path(peptides_tsv), path(hla_tsv)

    output:
    tuple val(meta), path("${meta.id}.netmhciipan.tsv"), emit: predictions

    script:
    """
    netMHCIIpan \\
        -f ${peptides_fasta} \\
        -a \$(python -m neoflux.cli hla-list --input ${hla_tsv} --class II) \\
        -xls \\
        -xlsfile ${meta.id}.netmhciipan.raw.xls \\
        > ${meta.id}.netmhciipan.log

    python -m neoflux.cli normalize-binding \\
        --tool netmhciipan \\
        --raw ${meta.id}.netmhciipan.raw.xls \\
        --peptides ${peptides_tsv} \\
        --hla ${hla_tsv} \\
        --output ${meta.id}.netmhciipan.tsv
    """

    stub:
    """
    printf "candidate_id\\tpeptide\\thla_allele\\tpredictor\\tic50\\tpercent_rank\\tbinding_level\\n" > ${meta.id}.netmhciipan.tsv
    """
}
