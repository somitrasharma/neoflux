// Discovery immunopeptidomics search with FragPipe/MSFragger.
// Expects centroided mzML and a generated or project-specific FASTA database.

process FRAGPIPE {
    tag "${meta.id}"

    input:
    tuple val(meta), path(mzml), path(ms_search_results)

    output:
    tuple val(meta), path("${meta.id}.fragpipe_peptides.tsv"), emit: search_results

    script:
    """
    if [ -s "${ms_search_results}" ]; then
        cp ${ms_search_results} ${meta.id}.fragpipe_peptides.tsv
    else
        fragpipe \\
            --headless \\
            --workflow immunopeptidomics.workflow \\
            --manifest ${mzml} \\
            --workdir fragpipe_${meta.id}

        cp fragpipe_${meta.id}/combined_peptide.tsv ${meta.id}.fragpipe_peptides.tsv
    fi
    """

    stub:
    """
    printf "Peptide\\tSpectralCount\\tProbability\\tBestScore\\n" > ${meta.id}.fragpipe_peptides.tsv
    """
}
