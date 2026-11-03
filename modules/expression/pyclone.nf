// PyClone-VI clonality inference with fallback to heuristic stub.
//
// Production use (requires conda):
//   conda install -c bioconda pyclone-vi
//
// Then uncomment the --purity parameter and set the pyclone-vi script path.
// Without pyclone-vi installed, falls back to depth-aware heuristic.

process PYCLONE {
    tag "${meta.id}"

    input:
    tuple val(meta), path(vaf_table)

    output:
    tuple val(meta), path("${meta.id}.clonality.tsv"), emit: clusters

    script:
    """
    python ${projectDir}/scripts/run_pyclone_vi.py \\
        --vaf-table ${vaf_table} \\
        --sample-id ${meta.id} \\
        --output ${meta.id}.clonality.tsv

    # Fallback if pyclone-vi not installed:
    # python -m neoflux.cli pyclone-stub \\
    #     --vaf-table ${vaf_table} \\
    #     --output ${meta.id}.clonality.tsv
    """

    stub:
    """
    printf "variant_id\\tclonal_status\\tcancer_cell_fraction\\n" > ${meta.id}.clonality.tsv
    """
}
