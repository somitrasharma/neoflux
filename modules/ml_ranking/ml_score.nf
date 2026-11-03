// Final NeoFlux ML ranking. Uses a trained model when provided, otherwise a
// transparent baseline scorer for development and dry-run benchmarking.

process ML_SCORE {
    tag "${meta.id}"

    input:
    tuple val(meta), path(feature_matrix), path(candidates_jsonl)
    val model_path

    output:
    tuple val(meta), path("${meta.id}.scored_candidates.jsonl"), emit: scored_candidates

    script:
    """
    python -m neoflux.cli score \\
        --features ${feature_matrix} \\
        --candidates ${candidates_jsonl} \\
        --model-path ${model_path ?: ""} \\
        --output ${meta.id}.scored_candidates.jsonl
    """

    stub:
    """
    cp ${candidates_jsonl} ${meta.id}.scored_candidates.jsonl
    """
}
