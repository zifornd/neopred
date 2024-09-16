process VATOOLS_REFTRANSCRIPTMISMATCHREPORTER {
    tag "$meta.id"
    label "process_medium"

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/griffithlab/vatools:5.1.1' :
        'griffithlab/vatools:5.1.1' }"

    input:
    tuple val(meta), path(vcf)
    val filter

    output:
    tuple val(meta), path("*.vcf")  , optional:true, emit: filter_vcf

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix = task.ext.prefix ?: "${meta.id}.filter"
    """
    ref-transcript-mismatch-reporter \\
        $vcf \\
        --filter ${filter} \\
        -o ${prfix}.vcf

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """
}
