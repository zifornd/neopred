process ARCASHLA_GENOTYPE {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/arcas-hla:0.5.0--hdfd78af_0':
        'docker.io/jfx319/arcashla:latest' }"
    input:
    tuple val(meta), path(read)

    output:
    tuple val(meta), path ("*.alignment.p")     , emit: alignment_p
    tuple val(meta), path ("*.genotype.log")    , emit: gt_log
    tuple val(meta), path ("*.genotype.json")   , emit: gt_json
    tuple val(meta), path ("*.genes.json")      , emit: genes_json
    path "versions.yml"                         , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def VERSION = "0.5.0" // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.

    """
    arcasHLA \\
        genotype \\
        $read \\
        -g A,B,C,DPB1,DQB1,DQA1,DRB1 \\
        -o . \\
        -t $task.cpus \\
        -v

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        arcashla: $VERSION
    END_VERSIONS
    """
}
