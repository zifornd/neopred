process ARCASHLA_CONVERT {
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/arcas-hla:0.5.0--hdfd78af_0':
        'docker.io/jfx319/arcashla:latest' }"

    input:
    path(gt)

    output:
    path ("*.p-group.tsv")   , emit: gt_group
    path "versions.yml"               , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def VERSION = "0.5.0" // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.

    """
    if [[ ${workflow.profile} == *"conda"* ]]; then
    arcasHLA reference --version 3.24.0
    fi
    arcasHLA \\
        convert \\
        -r p-group \\
        -o ./genotypes.p-group.tsv \\
        $gt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        arcashla: $VERSION
    END_VERSIONS
    """
}
