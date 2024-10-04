process GATK4_COUNTVARIANTS {
    tag "${meta.id}"
    label 'process_single'

    conda "bioconda::gatk4=4.3.0.0"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/gatk4:4.3.0.0--py36hdfd78af_0':
        'quay.io/biocontainers/gatk4:4.3.0.0--py36hdfd78af_0' }"

    input:
    tuple val(meta), path(vcf)     // channel: [ val(meta), vcf files]
    tuple val(meta), path(tbi)     // channel: [ val(meta), tbi file]

    output:
    tuple val(meta), path("*_counts"), emit: counts
    path "versions.yml"              , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    // Sample metaid
    def prefix = task.ext.prefix ?: "${meta.id}_counts"

    """
    gatk CountVariants \\
        -V $vcf \\
        -O ${prefix}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gatk4: \$(echo \$(gatk --version 2>&1) | sed 's/^.*(GATK) v//; s/ .*\$//')
    END_VERSIONS
    """
}
