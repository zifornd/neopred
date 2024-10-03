process PVACTOOLS_PVACSEQ {
    tag "$meta.id"
    label "process_high"

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/griffithlab/pvactools:4.1.1' :
        'docker.io/griffithlab/pvactools:4.1.1' }"

    input:
    tuple val(meta), path(vcf)
    tuple val(meta_1), val(hla)
    val algorithms
    val neoantigen_epitope1_lengths

    output:
    tuple val(meta), path("*.{tsv,R}")  , optional:true, emit: results
    path "versions.yml"            , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix = task.ext.prefix ?: "${meta.id}r"
    def callers = algorithms.split(',').join(' ') ?: "MHCflurry NetMHCcons MHCnuggetsII"
    def alleles = hla.join(',')
    """
    pvacseq run \\
        $vcf ${meta.id} $alleles $callers ./ \\
        -e1 ${neoantigen_epitope1_lengths} \\
        -t $task.cpus \\
        -r 5 --pass-only \\
        --fasta-size 200 \\
        -d 250 \\
        --iedb-install-directory /opt/iedb

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """
}
