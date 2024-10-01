process VATOOLS_VCFEXPRESSIONANNOTATOR {
    tag "$meta.id"
    label "process_medium"

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/griffithlab/vatools:5.1.1' :
        'docker.io/griffithlab/vatools:5.1.1' }"

    input:
    tuple val(meta), path(vcf)
    path csv

    output:
    tuple val(meta), path("*.vcf")  , optional:true, emit: expr_vcf
    path "versions.yml"            , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix = task.ext.prefix ?: "${meta.id}.mutect2.somatic.base.snp.Somatic.hc.filter.vep.gx"
    """
    tr ',' '\\t' < $csv | awk 'BEGIN {OFS="\\t"} \$1 != "Gene_ID" {gsub(/\\.[0-9]+/, "", \$1)} 1' > ${meta.id}.tmp
    vcf-expression-annotator \\
        $vcf ${meta.id}.tmp \\
        custom gene --id-column Gene_ID --expression-column ${meta.id} \\
        -s ${meta.id} \\
        -o ${prefix}.vcf
    rm ${meta.id}.tmp

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """
}
