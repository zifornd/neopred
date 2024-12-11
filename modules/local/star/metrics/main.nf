process STAR_METRICS {

    conda "perl=5.32.1=7_hd590300_perl5"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/perl:5.32' :
        'quay.io/biocontainers/perl:5.32' }"

    input:
    path(star_log)

    output:
    path("*.csv"), emit: csv
    path "versions.yml"           , emit: versions
	
    when:
    task.ext.when == null || task.ext.when

    script:

    """
    STAR_reports.pl -f $star_log 1> star_metrics.csv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        perl: \$(echo \$(perl --version 2>&1) | sed 's/.*v\\(.*\\)) built.*/\\1/')
    END_VERSIONS
    """

stub:

     """
   touch star_metrics.csv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        perl: \$(echo \$(perl --version 2>&1) | sed 's/.*v\\(.*\\)) built.*/\\1/')
    END_VERSIONS
    """
}
