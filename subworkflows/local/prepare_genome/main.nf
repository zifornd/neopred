//include indices

include { STAR_GENOMEGENERATE } from '../../../modules/nf-core/star/genomegenerate/main'


//Workflow

workflow PREPARE_GENOME {
    take:
    fasta                //      file: /path/to/genome.fasta
    gtf                  //      file: /path/to/genome.gtf

    main:
    ch_versions = Channel.empty()
    ch_fasta = Channel.value(file(fasta))
    ch_gtf = Channel.value(file(gtf))

    //use star genome generate
    ch_star_index = Channel.empty()
    ch_star_index = STAR_GENOMEGENERATE ( ch_fasta.map { [ [:], it ] }, ch_gtf.map { [ [:], it ] } ).index.map { it[1] }
    ch_versions   = ch_versions.mix(STAR_GENOMEGENERATE.out.versions)

    emit:
    fasta            = ch_fasta                  // channel: path(genome.fasta)
    gtf              = ch_gtf                    // channel: path(genome.gtf)
    star_index       = ch_star_index             // channel: path(star/index/)
    versions         = ch_versions.ifEmpty(null) // channel: [ versions.yml ]

}