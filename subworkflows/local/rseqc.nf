// TODO nf-core: If in doubt look at other nf-core/subworkflows to see how we are doing things! :)
//               https://github.com/nf-core/modules/tree/master/subworkflows
//               You can also ask for help via your pull request or on the #subworkflows channel on the nf-core Slack workspace:
//               https://nf-co.re/join
// TODO nf-core: A subworkflow SHOULD import at least two modules

include { RSEQC_TIN                 } from '../modules/nf-core/rseqc/tin/main'
include { RSEQC_READDISTRIBUTION    } from '../modules/nf-core/rseqc/readdistribution/main'
include { RSEQC_JUNCTIONSATURATION  } from '../modules/nf-core/rseqc/junctionsaturation/main'


workflow RSEQC {

    take:
    // TODO nf-core: edit input (take) channels
    ch_bam // channel: [ val(meta), [ bam ] ]

    main:

    ch_versions = Channel.empty()

    // TODO nf-core: substitute modules here for the modules of your subworkflow

    RSEQC_TIN ( ch_bam,
                ch_bed )
    ch_versions = ch_versions.mix(RSEQC_TIN.out.versions.first())

    RSEQC_READDISTRIBUTION ( SAMTOOLS_SORT.out.bam )
    ch_versions = ch_versions.mix(SAMTOOLS_INDEX.out.versions.first())

    emit:
    // TODO nf-core: edit emitted channels
    bam      = SAMTOOLS_SORT.out.bam           // channel: [ val(meta), [ bam ] ]
    bai      = SAMTOOLS_INDEX.out.bai          // channel: [ val(meta), [ bai ] ]
    csi      = SAMTOOLS_INDEX.out.csi          // channel: [ val(meta), [ csi ] ]

    versions = ch_versions                     // channel: [ versions.yml ]
}

