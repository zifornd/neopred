//
// Run Variant Annotation
//

include { ENSEMBLVEP_DOWNLOAD } from '../../../modules/nf-core/ensemblvep/download/main'
include { ENSEMBLVEP_VEP      } from '../../../modules/nf-core/ensemblvep/vep/main'

workflow VARIANT_ANNOTATION {
    take:
    ensemblvep_info
    snpeff_info

    main:
    versions = Channel.empty()

    ENSEMBLVEP_DOWNLOAD(ensemblvep_info)
    SNPEFF_DOWNLOAD(snpeff_info)

    // Gather versions of all tools used
    versions = versions.mix(ENSEMBLVEP_DOWNLOAD.out.versions)
    versions = versions.mix(SNPEFF_DOWNLOAD.out.versions)

    emit:
    ensemblvep_cache = ENSEMBLVEP_DOWNLOAD.out.cache.collect()  // channel: [ meta, cache ]
    snpeff_cache     = SNPEFF_DOWNLOAD.out.cache.collect()      // channel: [ meta, cache ]

    versions // channel: [ versions.yml ]
}