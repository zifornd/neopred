//
// Run Variant Annotation
//

include { PVACTOOLS_INSTALLVEPPLUGIN } from '../../../modules/nf-core/pvactools/installvepplugin/main'
include { ENSEMBLVEP_DOWNLOAD } from '../../../modules/nf-core/ensemblvep/download/main'
include { ENSEMBLVEP_VEP      } from '../../../modules/nf-core/ensemblvep/vep/main'

workflow VARIANT_ANNOTATION {
    take:
    ensemblvep_info
    vcf
    vcf_tbi
    fasta
    vep_genome
    vep_species
    vep_cache_version

    main:
    versions = Channel.empty()

    PVACTOOLS_INSTALLVEPPLUGIN()

    ENSEMBLVEP_DOWNLOAD(ensemblvep_info)
    //vep_cache = ENSEMBLVEP_DOWNLOAD.out.cache

    //vcf_for_vep = vcf.map{ meta, vcf -> [ meta, vcf, [] ] }
    vep_extra_files = []
    vep_cache = ENSEMBLVEP_DOWNLOAD.out.cache.map{ meta, cache -> cache}
    vcf_for_vep = vcf.map{vcf -> [ [], vcf, []]}
    vcf_tbi.view()

    ENSEMBLVEP_VEP(vcf_for_vep, vep_genome, vep_species, vep_cache_version, vep_cache, fasta, vep_extra_files)

    // Gather versions of all tools used
    versions = versions.mix(ENSEMBLVEP_DOWNLOAD.out.versions)
    versions = versions.mix(ENSEMBLVEP_VEP.out.versions)

    emit:
    ensemblvep_cache = ENSEMBLVEP_DOWNLOAD.out.cache.collect()  // channel: [ meta, cache ]
    //snpeff_cache     = SNPEFF_DOWNLOAD.out.cache.collect()      // channel: [ meta, cache ]

    versions // channel: [ versions.yml ]
}
