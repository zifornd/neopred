//
// Run Variant Annotation
//

include { PVACTOOLS_INSTALLVEPPLUGIN } from '../../../modules/local/pvactools/installvepplugin/main'
include { ENSEMBLVEP_DOWNLOAD } from '../../../modules/nf-core/ensemblvep/download/main'
include { ENSEMBLVEP_VEP      } from '../../../modules/nf-core/ensemblvep/vep/main'

workflow VARIANT_ANNOTATION {
    take:
    vcf
    vcf_tbi
    ensemblvep_info
    fasta
    vep_genome
    vep_species
    vep_cache_version

    main:
    versions = Channel.empty()

    PVACTOOLS_INSTALLVEPPLUGIN()
    versions = versions.mix(PVACTOOLS_INSTALLVEPPLUGIN.out.versions)

    ENSEMBLVEP_DOWNLOAD(ensemblvep_info)

    vcf_for_vep = vcf.map{ meta, vcf -> [ meta, vcf, [] ] }
    vep_extra_files = PVACTOOLS_INSTALLVEPPLUGIN.out.results
    vep_cache = ENSEMBLVEP_DOWNLOAD.out.cache.map{ meta, cache -> [cache] }.first()
    //vcf_for_vep = vcf.map{vcf -> [ [], vcf, []]}

    ENSEMBLVEP_VEP(vcf_for_vep, vep_genome, vep_species, vep_cache_version, vep_cache, fasta, vep_extra_files)

    // Gather versions of all tools used
    versions = versions.mix(ENSEMBLVEP_DOWNLOAD.out.versions)
    versions = versions.mix(ENSEMBLVEP_VEP.out.versions)

    emit:

    results = ENSEMBLVEP_VEP.out.vcf
    ensemblvep_cache = ENSEMBLVEP_DOWNLOAD.out.cache.collect()  // channel: [ meta, cache ]

    versions // channel: [ versions.yml ]
}
