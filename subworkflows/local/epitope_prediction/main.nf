//
// pVACseq epitope binding predictions
//

include { VATOOLS_REFTRANSCRIPTMISMATCHREPORTER             } from '../../../modules/local/vatools/reftranscriptmismatchreporter/main'
include { BCFTOOLS_VIEW           } from '../../../modules/nf-core/bcftools/view/main'
include { VATOOLS_VCFEXPRESSIONANNOTATOR                    } from '../../../modules/local/vatools/vcfexpressionannotator/main'
include { PVACTOOLS_PVACSEQ                             } from '../../../modules/local/pvactools/pvacseq/main'
include { SAMTOOLS_FLAGSTAT as SAMTOOLS_FLAGSTAT_BAMRG      } from '../../../modules/nf-core/samtools/flagstat/main'
include { SAMTOOLS_FLAGSTAT as SAMTOOLS_FLAGSTAT_BAMMARKDUP } from '../../../modules/nf-core/samtools/flagstat/main'

workflow EPITOPE_PREDICTION{

    take:
    vcf               // channel: [ val(meta), [ vcf ] ]
    filter_val
    csv           // channel: [ val(meta), [ csv_expression ] ]
    callers
    fasta             // channel: /path/to/genome.fa
    fasta_fai         // channel: /path/to/genome.fai

    main:

    ch_versions = Channel.empty()

    VATOOLS_REFTRANSCRIPTMISMATCHREPORTER(vcf,filter_val)
    ch_versions = ch_versions.mix(VATOOLS_REFTRANSCRIPTMISMATCHREPORTER.out.versions.first())

    ch_vcf_filt = VATOOLS_REFTRANSCRIPTMISMATCHREPORTER.out.filter_vcf
        .map{ meta, vcf -> [ meta, vcf, [] ] }

    BCFTOOLS_VIEW(ch_vcf_filt,[],[],[])
    ch_versions = ch_versions.mix(BCFTOOLS_VIEW.out.versions.first())

    VATOOLS_VCFEXPRESSIONANNOTATOR (BCFTOOLS_VIEW.out.vcf, csv)
    ch_versions = ch_versions.mix(VATOOLS_VCFEXPRESSIONANNOTATOR.out.versions.first())

    PVACTOOLS_PVACSEQ(VATOOLS_VCFEXPRESSIONANNOTATOR.out.expr_vcf,callers)
    ch_results = PVACTOOLS_PVACSEQ.out.results
    ch_versions = ch_versions.mix(PVACTOOLS_PVACSEQ.out.versions.first())

    emit:
    results      = ch_results           // channel: [ val(meta), bqsr_table ]
    versions     = ch_versions          // channel: [ versions.yml	]
}
