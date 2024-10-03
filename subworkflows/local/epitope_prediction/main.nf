//
// pVACseq epitope binding predictions
//

include { TX2GENE_PVACSEQ      } from '../../../modules/local/tx2gene_pvacseq'
include { TXIMPORT       } from '../../../modules/local/tximport'
include { BATCH_REMOVAL              } from '../../../modules/local/batch_removal'
include { GUNZIP             } from '../../../modules/nf-core/gunzip'
include { VATOOLS_REFTRANSCRIPTMISMATCHREPORTER             } from '../../../modules/local/vatools/reftranscriptmismatchreporter/main'
include { BCFTOOLS_VIEW           } from '../../../modules/nf-core/bcftools/view/main'
include { VATOOLS_VCFEXPRESSIONANNOTATOR                    } from '../../../modules/local/vatools/vcfexpressionannotator/main'
include { PVACTOOLS_PVACSEQ                               } from '../../../modules/local/pvactools/pvacseq/main'
include { SAMTOOLS_FLAGSTAT as SAMTOOLS_FLAGSTAT_BAMRG      } from '../../../modules/nf-core/samtools/flagstat/main'
include { SAMTOOLS_FLAGSTAT as SAMTOOLS_FLAGSTAT_BAMMARKDUP } from '../../../modules/nf-core/samtools/flagstat/main'

workflow EPITOPE_PREDICTION{

    take:
    salmon
    samplesheet
    design
    batch
    vcf               // channel: [ val(meta), [ vcf ] ]
    callers
    neoantigen_epitope1_lengths
    fasta             // channel: /path/to/genome.fa
    fasta_fai         // channel: /path/to/genome.fai
    gtf

    main:

    ch_versions = Channel.empty()

    TX2GENE_PVACSEQ(salmon.collect{it[1]}, 'salmon', gtf)
    ch_versions = ch_versions.mix(TX2GENE_PVACSEQ.out.versions.first())

    TXIMPORT(salmon.collect{it[1]}, TX2GENE_PVACSEQ.out.tsv.collect(), 'salmon')
    ch_versions = ch_versions.mix(TXIMPORT.out.versions.first())

    ch_samplesheet = Channel.value(file(samplesheet))

    ch_tpm = BATCH_REMOVAL(ch_samplesheet,batch,design,TXIMPORT.out.tpm_gene).after_br
    ch_versions = ch_versions.mix(BATCH_REMOVAL.out.versions.first())

    vcf.map{ meta,hla,vcf -> [meta,vcf]}.set{ch_vcf_res}
    vcf.map{ meta,hla,vcf -> [meta,hla]}.set{ch_hla_res}

    ch_vcf     = GUNZIP ( ch_vcf_res ).gunzip
    ch_versions = ch_versions.mix(GUNZIP.out.versions.first())

    VATOOLS_REFTRANSCRIPTMISMATCHREPORTER(ch_vcf)
    ch_versions = ch_versions.mix(VATOOLS_REFTRANSCRIPTMISMATCHREPORTER.out.versions.first())

    ch_vcf_filt = VATOOLS_REFTRANSCRIPTMISMATCHREPORTER.out.filter_vcf
        .map{ meta, vcf -> [ meta, vcf, [] ] }

    ch_in_vcfexp = BCFTOOLS_VIEW(ch_vcf_filt,[],[],[]).vcf
    ch_versions = ch_versions.mix(BCFTOOLS_VIEW.out.versions.first())

    VATOOLS_VCFEXPRESSIONANNOTATOR (ch_in_vcfexp, ch_tpm)
    ch_versions = ch_versions.mix(VATOOLS_VCFEXPRESSIONANNOTATOR.out.versions.first())

    PVACTOOLS_PVACSEQ(VATOOLS_VCFEXPRESSIONANNOTATOR.out.expr_vcf,ch_hla_res,callers,neoantigen_epitope1_lengths)
    ch_results = PVACTOOLS_PVACSEQ.out.results
    ch_versions = ch_versions.mix(PVACTOOLS_PVACSEQ.out.versions.first())

    emit:
    results      = ch_results           // channel: [ val(meta), bqsr_table ]
    versions     = ch_versions          // channel: [ versions.yml	]
}
