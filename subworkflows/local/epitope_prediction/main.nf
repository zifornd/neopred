//
// pVACseq epitope binding predictions
//

include { TX2GENE_PVACSEQ      } from '../../../modules/local/tx2gene_pvacseq'
include { TXIMPORT       } from '../../../modules/local/tximport'
include { BATCH_REMOVAL              } from '../../../modules/local/batch_removal'
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
    hla
    filter_val
    callers
    neoantigen_epitope1_lengths
    fasta             // channel: /path/to/genome.fa
    fasta_fai         // channel: /path/to/genome.fai
    gtf

    main:

    ch_versions = Channel.empty()

    TX2GENE_PVACSEQ(salmon.collect{it[1]}, 'salmon', gtf)
    ch_versions = ch_versions.mix(TX2GENE_PVACSEQ.out.versions.first())

    TXIMPORT(salmon.collect{it[1]}, TX2GENE.out.tsv.collect(), 'salmon')
    ch_versions = ch_versions.mix(TXIMPORT.out.versions.first())

    ch_samplesheet = Channel.value(file(samplesheet))

    BATCH_REMOVAL(ch_samplesheet,batch,design,TXIMPORT.out.tpm_gene)
    ch_versions = ch_versions.mix(VATOOLS_REFTRANSCRIPTMISMATCHREPORTER.out.versions.first())

    VATOOLS_REFTRANSCRIPTMISMATCHREPORTER(vcf,filter_val)
    ch_versions = ch_versions.mix(VATOOLS_REFTRANSCRIPTMISMATCHREPORTER.out.versions.first())

    ch_vcf_filt = VATOOLS_REFTRANSCRIPTMISMATCHREPORTER.out.filter_vcf
        .map{ meta, vcf -> [ meta, vcf, [] ] }

    BCFTOOLS_VIEW(ch_vcf_filt,[],[],[])
    ch_versions = ch_versions.mix(BCFTOOLS_VIEW.out.versions.first())

    VATOOLS_VCFEXPRESSIONANNOTATOR (BCFTOOLS_VIEW.out.vcf, BATCH_REMOVAL.out.after_br)
    ch_versions = ch_versions.mix(VATOOLS_VCFEXPRESSIONANNOTATOR.out.versions.first())

    PVACTOOLS_PVACSEQ(VATOOLS_VCFEXPRESSIONANNOTATOR.out.expr_vcf,hla,callers,neoantigen_epitope1_lengths)
    ch_results = PVACTOOLS_PVACSEQ.out.results
    ch_versions = ch_versions.mix(PVACTOOLS_PVACSEQ.out.versions.first())

    emit:
    results      = ch_results           // channel: [ val(meta), bqsr_table ]
    versions     = ch_versions          // channel: [ versions.yml	]
}
