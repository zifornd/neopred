//include indices

include { STAR_GENOMEGENERATE               } from '../../../modules/nf-core/star/genomegenerate/main'
include { GUNZIP as GUNZIP_FASTA            } from '../../../modules/nf-core/gunzip'
include { GUNZIP as GUNZIP_GTF              } from '../../../modules/nf-core/gunzip'
include { GUNZIP as GUNZIP_GFF              } from '../../../modules/nf-core/gunzip'
include { GFFREAD                           } from '../../../modules/nf-core/gffread'
include { GUNZIP as GUNZIP_TRANSCRIPT_FASTA } from '../../../modules/nf-core/gunzip'
include { SAMTOOLS_FAIDX                            } from '../../../modules/nf-core/samtools/faidx/main'
include { RSEM_PREPAREREFERENCE as MAKE_TRANSCRIPTS_FASTA       } from '../../../modules/nf-core/rsem/preparereference'

//Workflow

workflow PREPARE_GENOME {
    take:
    fasta                //      file: /path/to/genome.fasta
    gtf                  //      file: /path/to/genome.gtf
    transcript_fasta     //      file: /path/to/transcript.fasta

    main:
    ch_versions = Channel.empty()
    //
    // Uncompress genome fasta file if required
    //
    if (fasta.endsWith('.gz')) {
        ch_fasta    = GUNZIP_FASTA ( [ [:], fasta ] ).gunzip.map { it[1] }
        ch_versions = ch_versions.mix(GUNZIP_FASTA.out.versions)
    } else {
        ch_fasta = Channel.value(file(fasta))
    }

    //
    // Uncompress GTF annotation file or create from GFF3 if required
    //
    if (gtf || gff) {
        if (gtf) {
            if (gtf.endsWith('.gz')) {
                ch_gtf      = GUNZIP_GTF ( [ [:], gtf ] ).gunzip.map { it[1] }
                ch_versions = ch_versions.mix(GUNZIP_GTF.out.versions)
            } else {
                ch_gtf = Channel.value(file(gtf))
            }
        } else if (gff) {
            if (gff.endsWith('.gz')) {
                ch_gff      = GUNZIP_GFF ( [ [:], gff ] ).gunzip.map { it[1] }
                ch_versions = ch_versions.mix(GUNZIP_GFF.out.versions)
            } else {
                ch_gff = Channel.value(file(gff))
            }
            ch_gtf      = GFFREAD ( ch_gff ).gtf
            ch_versions = ch_versions.mix(GFFREAD.out.versions)
        }
    }

    //
    // Uncompress transcript fasta file / create if required
    //
    if (transcript_fasta) {
        if (transcript_fasta.endsWith('.gz')) {
            ch_transcript_fasta = GUNZIP_TRANSCRIPT_FASTA ( [ [:], transcript_fasta ] ).gunzip.map { it[1] }
            ch_versions         = ch_versions.mix(GUNZIP_TRANSCRIPT_FASTA.out.versions)
        } else {
            ch_transcript_fasta = Channel.value(file(transcript_fasta))
        }
    } else {
        ch_transcript_fasta = MAKE_TRANSCRIPTS_FASTA ( ch_fasta, ch_gtf ).transcript_fasta
        ch_versions         = ch_versions.mix(MAKE_TRANSCRIPTS_FASTA.out.versions)
    }

    //use star genome generate
    ch_star_index = Channel.empty()
    ch_star_index = STAR_GENOMEGENERATE ( ch_fasta.map { [ [:], it ] }, ch_gtf.map { [ [:], it ] } ).index.map { it[1] }
    ch_versions   = ch_versions.mix(STAR_GENOMEGENERATE.out.versions)

    SAMTOOLS_FAIDX(ch_fasta.map { [ [:], it ] }, [[:],[]] )
    ch_fasta_fai   = SAMTOOLS_FAIDX.out.fai.map{ meta, fai -> [fai] }
    ch_versions   = ch_versions.mix(SAMTOOLS_FAIDX.out.versions)

    emit:
    fasta            = ch_fasta                  // channel: path(genome.fasta)
    gtf              = ch_gtf                    // channel: path(genome.gtf)
    star_index       = ch_star_index             // channel: path(star/index/)
    fasta_fai        = ch_fasta_fai              // channel: path(genome.fasta.fai)
    transcript_fasta = ch_transcript_fasta       // channel: path(transcript.fasta)
    versions         = ch_versions.ifEmpty(null) // channel: [ versions.yml ]
}
