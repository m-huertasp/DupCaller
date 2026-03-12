#!/usr/bin/env nextflow

/*
 * DupCaller Nextflow Pipeline
 *
 * A pipeline for calling somatic mutations from barcoded error-corrected NGS data
 *
 */

nextflow.enable.dsl=2

/*
 * Pipeline parameters with default values
 */
params.help = false
params.outdir = "./results"

// Input files
params.sample_name = null
params.read1 = null
params.read2 = null
params.normal_read1 = null
params.normal_read2 = null
params.normal_bam = null
params.samplesheet = null

// Reference files
params.reference = null
params.germline_vcf = null
params.noise_mask = null
params.indel_bed = null

// Barcode trimming parameters
params.barcode_pattern = "NNNXXXX"

// Alignment parameters
params.threads = 4
params.read_group_id = null
params.platform = "ILLUMINA"

// DupCaller parameters
params.regions = "chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22 chrX chrY"
params.region_file = null
params.max_af = 1.0
params.germline_af_cutoff = 0.001
params.min_n_depth = 10
params.max_zero_qual_frac = 0.5
params.trim_template = 7
params.trim_read = 7
params.mapq = 40
params.window_size = 100000

// Burden estimation parameters
params.estimate_clonal = false
params.estimate_dilute = false
params.gene_bed = null

// Skip options
params.skip_trim = false
params.skip_align = false
params.skip_normal = false

def helpMessage() {
    log.info"""
    ================================================================
    DupCaller Pipeline
    ================================================================

    Usage:
    nextflow run DupCaller.nf --sample_name SAMPLE --read1 R1.fq.gz --read2 R2.fq.gz --reference ref.fa [options]

    Required arguments:
      --sample_name         Sample name for output files
      --reference           Reference genome fasta file (must be indexed with DupCaller.py index)

    Input options (choose one):
      --read1               Read 1 fastq file (paired with --read2)
      --read2               Read 2 fastq file (paired with --read1)
      --skip_trim           Skip trimming step if BAM already exists
      --skip_align          Skip alignment step if BAM already exists

    Normal sample:
      --normal_read1        Normal sample read 1 fastq
      --normal_read2        Normal sample read 2 fastq
      --normal_bam          Pre-aligned normal BAM file
      --skip_normal         Run without matched normal (requires --max_af)

    Reference resources:
      --germline_vcf        Indexed germline VCF with AF field
      --noise_mask          BED file(s) for noise masking (space-separated if multiple)
      --indel_bed           Enhanced panel of normal for indels

    Barcode trimming:
      --barcode_pattern     Barcode pattern (default: NNNXXXX for NanoSeq)

    Alignment:
      --threads             Number of threads (default: 4)
      --platform            Sequencing platform (default: ILLUMINA)

    Variant calling:
      --regions             Contigs for variant calling (default: chr1-22,X,Y)
      --max_af              Maximum allele frequency (default: 1.0, set to 0.1 for panel without normal)
      --germline_af_cutoff  Germline AF cutoff (default: 0.001)
      --min_n_depth         Minimum normal depth (default: 10)
      --max_zero_qual_frac  Maximum zero quality fraction (default: 0.5)
      --trim_template       Template end trim distance (default: 7)
      --trim_read           Read end trim distance (default: 7)
      --mapq                Minimum mapping quality (default: 40)

    Burden estimation:
      --estimate_clonal     Consider clonal mutations (default: false)
      --estimate_dilute     Dilute mode for same starting material (default: false)
      --gene_bed            Gene BED file for coverage calculation

    Output:
      --outdir              Output directory (default: ./results)

    Example:
      nextflow run DupCaller.nf \\
        --sample_name sample1 \\
        --read1 sample1_R1.fq.gz \\
        --read2 sample1_R2.fq.gz \\
        --normal_bam normal.bam \\
        --reference hg38.fa \\
        --germline_vcf gnomad.vcf.gz \\
        --noise_mask noise.bed.gz \\
        --threads 8

    """.stripIndent()
}

if (params.help) {
    helpMessage()
    exit 0
}

// Validate required parameters
if (!params.sample_name) {
    error "Error: --sample_name is required"
}

if (!params.reference) {
    error "Error: --reference is required"
}


include { INPUT_CHECK                                                           } from './subworkflows/input_check'
include { TRIM_BARCODES                                                         } from './modules/trim_barcodes'
include { BWA_ALIGN                                                             } from './modules/bwa_align'
include { MARK_DUPLICATES                                                       } from './modules/mark_duplicates'
include { CALL_VARIANTS                                                         } from './modules/call_variants'
include { ESTIMATE_BURDEN                                                       } from './modules/estimate_burden'



/*
 * Main workflow
 */
workflow {
    // Check reference files
    ref_ch = channel.fromPath(params.reference, checkIfExists: true)
    ref_h5 = params.reference+".ref.h5"
    tn_h5 = params.reference+".tn.h5"
    hp_h5 = params.reference+".hp.h5"
    indexes = channel.fromPath(params.reference+".{amb,ann,bwt,pac,sa}", checkIfExists: true).collect()

    // Handle optional files
    germline_ch = params.germline_vcf ? channel.fromPath(params.germline_vcf, checkIfExists: true) : channel.value(file('NO_GERMLINE'))
    germline_index_ch = params.germline_vcf ? channel.fromPath(params.germline_vcf + ".tbi", checkIfExists: true) : channel.value(file('NO_GERMLINE_IDX'))
    
    noise_ch = params.noise_mask ? channel.fromPath(params.noise_mask, checkIfExists: true) : channel.value(file('NO_NOISE'))
    noise_index_ch = params.noise_mask ? channel.fromPath(params.noise_mask + ".tbi", checkIfExists: true) : channel.value(file('NO_NOISE_IDX'))
    
    indel_ch = params.indel_bed ? channel.fromPath(params.indel_bed, checkIfExists: true) : channel.value(file('NO_INDEL'))
    gene_ch = params.gene_bed ? channel.fromPath(params.gene_bed, checkIfExists: true) : channel.value(file('NO_GENE')) 

    // Tumor sample processing
    if (!params.skip_trim && params.read1 && params.read2) {
        // Start from FASTQ
        reads_ch = INPUT_CHECK(channel.value(file(params.samplesheet)))
        trimmed = TRIM_BARCODES(reads_ch)

        // Combine with reference
        align_input = trimmed.combine(ref_ch).combine(indexes)
        aligned = BWA_ALIGN(align_input)
    } else if (params.skip_trim) {
        error "BAM file input not yet implemented. Please provide --read1 and --read2"
    }

    // Mark duplicates
    markdup = MARK_DUPLICATES(aligned)

    // Normal sample handling
    if (params.normal_bam) {
        normal_ch = channel.fromPath(params.normal_bam, checkIfExists: true)
    } else if (params.normal_read1 && params.normal_read2 && !params.skip_normal) {
        normal_reads = channel.of(["normal_${params.sample_name}", file(params.normal_read1), file(params.normal_read2)])
        normal_trimmed = TRIM_BARCODES(normal_reads)
        normal_align_input = normal_trimmed.combine(ref_ch).combine(indexes)
        normal_aligned = BWA_ALIGN(normal_align_input)
        normal_markdup = MARK_DUPLICATES(normal_aligned)
        normal_ch = normal_markdup.map { it[1] }
    } else {
        normal_ch = channel.value(file('NO_NORMAL'))
    }

    // Variant calling
    // Add h5 files to the call_input tuple
    ref_files_ch = channel.of([file(ref_h5), file(tn_h5), file(hp_h5)])
    call_input = markdup.combine(ref_ch).combine(ref_files_ch)

    variants = CALL_VARIANTS(
        call_input,
        normal_ch,
        germline_ch.first(),
        germline_index_ch.first(),
        noise_ch.first(),
        noise_index_ch.first(),
        indel_ch
    )

    // Burden estimation
    burden_input = variants.combine(ref_ch).combine(ref_files_ch)
    ESTIMATE_BURDEN(burden_input, gene_ch)
}

workflow.onComplete {
    log.info """
    ================================================================
    Pipeline completed at: ${workflow.complete}
    Execution status: ${workflow.success ? 'SUCCESS' : 'FAILED'}
    Duration: ${workflow.duration}
    Output directory: ${params.outdir}
    ================================================================
    """.stripIndent()
}
