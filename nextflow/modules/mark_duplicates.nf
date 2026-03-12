process MARK_DUPLICATES {
    tag "${sample_name}"
    publishDir "${params.outdir}/${sample_name}/markdup", mode: 'copy'

    label 'dupcaller_nextflow'

    input:
    tuple val(sample_name), path(bam), path(bai)

    output:
    tuple val(sample_name), path("${sample_name}.mkdped.bam"), path("${sample_name}.mkdped.bam.bai"), path("${sample_name}.mkdp_metrics.txt")

    script:
    """
    gatk MarkDuplicates \\
        -I ${bam} \\
        -O ${sample_name}.mkdped.bam \\
        -M ${sample_name}.mkdp_metrics.txt \\
        --READ_NAME_REGEX "(?:.*:)?([0-9]+)[^:]*:([0-9]+)[^:]*:([0-9]+)[^:]*\$" \\
        --DUPLEX_UMI \\
        --TAGGING_POLICY OpticalOnly \\
        --BARCODE_TAG DB

    samtools index ${sample_name}.mkdped.bam
    """
}
