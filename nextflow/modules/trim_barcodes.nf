process TRIM_BARCODES {
    tag "${sample_name}"
    publishDir "${params.outdir}/${sample_name}/trimmed", mode: 'copy'

    label 'dupcaller_nextflow'

    input:
    tuple val(sample_name), path(read1), path(read2), val(barcode_pattern)

    output:
    tuple val(sample_name), path("${sample_name}_1.fastq"), path("${sample_name}_2.fastq")

    script:
    """
    DupCaller.py trim \\
        -i ${read1} \\
        -i2 ${read2} \\
        -p ${barcode_pattern} \\
        -o ${sample_name}
    """
}
