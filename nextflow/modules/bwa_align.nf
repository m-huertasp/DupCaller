process BWA_ALIGN {
    tag "${sample_name}"
    publishDir "${params.outdir}/${sample_name}/aligned", mode: 'copy'
    cpus params.threads

    label 'dupcaller_nextflow'

    input:
    tuple val(sample_name), path(read1), path(read2), path(reference), path(ref_amb), path(ref_ann), path(ref_bwt), path(ref_pac), path(ref_sa)

    output:
    tuple val(sample_name), path("${sample_name}.bam"), path("${sample_name}.bam.bai")

    script:
    // def rg_id = params.read_group_id ?: sample_name
    def rg_id =  "${sample_name}"

    """
    bwa mem -C \\
        -t ${task.cpus} \\
        -R "@RG\\tID:${rg_id}\\tSM:${sample_name}\\tPL:${params.platform}" \\
        ${reference} \\
        ${read1} ${read2} | \\
    samtools sort -@ ${task.cpus-1} -o ${sample_name}.bam

    samtools index -@ ${task.cpus-1} ${sample_name}.bam
    """
}
