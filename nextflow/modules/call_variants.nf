process CALL_VARIANTS {
    tag "${sample_name}"
    publishDir "${params.outdir}/${sample_name}/variants", mode: 'copy'
    cpus params.threads

    label 'dupcaller_nextflow'

    input:
    tuple val(sample_name), path(bam), path(bai), path(metrics), path(reference), path(ref_h5), path(tn_h5), path(hp_h5)
    path(normal_bam)
    path(germline_vcf)
    path(germline_index)
    path(noise_mask)
    path(noise_index)
    path(indel_bed)

    output:
    tuple val(sample_name),
          path("${sample_name}_snv.vcf"),
          path("${sample_name}_indel.vcf"),
          path("${sample_name}_coverage.bed.gz"),
          path("${sample_name}_coverage.bed.gz.tbi"),
          path("${sample_name}_trinuc_by_duplex_group.txt"),
          path("${sample_name}_duplex_group_stats.txt"),
          path("${sample_name}_stats.txt"),
          path("${sample_name}.amp.tn.txt"),
          path("${sample_name}.amp.id.txt"),
          path("${sample_name}.dmg.tn.txt"),
          path("${sample_name}.dmg.id.txt")

    script:
    def normal_arg = normal_bam.name != 'NO_NORMAL' ? "-n ${normal_bam}" : ""
    def germline_arg = germline_vcf.name != 'NO_GERMLINE' ? "-g ${germline_vcf}" : ""
    def noise_arg = noise_mask.name != 'NO_NOISE' ? "-m ${noise_mask}" : ""
    def indel_arg = indel_bed.name != 'NO_INDEL' ? "-id ${indel_bed}" : ""

    """
    DupCaller.py call \\
        -b ${bam} \\
        -f ${reference} \\
        -o ${sample_name} \\
        -p ${task.cpus} \\
        -r ${params.regions} \\
        ${normal_arg} \\
        ${germline_arg} \\
        ${noise_arg} \\
        ${indel_arg} \\
        -maf ${params.max_af} \\
        -gaf ${params.germline_af_cutoff} \\
        -d ${params.min_n_depth} \\
        -z ${params.max_zero_qual_frac} \\
        -tt ${params.trim_template} \\
        -tr ${params.trim_read} \\
        -mq ${params.mapq} \\
        -w ${params.window_size}
    """
}
