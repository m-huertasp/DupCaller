process ESTIMATE_BURDEN {
    tag "${sample_name}"
    publishDir "${params.outdir}/${sample_name}/burden", mode: 'copy'

    label 'dupcaller_nextflow'

    input:
    tuple val(sample_name),
          path(snv_vcf),
          path(indel_vcf),
          path(coverage),
          path(coverage_index),
          path(trinuc),
          path(duplex_stats),
          path(stats),
          path(amp_tn),
          path(amp_id),
          path(dmg_tn),
          path(dmg_id),
          path(reference),
          path(ref_h5),
          path(tn_h5),
          path(hp_h5)
    path(gene_bed)

    output:
    tuple val(sample_name),
          path("**${sample_name}_sbs_burden.txt"),
          path("**${sample_name}_indel_burden.txt"),
          path("**${sample_name}_sbs_96_corrected.txt"),
          path("**${sample_name}_sbs_96_corrected.png"),
          path("**${sample_name}_sbs_burden_by_min_read_group_size.txt"),
          path("**${sample_name}_sbs_burden_by_min_read_group_size.png"),
          path("**${sample_name}_duplex_allele_counts.txt")

    script:
    def gene_arg = gene_bed.name != 'NO_GENE' ? "-gb ${gene_bed}" : ""
    def clonal_arg = params.estimate_clonal ? "-c true" : ""
    def dilute_arg = params.estimate_dilute ? "-d true" : ""

    """
    # Copy variant call outputs to current directory for estimate command

    mkdir -p ${sample_name}
    ln -s "\$(readlink -f ${snv_vcf})" ${sample_name}/.
    ln -s "\$(readlink -f ${indel_vcf})" ${sample_name}/.
    ln -s "\$(readlink -f ${coverage})" ${sample_name}/.
    ln -s "\$(readlink -f ${coverage_index})" ${sample_name}/.
    ln -s "\$(readlink -f ${trinuc})" ${sample_name}/.
    ln -s "\$(readlink -f ${duplex_stats})" ${sample_name}/.
    ln -s "\$(readlink -f ${stats})" ${sample_name}/.
    ln -s "\$(readlink -f ${amp_tn})" ${sample_name}/.
    ln -s "\$(readlink -f ${amp_id})" ${sample_name}/.
    ln -s "\$(readlink -f ${dmg_tn})" ${sample_name}/.
    ln -s "\$(readlink -f ${dmg_id})" ${sample_name}/.

    DupCaller.py estimate \\
        -i ${sample_name} \\
        -f ${reference} \\
        -r ${params.regions} \\
        ${gene_arg} \\
        ${clonal_arg} \\
        ${dilute_arg}
    """
}
