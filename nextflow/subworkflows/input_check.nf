//
// Check input samplesheet and get read channels
//


workflow INPUT_CHECK {
    take:
    samplesheet // file: /path/to/samplesheet.csv

    main:
    samplesheet
        .splitCsv ( header:true, sep:',' )
        .map { it -> create_input_channel(it) }
        .set { sample_inputs }

    emit:
    sample_inputs                       // channel: [ val(sample_name), file(row.vcf), file(row.bam) ]
}

// Function to get list of [ sample_name, fastq_1, fastq_2 ]
def create_input_channel(LinkedHashMap row) {

    // add path(s) of the fastq file(s) to the meta map
    def input_meta = []
    if (!file(row.fastq_1).exists()) {
        exit 1, "ERROR: Please check input samplesheet -> Read 1 FastQ file does not exist!\n${row.fastq_1}"
    }
    if (!file(row.fastq_2).exists()) {
        exit 1, "ERROR: Please check input samplesheet -> Read 2 FastQ file does not exist!\n${row.fastq_2}"
    }
    input_meta = [ row.sample_name, file(row.fastq_1), file(row.fastq_2), row.barcode_pattern ]

    return input_meta
}
