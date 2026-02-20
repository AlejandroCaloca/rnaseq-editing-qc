/*
========================================================================================
    SUBWORKFLOW: INPUT_CHECK
    Validates samplesheet and creates channels for downstream processing
========================================================================================
*/

nextflow.enable.dsl = 2

workflow INPUT_CHECK {

    take:
    samplesheet // file: /path/to/samplesheet.csv

    main:

    Channel
        .fromPath(samplesheet)
        .splitCsv(header: true, sep: ',')
        .map { row -> validate_and_parse_row(row) }
        .set { ch_reads }

    emit:
    reads    = ch_reads      // channel: [ val(meta), [ reads ] ]
    versions = Channel.empty()
}

// Function to validate and parse each samplesheet row
def validate_and_parse_row(LinkedHashMap row) {
    // Check required columns
    def required_cols = ['sample', 'fastq_1', 'fastq_2', 'condition', 'replicate']
    required_cols.each { col ->
        if (!row.containsKey(col) || !row[col]) {
            exit 1, "ERROR: Samplesheet missing required column '${col}' or value is empty.\n" +
                     "       Check your samplesheet format against the template in assets/samplesheet_template.csv"
        }
    }

    // Build meta map
    def meta = [
        id         : "${row.sample}",
        condition  : row.condition,
        replicate  : row.replicate.toInteger(),
        single_end : false  // Pipeline requires paired-end
    ]

    // Check FASTQ files exist
    def fastq_1 = file(row.fastq_1)
    def fastq_2 = file(row.fastq_2)

    if (!fastq_1.exists()) {
        exit 1, "ERROR: FASTQ file does not exist: ${row.fastq_1}"
    }
    if (!fastq_2.exists()) {
        exit 1, "ERROR: FASTQ file does not exist: ${row.fastq_2}"
    }

    return [ meta, [ fastq_1, fastq_2 ] ]
}
