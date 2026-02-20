/*
========================================================================================
    SUBWORKFLOW: INPUT_CHECK_COUNTS
    Validates counts samplesheet and emits [ meta, counts_file ] tuples
========================================================================================
*/

nextflow.enable.dsl = 2

workflow INPUT_CHECK_COUNTS {

    take:
    samplesheet

    main:

    Channel
        .fromPath(samplesheet)
        .splitCsv(header: true, sep: ',')
        .map { row -> parse_counts_row(row) }
        .set { ch_counts }

    emit:
    counts   = ch_counts       // [ val(meta), path(counts_file) ]
    versions = Channel.empty()
}

def parse_counts_row(LinkedHashMap row) {
    def required = ['sample', 'counts_file', 'condition', 'replicate']
    required.each { col ->
        if (!row.containsKey(col) || !row[col]) {
            exit 1, """
ERROR: Counts samplesheet missing required column '${col}' or value is empty.
       Expected columns: sample, counts_file, condition, replicate
       Check assets/samplesheet_counts_template.csv for the correct format.
"""
        }
    }

    def meta = [
        id        : row.sample,
        condition : row.condition,
        replicate : row.replicate.toInteger()
    ]

    def counts_file = file(row.counts_file)
    if (!counts_file.exists()) {
        exit 1, "ERROR: Counts file does not exist: ${row.counts_file}  (sample: ${row.sample})"
    }

    return [ meta, counts_file ]
}
