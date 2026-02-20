/*
========================================================================================
    SUBWORKFLOW: INPUT_CHECK_FASTQ
    Validates FASTQ samplesheet and emits [ meta, [ fastq_1, fastq_2 ] ] tuples
========================================================================================
*/

nextflow.enable.dsl = 2

workflow INPUT_CHECK_FASTQ {

    take:
    samplesheet

    main:

    Channel
        .fromPath(samplesheet)
        .splitCsv(header: true, sep: ',')
        .map { row -> parse_fastq_row(row) }
        .set { ch_reads }

    emit:
    reads    = ch_reads
    versions = Channel.empty()
}

def parse_fastq_row(LinkedHashMap row) {
    def required = ['sample', 'fastq_1', 'fastq_2', 'condition', 'replicate']
    required.each { col ->
        if (!row.containsKey(col) || !row[col]) {
            exit 1, """
ERROR: FASTQ samplesheet missing required column '${col}' or value is empty.
       Expected columns: sample, fastq_1, fastq_2, condition, replicate
       Check assets/samplesheet_template.csv for the correct format.
"""
        }
    }

    def meta = [
        id        : row.sample,
        condition : row.condition,
        replicate : row.replicate.toInteger(),
        single_end: false
    ]

    def fq1 = file(row.fastq_1)
    def fq2 = file(row.fastq_2)

    if (!fq1.exists()) exit 1, "ERROR: FASTQ R1 does not exist: ${row.fastq_1}  (sample: ${row.sample})"
    if (!fq2.exists()) exit 1, "ERROR: FASTQ R2 does not exist: ${row.fastq_2}  (sample: ${row.sample})"

    return [ meta, [ fq1, fq2 ] ]
}
