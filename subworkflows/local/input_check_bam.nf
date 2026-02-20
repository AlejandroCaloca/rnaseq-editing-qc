/*
========================================================================================
    SUBWORKFLOW: INPUT_CHECK_BAM
    Validates BAM samplesheet and emits [ meta, bam, bai ] tuples
========================================================================================
*/

nextflow.enable.dsl = 2

workflow INPUT_CHECK_BAM {

    take:
    samplesheet  // file: path to CSV

    main:

    Channel
        .fromPath(samplesheet)
        .splitCsv(header: true, sep: ',')
        .map { row -> parse_bam_row(row) }
        .set { ch_bam }

    emit:
    bam      = ch_bam       // channel: [ val(meta), path(bam), path(bai) ]
    versions = Channel.empty()
}

def parse_bam_row(LinkedHashMap row) {

    // ── Required columns ─────────────────────────────────────────────────────
    def required = ['sample', 'bam', 'condition', 'replicate']
    required.each { col ->
        if (!row.containsKey(col) || !row[col]) {
            exit 1, """
ERROR: BAM samplesheet missing required column '${col}' or value is empty.
       Expected columns: sample, bam, bai (optional), condition, replicate
       Check assets/samplesheet_bam_template.csv for the correct format.
"""
        }
    }

    // ── Build meta map ───────────────────────────────────────────────────────
    def meta = [
        id        : row.sample,
        condition : row.condition,
        replicate : row.replicate.toInteger(),
        single_end: false
    ]

    // ── Resolve BAM file ─────────────────────────────────────────────────────
    def bam_file = file(row.bam)
    if (!bam_file.exists()) {
        exit 1, "ERROR: BAM file does not exist:\n       ${row.bam}\n       Sample: ${row.sample}"
    }
    if (!bam_file.name.endsWith('.bam')) {
        exit 1, "ERROR: File does not appear to be a BAM file:\n       ${row.bam}\n       Sample: ${row.sample}"
    }

    // ── Resolve BAI index ────────────────────────────────────────────────────
    // If 'bai' column is present and filled, use it.
    // Otherwise look for <bam>.bai or <bam_noext>.bai next to the BAM.
    def bai_file
    if (row.containsKey('bai') && row.bai) {
        bai_file = file(row.bai)
        if (!bai_file.exists()) {
            exit 1, "ERROR: BAI index file does not exist:\n       ${row.bai}\n       Sample: ${row.sample}"
        }
    } else {
        // Try common naming conventions
        def bai_option1 = file("${row.bam}.bai")          // file.bam.bai
        def bai_option2 = file("${row.bam}".replaceAll(/\.bam$/, '.bai'))  // file.bai
        if (bai_option1.exists()) {
            bai_file = bai_option1
        } else if (bai_option2.exists()) {
            bai_file = bai_option2
        } else {
            log.warn "WARNING: No BAI index found for ${row.sample}. " +
                     "Add a 'bai' column to your samplesheet, or index your BAM first:\n" +
                     "         samtools index ${row.bam}"
            bai_file = file("${row.bam}.bai")  // Placeholder — will fail downstream if missing
        }
    }

    return [ meta, bam_file, bai_file ]
}
