#!/usr/bin/env nextflow

/*
========================================================================================
    rnaseq-editing-qc
========================================================================================
    IGF2BP3 Knockout RNA-seq Analysis Pipeline
    Based on nf-core/rnaseq framework
    Author: Bioinformatics Research Lab
    GitHub: https://github.com/your-org/rnaseq-editing-qc
----------------------------------------------------------------------------------------
*/

nextflow.enable.dsl = 2

// Print pipeline header
log.info """
========================================================================================
    R N A S E Q - E D I T I N G - Q C   P I P E L I N E
========================================================================================
    Input samplesheet : ${params.input}
    Genome            : ${params.genome}
    Outdir            : ${params.outdir}
    KO gene           : ${params.ko_gene}
    KO condition      : ${params.ko_condition}
    Control condition : ${params.control_condition}
========================================================================================
""".stripIndent()

// Validate mandatory parameters
if (!params.input) { exit 1, "ERROR: '--input' samplesheet not specified." }
if (!params.genome && !params.fasta) { exit 1, "ERROR: '--genome' or '--fasta' must be specified." }

// Import workflows
include { RNASEQ_PIPELINE } from './workflows/rnaseq'

/*
========================================================================================
    MAIN WORKFLOW
========================================================================================
*/

workflow {
    RNASEQ_PIPELINE ()
}

/*
========================================================================================
    ON COMPLETE
========================================================================================
*/

workflow.onComplete {
    log.info ( workflow.success ? """
========================================================================================
    Pipeline completed successfully!
    Results: ${params.outdir}
    Duration: ${workflow.duration}
========================================================================================
""" : """
========================================================================================
    Pipeline FAILED
    Exit status: ${workflow.exitStatus}
    Error message: ${workflow.errorMessage}
========================================================================================
""" )
}
