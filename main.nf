#!/usr/bin/env nextflow

nextflow.enable.dsl=2

// include definitions
include  { helpMessage; Version } from './modules/messages.nf'

// include processes
include { KRAKEN; EXTRACT_TAXON_SPECIFIC_INFO; COMBINE_KRAKEN_REPORT_FROM_TAXA } from './modules/processes.nf'

log.info """\
    ======================================
    KRAKEN  - TAPIR  P I P E L I N E
    ======================================
    output_dir      : ${params.output_dir}
    """
    .stripIndent()


workflow  {
         reads_ch = channel
                          .fromPath( params.reads, checkIfExists: true )
                          .map { file -> tuple(file.simpleName, file) }


         KRAKEN(reads_ch, params.database)

         EXTRACT_TAXON_SPECIFIC_INFO(KRAKEN.out.kraken_report_ch, params.taxon)

         collected_kraken_taxon_ch = EXTRACT_TAXON_SPECIFIC_INFO.out.taxon_kraken_ch.collect( sort: {a, b -> a[0].getBaseName() <=> b[0].getBaseName()} )

         COMBINE_KRAKEN_REPORT_FROM_TAXA(collected_kraken_taxon_ch, params.taxon, params.sequencing_date)

}
