/*
 * pipeline input parameters
 */

params.reads = "$projectDir/data/ggal/gut_{1,2}.fastq"
params.output_dir = "$PWD/results"
params.database = "/KRAKEN_DB/"
nextflow.enable.dsl=2


log.info """\
    KRAKEN  - TAPIR   P I P E L I N E
    ============================================
    output_dir       : ${params.output_dir}
    """
    .stripIndent()

/*
 * define the `index` process that creates a binary index
 * given the transcriptome file
 */
 
process KRAKEN {
    publishDir "${params.output_dir}", mode:'copy'
    tag "speciation of ${meta}"
    
    conda "bioconda::kraken=1.1.1"
    
    input:
    tuple val(meta), path(reads)
    val database

    output:
    path('*.out'),       emit: kraken_out_ch
    path('*report.txt'), emit: kraken_report_ch
    path "versions.yml", emit: versions_ch

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta}"
    """
    # speciation step using kraken
    kraken --db ${database} --threads $task.cpus --fastq-input --gzip-compressed ${prefix}.fastq.gz > ${prefix}.kraken.out
    
    # generate report from kraken output of microbial reads
    kraken-report --db ${database} ${prefix}.kraken.out > ${prefix}.kraken.report.txt
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        kraken: \$(echo \$(kraken --version 2>&1) | sed 's/^.*Kraken version //; s/ .*\$//')
    END_VERSIONS
    """
}




workflow  {
         reads_ch = channel
                          .fromPath( params.reads, checkIfExists: true )
			  .map { file -> tuple(file.simpleName, file) }		
	 
       
	 KRAKEN(reads_ch, params.database)	     
}

