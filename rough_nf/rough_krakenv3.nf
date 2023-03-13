/*
 * pipeline input parameters
 */

params.reads = "$projectDir/data/ggal/gut_{1,2}.fastq"
params.output_dir = "$PWD/results"
params.database = "/KRAKEN_DB/"
params.taxon = 'S'
params.sequencing_date = 'GYYMMDD'
nextflow.enable.dsl=2


log.info """\
    KRAKEN  - TAPIR   P I P E L I N E
    ============================================
    output_dir       : ${params.output_dir}
    database         : ${params.database}
    taxon            : ${params.taxon}
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
    path('*.out'),                        emit: kraken_out_ch
    tuple val(meta), path('*report.txt'), emit: kraken_report_ch
    path "versions.yml",                  emit: versions_ch

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


process EXTRACT_TAXON_SPECIFIC_INFO {
    //publishDir "${params.output_dir}", mode:'copy'
    tag "extract $taxon from ${sample_id}.kraken.report.txt"


    input:
    tuple val(sample_id), path(kraken_report)
    val(taxon)


    output:
    path("*.${taxon}.kraken.txt"), emit: taxon_kraken_ch


    script:
    """
    
    echo "percentage\tcladeReads\ttaxonReads\ttaxRank\ttaxID\tspecies" > ${sample_id}.${taxon}.kraken.txt

    grep "\\s${taxon}\\s" ${kraken_report} >> ${sample_id}.${taxon}.kraken.txt

    """
}


process COMBINE_KRAKEN_REPORT_FROM_TAXA {
    publishDir "${params.output_dir}", mode:'copy'
    tag "combine kraken.report.txt for $taxon"


    input:
    path(kraken_taxon_report_files)
    val(taxon)
    val(date)


    output:
    path("combined_kraken_report_${taxon}_${date}.txt"), emit: kraken_comb_ch


    script:
    """
    KRAKEN_TAXON_REPORT_FILES=(${kraken_taxon_report_files})

    for index in \${!KRAKEN_TAXON_REPORT_FILES[@]}; do
    KRAKEN_TAXON_REPORT_FILE=\${KRAKEN_TAXON_REPORT_FILES[\$index]}
    sample_id=\${KRAKEN_TAXON_REPORT_FILE%.${taxon}.kraken.txt}

    # add header line if first file
    if [[ \$index -eq 0 ]]; then
      echo "samplename\t\$(head -1 \${KRAKEN_TAXON_REPORT_FILE})" >> combined_kraken_report_${taxon}_${date}.txt
    fi
    #echo "\${sample_id}\t\$(awk 'FNR>=2 {print}' \${KRAKEN_TAXON_REPORT_FILE})" >> combined_kraken_report_${taxon}_${date}.txt
    #awk '{print \${sample_id}, \$0' \${KRAKEN_TAXON_REPORT_FILE} >> combined_kraken_report_${taxon}_${date}.txt
    awk -F '\\t' 'FNR>=2 { print FILENAME, \$0 }' \${KRAKEN_TAXON_REPORT_FILE} |  sed 's/\\.${taxon}\\.kraken\\.txt//g' >> combined_kraken_report_${taxon}_${date}.txt
    done


    """
}



workflow  {
         reads_ch = channel
                          .fromPath( params.reads, checkIfExists: true )
			  .map { file -> tuple(file.simpleName, file) }		
	 
       
	 KRAKEN(reads_ch, params.database)
	 
	 EXTRACT_TAXON_SPECIFIC_INFO(KRAKEN.out.kraken_report_ch, params.taxon)	     
	 
	 collected_kraken_taxon_ch = EXTRACT_TAXON_SPECIFIC_INFO.out.taxon_kraken_ch.collect( sort: {a, b -> a[0].getBaseName() <=> b[0].getBaseName()} )

         COMBINE_KRAKEN_REPORT_FROM_TAXA(collected_kraken_taxon_ch, params.taxon, params.sequencing_date)

}

