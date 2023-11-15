## Workflow to predict taxonomy from uncorrected long reads using kraken.
### Usage

```

===============================================
 KRAKEN TAPIR Pipeline version 1.0dev
===============================================
 Mandatory arguments:
         --reads                        Query fastqz file of sequences you wish to supply as input (e.g., "/MIGE/01_DATA/01_FASTQ/T055-8-*.fastq.gz")
         --database                     KRAKEN database directory (full path required, e.g., "/KRAKEN_DB")
         --output_dir                   Output directory to place final combined kraken output (e.g., "/MIGE/01_DATA/06_SPECIES_CHECK")
         --sequencing_date              Sequencing Date (for TAPIR, must start with G e.g., "G230223")

        Optional arguments:
         --taxon                        value must be one of: "S", "G", "F", "O", "C", "P", "K", "D", (Default: "S")
         --help                         This usage statement.
         --version                      Version statement



```


## Introduction
This pipeline calls taxonomy from uncorrected long reads derived from pure cultures or metagenomic samples. This Nextflow pipeline was adapted from the Kraken [github page](https://github.com/DerrickWood/kraken), and the NF Core Modules [github page](https://github.com/nf-core/modules/tree/master/modules/nf-core). Inputs are fastqs specified using `--reads`. 


## Sample command
An example of a command to run this pipeline is:

```
nextflow run main.nf --reads "Sample_files/*.fastq.gz" --output_dir "test2" --database "/KRAKEN_DB" --sequencing_date "G230202" --taxon "S"
```

## Word of Note
This is an ongoing project at the Microbial Genome Analysis Group, Institute for Infection Prevention and Hospital Epidemiology, Üniversitätsklinikum, Freiburg. The project is funded by BMBF, Germany, and is led by [Dr. Sandra Reuter](https://www.uniklinik-freiburg.de/institute-for-infection-prevention-and-control/microbial-genome-analysis.html).


## Authors and acknowledgment
The TAPIR (Track Acquisition of Pathogens In Real-time) team.
