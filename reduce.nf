
nextflow.enable.dsl=2

include {RNAseq_only} from './workflows/RNAseq_only.nf'
include {Exome_only_WF} from './workflows/Exome_only_WF.nf'
include {Tumor_RNAseq_WF} from './workflows/Tumor_RNAseq_WF.nf'


process PREPARE_SAMPLESHEET {

    input:
    path samplesheet
    val genome_version

    output:
    path("*csv")

    script:
    """

    if [ ${genome_version} == "hg19" ]; then
        python ${workflow.projectDir}/bin/split_samplesheet.py ${samplesheet} .
    elif [ ${genome_version} == "mm39" ]; then
        cp ${samplesheet} mouse_rnaseq.csv
    else
        echo "Error: Unknown genome: ${genome_version}"
        exit 1
    fi
    """
}

// Launch workflow by checking the samplesheet availability
workflow {

    samplesheet_ch = Channel.fromPath(params.samplesheet)
    genome_version_ch = Channel.value(params.genome_version)

    // Run PREPARE_SAMPLESHEET process
    prepared_samplesheets = PREPARE_SAMPLESHEET(samplesheet_ch, genome_version_ch)

    prepared_samplesheets.view()

    case_type = prepared_samplesheets.branch {
        rnaseq: it.name == 'RNAseq.csv'
        exome: it.name == 'Exome.csv'
        TR: it.name == 'Tumor_RNAseq.csv'
    }
    case_type.TR|Tumor_RNAseq_WF
    case_type.exome|Exome_only_WF
    case_type.rnaseq|RNAseq_only

}