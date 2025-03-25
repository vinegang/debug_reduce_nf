


workflow Tumor_RNAseq_WF {

take:
    tumor_rnaseq_samplesheet


main:
// Parse the samplesheet to generate fastq tuples
//samples = Channel.fromPath("Tumor_RNAseq.csv")
samples = tumor_rnaseq_samplesheet
.splitCsv(header:true)
.filter { row -> row.type == "tumor_DNA" || row.type == "cell_line_DNA" || row.type == "tumor_RNA" || row.type == "cell_line_RNA" || row.type == "normal_DNA" }
.map { row ->
    def meta = [:]
    meta.id    =  row.sample
    meta.lib   =  row.library
    meta.sc    =  row.sample_captures
    meta.casename  = row.casename
    meta.type     = row.type
    meta.diagnosis =row.Diagnosis
    def fastq_meta = []
    fastq_meta = [ meta,  file(row.read1), file(row.read2)  ]

    return fastq_meta
}

samples_branch = samples.branch{
        exome: it[0].type == "tumor_DNA" || it[0].type == "cell_line_DNA" || it[0].type == "normal_DNA"
        rnaseq: it[0].type == "tumor_RNA" || it[0].type == "cell_line_RNA"
}
samples_branch.exome.view()


}