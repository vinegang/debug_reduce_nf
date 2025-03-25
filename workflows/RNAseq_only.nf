def metadatareducer(inputChannel) {
    return inputChannel.map { meta, file ->
            [ meta.id, meta.casename, meta.diagnosis, meta.lib, file ]}
        .reduce([[:], [], []]) { result, item ->
            def (patient, casename, diagnosis, lib, filePath) = item

            // Dynamically set the id from the meta data
            result[0] = [id: patient, casename: casename, diagnosis: diagnosis]

            // Append lib names
            result[1].add(lib)

            // Append file paths
            result[2].add(filePath)

            return result
        }
}

workflow RNAseq_only {


take:
    rnaseq_samplesheet

main:

//create a sample channel using meta hashmap
samples_rnaseq = rnaseq_samplesheet
.splitCsv(header:true)
.filter { row -> row.type == "tumor_RNA" || row.type == "cell_line_RNA" || row.type == "xeno_RNA"}
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
samples_rnaseq.view()

    inputChannels = Channel.from(
        [[id: 'Test1', lib: 'Test9_T1D_E', sc: 'clin.ex.v1', casename: 'NFtest', type: 'tumor_DNA', diagnosis: 'Osteosarcoma'], 'assets/Test9.txt'],
        [[id: 'Test1', lib: 'Test5_T1D_E', sc: 'clin.ex.v1', casename: 'NFtest', type: 'tumor_DNA', diagnosis: 'Osteosarcoma'], 'assets/Test5.txt'],
        [[id: 'Test1', lib: 'Test1_R_T', sc: 'access', casename: 'NFtest', type: 'tumor_RNA', diagnosis: 'Osterosarcoma'], 'assets/Test1.txt']
    )
  
//uncomment this block and reduce function gets executed even though input RNAseq.csv is missing.  
/*
Combined_snpeff_vcf2txt_ch = metadatareducer(inputChannels)
Combined_snpeff_vcf2txt_ch.view()
*/
}