params.study = 'calderon_et_al'
params.assay = 'rna_seq'
params.sra_table = '../metadata/reads_calderon_et_al_rna_seq.txt'

Channel
    .fromPath(params.sra_table)
    .splitCsv(header:true)
    .map{row -> row.Run}
    .set{accessions}

process dumpFastq {

    input:
    val accession from accessions

    """
    mkdir -p ~/gains_team282/epigenetics/$params.study/$params.assay/
    for i in 1 2 3 4 5; do fastq-dump --outdir ~/gains_team282/epigenetics/$params.study/$params.assay/ $accession && break || sleep 15; done
    chmod 444 ~/gains_team282/epigenetics/$params.study/$params.assay/${accession}.fastq
    """
}
