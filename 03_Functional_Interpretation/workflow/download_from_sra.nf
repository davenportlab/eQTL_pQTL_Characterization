params.study = "calderon_et_al"
params.assay = "rna_seq"
params.library_layout = "paired"
params.sra_table = "../metadata/reads_calderon_et_al_rna_seq.txt"

Channel
    .fromPath(params.sra_table)
    .splitCsv(header:true)
    .map{row -> row.Run}
    .set{accessions}

process dumpFastq {

    label "dump_fastq"

    input:
        val accession from accessions

    script:

        """
        mkdir -p ~/gains_team282/epigenetics/$params.study/raw/$params.assay/
        """

        if (params.library_layout == "single")

            """
            for i in {1..5}
            do
                fastq-dump \
                    --outdir ~/gains_team282/epigenetics/$params.study/raw/$params.assay/ \
                    --gzip \
                    --read-filter pass \
                    --dumpbase \
                    --clip \
                    $accession && break || sleep 15
            done

            mv ~/gains_team282/epigenetics/$params.study/raw/$params.assay/${accession}_pass.fastq.gz ~/gains_team282/epigenetics/$params.study/raw/$params.assay/${accession}.fastq.gz

            chmod 444 ~/gains_team282/epigenetics/$params.study/raw/$params.assay/${accession}.fastq.gz
            """
        
        else if (params.library_layout == "paired")

            """
            for i in {1..5}
            do
                fastq-dump \
                    --outdir ~/gains_team282/epigenetics/$params.study/raw/$params.assay/ \
                    --gzip \
                    --skip-technical \
                    --readids \
                    --read-filter pass \
                    --dumpbase \
                    --split-3 \
                    --clip \
                    $accession && break || sleep 15
            done

            mv ~/gains_team282/epigenetics/$params.study/raw/$params.assay/${accession}_pass_1.fastq.gz ~/gains_team282/epigenetics/$params.study/raw/$params.assay/${accession}_1.fastq.gz
            mv ~/gains_team282/epigenetics/$params.study/raw/$params.assay/${accession}_pass_2.fastq.gz ~/gains_team282/epigenetics/$params.study/raw/$params.assay/${accession}_2.fastq.gz

            chmod 444 ~/gains_team282/epigenetics/$params.study/raw/$params.assay/${accession}_1.fastq.gz
            chmod 444 ~/gains_team282/epigenetics/$params.study/raw/$params.assay/${accession}_2.fastq.gz
            """
}
