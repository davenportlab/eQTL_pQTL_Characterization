nextflow.enable.dsl=2

params.genome = "/lustre/scratch118/humgen/resources/rna_seq_genomes/Homo_sapiens.GRCh38.dna.primary_assembly.fa"
params.genome_annotation = "/lustre/scratch118/humgen/resources/rna_seq_genomes/Homo_sapiens.GRCh38.99.gtf"
params.star_genome_dir = "~/gains_team282/epigenetics/star_genome_index/"
params.bowtie2_genome_dir = "~/gains_team282/epigenetics/bowtie2_genome_index/"
params.star = true
params.bowtie2 = true


process IndexGenomeSTAR {

    label "index_star"

    """
    mkdir -p $params.genome_dir

    STAR \
        --runThreadN 16 \
        --runMode genomeGenerate \
        --genomeDir $params.star_genome_dir \
        --genomeFastaFiles $params.genome \
        --sjdbGTFfile $params.genome_annotation \
        --sjdbOverhang 99
    """
}

process IndexGenomeBowtie2 {

    label "index_bowtie2"

    """
    mkdir -p $params.bowtie2_genome_dir

    bowtie2-build $params.genome $params.bowtie2_genome_dir/Homo_sapiens.GRCh38.99
    """
}

workflow {

    if (params.star)
        IndexGenomeSTAR()

    if (params.bowtie2)
        IndexGenomeBowtie2()
}