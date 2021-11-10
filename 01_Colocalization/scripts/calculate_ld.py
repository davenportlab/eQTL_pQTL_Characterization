import subprocess


def calculate_ld(vcf_dir, samples_file, snps_file, chromosome):

    if chromosome == 'X':
        vcf_file = f'{vcf_dir}/CCDG_14151_B01_GRM_WGS_2020-08-05_chr{chromosome}.filtered.eagle2-phased.vcf.gz'
    else:
        vcf_file = f'{vcf_dir}/CCDG_14151_B01_GRM_WGS_2020-08-05_chr{chromosome}.filtered.shapeit2-duohmm-phased.vcf.gz'

    bcf_proc = subprocess.Popen([
        'bcftools', 'view', '-R', snps_file, '-S', samples_file, '--force-samples', vcf_file
    ], stdout=subprocess.PIPE)

    vcf_proc = subprocess.Popen([
        'vcftools', '--vcf', '-', '--geno-r2', '--stdout'
    ], stdin=bcf_proc.stdout, stdout=subprocess.PIPE)
    
    vcf_output, _ = vcf_proc.communicate()

    output = vcf_output.decode('UTF-8').strip().split('\n')
    del output[0]  # Remove the header from the output

    ld_pairs = [ld_pair.strip().split('\t') for ld_pair in output]

    # Some SNPs are multi-allelic and represented with separate recores in the VCF file
    # Remove these by identifying duplicates
    snp_pairs = dict()
    for i, ld_pair in enumerate(ld_pairs):
        key = f'{ld_pair[1]}-{ld_pair[2]}'
        if key not in snp_pairs:
            snp_pairs[key] = list()
        snp_pairs[key].append(i)

    idx_to_del = list()
    for key, value in snp_pairs.items():
        if len(value) > 1:
            idx_to_del += value
        
    for i in sorted(idx_to_del, reverse=True):
        del ld_pairs[i]

    return ld_pairs
