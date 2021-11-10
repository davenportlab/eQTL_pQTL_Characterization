#!/usr/bin/env python

import click
import pandas as pd

@click.group()
def cli():
    pass

#----------------------------------------------------------
# FastQC
#----------------------------------------------------------

@cli.command()
@click.argument('input_paths', nargs=-1)
@click.argument('output_path', nargs=1)
@click.argument('base_names', nargs=1)
def fastqc(input_paths, output_path, base_names):

    base_names = base_names.split(",")

    fastqc_data = dict()

    for fastqc_log, base_name in zip(input_paths, base_names):

        with open(fastqc_log, 'r') as f_in:

            for line in f_in:

                if '>>' in line and '>>END_MODULE' not in line:

                    components = line[2:].strip().split('\t')
                    key = components[0].strip()
                    value = components[1].strip()

                    if key not in fastqc_data:
                        fastqc_data[key] = dict()

                    fastqc_data[key][base_name] = value

    fastqc_df = pd.DataFrame(fastqc_data)
    fastqc_df.to_csv(output_path)

#----------------------------------------------------------
# Trim Galore
#----------------------------------------------------------

@cli.command()
@click.argument('input_paths', nargs=-1)
@click.argument('output_path', nargs=1)
@click.argument('base_names', nargs=1)
def trim_galore(input_paths, output_path, base_names):

    base_names = base_names.split(",")

    trim_galore_data = dict()
    trim_galore_data['Total Reads Processed'] = dict()
    trim_galore_data['Reads with Adapters'] = dict()
    trim_galore_data['Reads Written'] = dict()
    trim_galore_data['Total Basepairs Processed'] = dict()
    trim_galore_data['Quality-Trimmed Basepairs'] = dict()
    trim_galore_data['Total Basepairs Written'] = dict()
    trim_galore_data['Short Sequences Removed'] = dict()

    for trim_galore_log, base_name in zip(input_paths, base_names):

        with open(trim_galore_log, 'r') as f_in:

            for line in f_in:

                if 'Total reads processed' in line:
                    components = line.split()
                    value = components[len(components) - 1]
                    trim_galore_data['Total Reads Processed'][base_name] = value
                
                if 'Reads with adapters' in line:
                    components = line.split()
                    value = ' '.join(components[len(components) - 2:len(components)])
                    trim_galore_data['Reads with Adapters'][base_name] = value

                if 'Reads written (passing filters)' in line:
                    components = line.split()
                    value = ' '.join(components[len(components) - 2:len(components)])
                    trim_galore_data['Reads Written'][base_name] = value

                if 'Total basepairs processed' in line:
                    components = line.split()
                    value = ' '.join(components[len(components) - 2:len(components)])
                    trim_galore_data['Total Basepairs Processed'][base_name] = value

                if 'Quality-trimmed' in line:
                    components = line.split()
                    value = ' '.join(components[len(components) - 3:len(components)])
                    trim_galore_data['Quality-Trimmed Basepairs'][base_name] = value

                if 'Total written (filtered)' in line:
                    components = line.split()
                    value = ' '.join(components[len(components) - 3:len(components)])
                    trim_galore_data['Total Basepairs Written'][base_name] = value

                if 'Sequences removed because they became shorter than the length cutoff' in line:
                    components = line.split()
                    value = ' '.join(components[len(components) - 2:len(components)])
                    trim_galore_data['Short Sequences Removed'][base_name] = value

    trim_galore_df = pd.DataFrame(trim_galore_data)
    trim_galore_df.to_csv(output_path)

#----------------------------------------------------------
# STAR Alignment
#----------------------------------------------------------

@cli.command()
@click.argument('input_paths', nargs=-1)
@click.argument('output_path', nargs=1)
@click.argument('base_names', nargs=1)
def star_alignment(input_paths, output_path, base_names):

    base_names = base_names.split(",")

    align_final_data = dict()

    for align_final_log, base_name in zip(input_paths, base_names):

        with open(align_final_log, 'r') as f_in:

            for line in f_in:

                if '|' in line:

                    components = line.split('|')
                    key = components[0].strip()
                    value = components[1].strip()

                    if key not in align_final_data:
                        align_final_data[key] = dict()
                    
                    align_final_data[key][base_name] = value

    align_final_df = pd.DataFrame(align_final_data)
    align_final_df.to_csv(output_path)

#----------------------------------------------------------
# Mark Duplicates
#----------------------------------------------------------

@cli.command()
@click.argument('input_paths', nargs=-1)
@click.argument('output_path', nargs=1)
@click.argument('base_names', nargs=1)
def mark_duplicates(input_paths, output_path, base_names):

    base_names = base_names.split(",")

    mark_duplicates_data = dict()

    for mark_duplicates_log, base_name in zip(input_paths, base_names):

        components = list()

        with open(mark_duplicates_log, 'r') as f_in:

            for line in f_in:

                if 'LIBRARY' in line:

                    components = line.strip().split()
                    for component in components:
                        if component not in mark_duplicates_data:
                            mark_duplicates_data[component] = dict()
                
                if 'Unknown Library' in line:
                    
                    values = line.strip().split()
                    for component, value in zip(components, values):
                        mark_duplicates_data[component][base_name] = value

    mark_duplicates_df = pd.DataFrame(mark_duplicates_data)
    mark_duplicates_df.to_csv(output_path)

#----------------------------------------------------------
# BAM Statistics
#----------------------------------------------------------

@cli.command()
@click.argument('input_paths', nargs=-1)
@click.argument('output_path', nargs=1)
@click.argument('base_names', nargs=1)
def bam_statistics(input_paths, output_path, base_names):

    base_names = base_names.split(",")

    bam_statistics_data = dict()

    for bam_statistic, base_name in zip(input_paths, base_names):

        with open(bam_statistic, 'r') as f_in:

            for line in f_in:

                if ':' in line:

                    components = line.split(':')
                    key = components[0].strip()
                    value = components[1].strip()

                    if key not in bam_statistics_data:
                        bam_statistics_data[key] = dict()
                    
                    bam_statistics_data[key][base_name] = value

    bam_statistics_df = pd.DataFrame(bam_statistics_data)
    bam_statistics_df.to_csv(output_path)

#----------------------------------------------------------
# BAM GC
#----------------------------------------------------------

@cli.command()
@click.argument('input_paths', nargs=-1)
@click.argument('output_path', nargs=1)
@click.argument('base_names', nargs=1)
def bam_gc(input_paths, output_path, base_names):

    base_names = base_names.split(",")

    bam_gc_data = dict()
    bam_gc_data['GC Mean'] = dict()
    bam_gc_data['GC SD'] = dict()

    for bam_gc, base_name in zip(input_paths, base_names):

        bam_gc_df = pd.read_csv(bam_gc, sep='\t')

        total_reads = bam_gc_df['read_count'].sum()

        bam_gc_data['GC Mean'][base_name] = (bam_gc_df['GC%'] * bam_gc_df['read_count']).sum() / total_reads
        bam_gc_data['GC SD'][base_name] = ((bam_gc_df['GC%'] - bam_gc_data['GC Mean'][base_name]) ** 2).sum() / (total_reads - 1)

    bam_gc_df = pd.DataFrame(bam_gc_data)
    bam_gc_df.to_csv(output_path)

#----------------------------------------------------------
# Feature Counts
#----------------------------------------------------------

@cli.command()
@click.argument('input_path', nargs=1)
@click.argument('output_path', nargs=1)
@click.argument('base_names', nargs=1)
def feature_counts(input_path, output_path, base_names):

    base_names = base_names.split(",")

    feature_counts_df = pd.read_csv(input_path, sep='\t')
    colnames = feature_counts_df['Status']
    feature_counts_df = feature_counts_df.transpose()
    feature_counts_df = feature_counts_df.reindex(base_names)
    feature_counts_df.columns = colnames
    feature_counts_df.to_csv(output_path)

#----------------------------------------------------------
# Main
#----------------------------------------------------------

if __name__ == '__main__':

    cli()
