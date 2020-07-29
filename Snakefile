#!/usr/bin/env python3

from pathlib import Path
import pandas

# demuxbyname.sh \
# in=r1.fastq \
# in2=r2.fastq \
# delimiter=: \
# column=10 \
# prefixmode=f \
# zl=9 \
# out=%_r1.fastq.gz \
# out2=%_r2.fastq.gz


bbmap = 'shub://TomHarrop/seq-utils:bbmap_38.76'
bioconductor = 'shub://TomHarrop/r-containers:bioconductor_3.11'


def demux_target(wildcards):
    '''
    glob the demux output directory to see which barcodes worked
    '''
    cdir = checkpoints.demultiplex.get(**wildcards).output[0]
    my_bcs = glob_wildcards(Path(cdir, '{bc}_r1.fastq.gz')).bc
    output_dict = {
        'files': expand('output/000_tmp/reads/{barcode}_r{r}.fastq.gz',
                        barcode=my_bcs,
                        r=['1', '2']),
        'directory': cdir}
    return(output_dict)


def kept_indiv_reads(wildcards):
    '''
    Parse the barcode_csv to find which indivs we want to keep
    '''
    # barcode_csv = 'output/010_demux/barcode_distance.csv'
    barcode_csv = checkpoints.plot_barcode_distance.get(
        **wildcards).output['report']
    read_path = 'output/020_reads/{indiv}_r{r}.fastq.gz'
    found_bc = pandas.read_csv(barcode_csv)
    kept_df = found_bc.loc[(found_bc.dist == 0) & (found_bc.Reads > 0)]
    kept_indivs = sorted(set(kept_df['sample']))
    return(expand(read_path, indiv=kept_indivs, r=['1', '2']))


rule target:
    input:
        # unpack(demux_target),
        kept_indiv_reads,
        'output/010_demux/barcode_content.pdf',
        'output/010_demux/barcode_distance.pdf'


def get_indiv_reads(wildcards, params):
    '''
    get a list of reads to concat
    '''
    # barcode_csv = 'output/010_demux/barcode_distance.csv'
    # my_indiv = 'indiv83_mhyp_lincoln'
    # my_dist = 1
    read_path = 'output/000_tmp/reads/{bc}_r{r}.fastq.gz'
    my_dist = int(params.max_dist)
    my_indiv = wildcards.indiv
    barcode_csv = checkpoints.plot_barcode_distance.get(
        **wildcards).output['report']
    found_bc = pandas.read_csv(barcode_csv)
    my_bcs = sorted(set(
        found_bc.loc[(found_bc['sample'] == my_indiv) & (found_bc.dist <= my_dist)]['found_barcode']))
    return(
        {'r1': expand(read_path, bc=my_bcs, r=['1']),
         'r2': expand(read_path, bc=my_bcs, r=['2'])})


rule mv_reads:
    input:
        unpack(get_indiv_reads)
    output:
        r1 = 'output/020_reads/{indiv}_r1.fastq.gz'
        r2 = 'output/020_reads/{indiv}_r2.fastq.gz'
    params:
        max_dist = 1
    shell:
        'cat {input.r1} > {output.r1} & '
        'cat {input.r2} > {output.r2} & '
        'wait'

checkpoint plot_barcode_distance:
    input:
        stats = 'output/010_demux/stats.txt',
        barcodes = 'data/samples.csv',
        foundbc = 'output/010_demux/hamming_distances.Rds'
    output:
        plot = 'output/010_demux/barcode_distance.pdf',
        report = 'output/010_demux/barcode_distance.csv'
    params:
        read_min_freq = 1e-4,
        max_mismatches = 2
    log:
        'output/logs/plot_barcode_distance.log'
    singularity:
        bioconductor
    script:
        'src/plot_barcode_distance.R'

rule calculate_hamming_distance:
    input:
        stats = 'output/010_demux/stats.txt',
        barcodes = 'data/samples.csv'
    output:
        foundbc = 'output/010_demux/hamming_distances.Rds'
    log:
        'output/logs/calculate_hamming_distance.log'
    threads:
        workflow.cores
    singularity:
        bioconductor
    script:
        'src/calculate_hamming_distance.R'

rule plot_barcode_content:
    input:
        stats = 'output/010_demux/stats.txt',
    output:
        plot = 'output/010_demux/barcode_content.pdf',
    log:
        'output/logs/plot_barcode_content.log'
    singularity:
        bioconductor
    script:
        'src/plot_barcode_content.R'

checkpoint demultiplex:
    input:
        r1 = 'data/muxed2/Undetermined_S0_L006_R1_001.fastq.gz',
        r2 = 'data/muxed2/Undetermined_S0_L006_R2_001.fastq.gz'
    output:
        directory('output/000_tmp/reads'),
        stats = 'output/010_demux/stats.txt'
    log:
        'output/logs/demultiplex.log'
    params:
        out = 'output/000_tmp/reads/%_r1.fastq.gz',
        out2 = 'output/000_tmp/reads/%_r2.fastq.gz',
        streams = min(workflow.cores, 40) // 2
    threads:
        min(workflow.cores, 40)
    singularity:
        bbmap
    shell:
        'demuxbyname.sh '
        'in={input.r1} '
        'in2={input.r2} '
        'out={params.out} '
        'out2={params.out2} '
        'stats={output.stats} '
        'streams={params.streams} '
        'delimiter=: '
        'column=10 '
        'prefixmode=f '
        'zl=9 '
        '-Xmx100g '
        '2> {log}'
