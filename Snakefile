#!/usr/bin/env python3

from pathlib import Path

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


rule target:
    input:
        unpack(demux_target)

checkpoint demultiplex:
    input:
        r1 = 'data/muxed/Undetermined_S0_L002_R1_001.fastq.gz',
        r2 = 'data/muxed/Undetermined_S0_L002_R2_001.fastq.gz'
    output:
        directory('output/000_tmp/reads'),
        stats = 'output/010_demux/stats.txt'
    log:
        'output/logs/demultiplex.log'
    params:
        out = 'output/000_tmp/reads/%_r1.fastq.gz',
        out2 = 'output/000_tmp/reads/%_r2.fastq.gz'
        streams = workflow.cores // 2
    threads:
        workflow.cores
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
