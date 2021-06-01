# csdemux

Preconfigured pipeline for demultiplexing CS's Illumina reads

1. Demultiplex on the barcode field in the fastq header.
2. Analyse the barcode content and generate plots.
3. Filter and trim adaptors, and plot adaptor content.
4. Collate reads into a pair of fastq files for each sample.

## Install

Use the singularity container from the [Releases](https://github.com/TomHarrop/csdemux/releases) tab, or the Docker container from [GHCR](https://github.com/users/TomHarrop/packages/container/package/csdemux).

## Usage

- `threads`: Number of threads to use.
- `mem_gb`: Amount of RAM in GB.
- `restart_times`: Number of times to restart failing jobs.
- `r1`: Muxed R1 file.
- `r2`: Muxed R2 file.
- `samples_csv`: Demux config file. See `data/example.samples.csv` for an example. A csv with the following columns:
    - `sample`: sample name (will be propagated to output files);
    - `barcode`: sample barcode, will be used to demultiplex;
    - `r1_path` (optional): currently not used;
    - `r2_path` (optional): currently not used;
    - `metadata` (optional): currently not used.
- `outdir`: Output directory.

```
csdemux [-h] [-n] [--threads int] [--mem_gb int]
               [--restart_times RESTART_TIMES] --r1 R1_FILE --r2
               R2_FILE --samples_csv SAMPLES_CSV --outdir OUTDIR

optional arguments:
  -h, --help            show this help message and exit
  -n                    Dry run
  --threads int         Number of threads. Default: 4
  --mem_gb int          Amount of RAM in GB. Default: 15
  --restart_times RESTART_TIMES
                        number of times to restart failing jobs
                        (default 0)
  --r1 R1_FILE          Muxed R1 file
  --r2 R2_FILE          Muxed R2 file
  --samples_csv SAMPLES_CSV
                        Sample csv (see README)
  --outdir OUTDIR       Output directory
```

## Graph

![](graph.svg)

