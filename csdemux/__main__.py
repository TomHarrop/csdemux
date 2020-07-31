#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from pkg_resources import resource_filename
import argparse
import logging
import psutil
import snakemake


#############
# FUNCTIONS #
#############

def parse_arguments():
    parser = argparse.ArgumentParser(
        prog='csdemux')
    parser.add_argument(
        '-n',
        help='Dry run',
        dest='dry_run',
        action='store_true')
    default_threads = 4
    parser.add_argument(
        '--threads',
        help=('Number of threads. Default: %i' % default_threads),
        metavar='int',
        type=int,
        dest='threads',
        default=default_threads)
    default_mem_gb = int(psutil.virtual_memory().available * .8 // 1e9)
    parser.add_argument(
        '--mem_gb',
        help=('Number of threads. Default: %i' % default_threads),
        metavar='int',
        type=int,
        dest='mem_gb',
        default=default_mem_gb)
    parser.add_argument(
        '--restart_times',
        required=False,
        help='number of times to restart failing jobs (default 0)',
        type=int,
        dest='restart_times',
        default=0)
    parser.add_argument(
        '--r1',
        required=True,
        help='Muxed R1 file',
        type=str,
        dest='r1_file')
    parser.add_argument(
        '--r2',
        required=True,
        help='Muxed R2 file',
        type=str,
        dest='r2_file')
    parser.add_argument(
        '--samples_csv',
        required=True,
        help='Sample csv (see README)',
        type=str,
        dest='samples_csv')
    parser.add_argument(
        '--outdir',
        required=True,
        help='Output directory',
        type=str,
        dest='outdir')

    args = vars(parser.parse_args())
    return(args)


########
# MAIN #
########

def main():
    # set up log
    logging.basicConfig(
        format='%(asctime)s %(levelname)-8s %(message)s',
        datefmt='%Y-%m-%d %H:%M:%S',
        level=logging.INFO)

    # get the snakefile
    snakefile = resource_filename(__name__, 'Snakefile')
    logging.debug(f'Using snakefile {snakefile}')

    # get args
    args = parse_arguments()
    logging.debug(f'Entrypoint args\n{args}')

    # run the pipeline
    snakemake.snakemake(
        snakefile=snakefile,
        config=args,
        cores=args['threads'],
        lock=False,
        printreason=True,
        printshellcmds=True,
        dryrun=True if args['dry_run'] else False,
        restart_times=args['restart_times'])


if __name__ == '__main__':
    main()
