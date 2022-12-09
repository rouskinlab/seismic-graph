import click
import pandas as pd
from main import run
from dreem.util.cli_args import *

@click.command()

@optgroup.group('I/O')
@fasta 
@input_dir
@out_dir
@library

@optgroup.group('Selection')
@coords
@primers
@fill

@optgroup.group('Clustering')
@n_clusters
@max_clusters
@signal_thresh
@info_thresh
@include_g_u
@include_del
@min_reads
@convergence_cutoff
@num_runs

@optgroup.group('Miscellaneous')
@verbose


def cli(**args):
    run(**args)

if __name__ == '__main__':
    run()