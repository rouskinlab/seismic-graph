import logging
import os

import numpy as np
import pandas as pd

B_MATCH = 0
B_DELET = 1
B_INS_5 = 2
B_INS_3 = 3
B_SUB_A = 4
B_SUB_C = 5
B_SUB_G = 6
B_SUB_T = 7


def get_files(folder: str, ext: str):
    paths = []
    for reference in os.listdir(folder):
        if os.path.isdir(os.path.join(folder, reference)):
            for section in os.listdir(os.path.join(folder, reference)):
                if not section.endswith(ext):
                    continue
                paths.append(os.path.join(folder, reference, section))
    else:
        if reference.endswith(ext):
            paths.append(os.path.join(folder, reference))
    return paths


def fastq_to_df(fastq_file):    
    df, data = pd.DataFrame(), pd.read_csv(fastq_file, sep='\t', header=None)
    df['reference'] = data.iloc[np.arange(0,len(data),4)].reset_index(drop=True)
    df['sequence'] = data.iloc[np.arange(1,len(data),4)].reset_index(drop=True)
    df['quality'] = data.iloc[np.arange(3,len(data),4)].reset_index(drop=True)
    return df


def sam_to_df(path):
    skiprows = 0
    with open(path) as f:
        while f.readline().startswith('@'):
            skiprows += 1
    df = pd.read_csv(path, sep='\t', skiprows=skiprows, header=None)
    df = df[df.columns[:11]]
    df.columns = ['QNAME', 'FLAG', 'RNAME', 'POS', 'MAPQ', 'CIGAR', 'RNEXT', 'PNEXT','TLEN', 'SEQ', 'QUAL']
    return df


def sort_fastq_pairs(fq1s, fq2s):
    if type(fq1s) is str:
        fq1s = [fq1s]
    if type(fq2s) is str:
        fq2s = [fq2s]
    assert len(fq1s) >= len(fq2s), 'More fq2s than fq1s'
    for f2 in fq2s:
        if f2.replace('_R2', '_R1') not in fq1s:
            raise ValueError(f'No matching pair for {f2}')
    fq1s = [f2.replace('_R2', '_R1') for f2 in fq2s] + [f for f in fq1s if f.replace('_R1','_R2') not in fq2s]
    fq2s += [None for f in fq1s if f.replace('_R1','_R2') not in fq2s]
    samples = [f.split('/')[-1].split('.')[0].split('_')[:-1] for f in fq1s]
    samples = ['_'.join(s) for s in samples]
    return fq1s, fq2s, samples


def query_muts(muts: np.ndarray, bits: int, sum_up = True, axis=0, set_type = 'superset'):
    """
    Count the number of times a query mutation occurs in each column
    or one column of a set of mutation mut_vectors.
    The counting operation comprises three steps:
    1. bitwise AND to confirm at least one "1" bit is shared, e.g.
       bits: 11110000 & mut_vectors: 00100000 -> 00100000 (True)
       bits: 11110000 & mut_vectors: 00100010 -> 00100000 (True)
       bits: 11110000 & mut_vectors: 00000000 -> 00000000 (False)
    2. bitwise OR to confirm no "1" bit in mut_vectors is not in bits, e.g.
       bits: 11110000 | mut_vectors: 00100000 -> 11110000 =? 11110000 (True)
       bits: 11110000 | mut_vectors: 00100010 -> 11110010 =? 11110000 (False)
       bits: 11110000 | mut_vectors: 00000000 -> 11110000 =? 11110000 (True)
    3. logical AND to confirm that both tests pass, e.g.
       bits: 11110000, mut_vectors: 00100000 -> True  AND True  (True)
       bits: 11110000, mut_vectors: 00100010 -> True  AND False (False)
       bits: 11110000, mut_vectors: 00000000 -> False AND True  (False)

    Arguments
    mut_vectors: NDArray of a set of mutation mut_vectors (2-dimensional)
          or one column in a set of mutation mut_vectors (1-dimensional).
          Data type must be uint8.
    bits: One-byte int in the range [0, 256) representing the mutation
          to be queried. The bits in the int encode the mutation as
          defined above, e.g.
          - 00000010 (int 2) is a deletion
          - 11010001 (int 209) is either substitution to A, G, or T
                               or a match to C
    
    Returns
    if sum_up: 
        int of the number of times the query mutation occurs in mut_vectors
        count: If mut_vectors is 1-dimensional, int of the number of times the
            query mutation occurs in mut_vectors.
            If mut_vectors is 2-dimensional, NDArray with one int for each
            column in mut_vectors.
    if not sum_up:
        bool NDArray with one bool for each column in mut_vectors.
        True if the query mutation occurs in the column, False otherwise.
    """
    if not muts.dtype == np.uint8:
        raise TypeError('mut_vectors must be of type uint8 and not {}'.format(muts.dtype))
    #print(np.logical_and(mut_vectors & bits, (mut_vectors | bits) == bits).sum(axis=0))
    assert isinstance(bits, int) and 0 <= bits < 256

        
    if set_type == 'subset':
        logic_fun = lambda m, b: np.logical_and(m & b, (m | b) == b)
    if set_type == 'superset':
        logic_fun = lambda m, b: np.array(m & b, dtype=bool)
    
    if sum_up:
        return logic_fun(muts, bits).sum(axis=axis)
    else:
        return logic_fun(muts, bits)


def get_num_parallel(n_tasks: int,
                     max_procs: int,
                     parallel: bool,
                     hybrid: bool = False) -> tuple[int, int]:
    """ Determine how to parallelize the tasks.

    Parameters
    ----------
    n_tasks: int (≥ 1)
        Number of tasks to parallelize
    max_procs: int (≥ 1)
        Maximum number of processes to run at one time
    parallel: bool
        Whether to permit multiple tasks to be run in parallel
    hybrid: bool (default: False)
        Whether to allow both multiple tasks to run in parallel and,
        at the same, each task to run multiple processes in parallel

    Returns
    -------
    int (≥ 1)
        Number of tasks to run in parallel
    int (≥ 1)
        Number of processes to run for each task
    """
    if n_tasks >= 1 and max_procs >= 1:
        # This function only works if there is at least one task to
        # parallelize, at least one process is allowed, and parallel
        # is a valid option.
        if parallel:
            # Multiple tasks may be run in parallel. The number of tasks
            # run in parallel cannot exceed 1) the total number of tasks
            # and 2) the user-specified maximum number of processes.
            n_tasks_parallel = min(n_tasks, max_procs)
        else:
            # Otherwise, only one task at a time can be run.
            n_tasks_parallel = 1
        if n_tasks_parallel == 1 or hybrid:
            # Each individual task can be run by multiple processes in
            # parallel, as long as either 1) multiple tasks are not run
            # simultaneously in parallel (i.e. n_tasks_parallel == 1)
            # or 2) the calling function sets hybrid=True, which lets
            # multiple tasks run in parallel and each run with multiple
            # processes. Only the alignment module can simultaneously
            # run multiple tasks and multiple processes for each task
            # because its two most computation-heavy processes (cutadapt
            # and bowtie2) come with their own parallelization abilities
            # that can work independently of Python's multiprocessing
            # module. However, the other modules (e.g. vectoring) are
            # parallelized using the multiprocessing module, which does
            # not support "nesting" parallelization in multiple layers.
            # Because n_tasks_parallel is either 1 or the smaller of
            # n_tasks and n_procs (both of which are ≥ 1), it must be
            # that 1 ≤ n_tasks_parallel ≤ n_procs, and therefore that
            # 1 ≤ n_procs / n_tasks_parallel ≤ n_procs, so the
            # integer quotient must be a valid number of processes.
            n_procs_per_task = max_procs // n_tasks_parallel
        else:
            # Otherwise, only one process can work on each task.
            n_procs_per_task = 1
    else:
        logging.warning("Defaulting to 1 process due to invalid number of "
                        f"tasks ({n_tasks}) and/or processes ({max_procs}).")
        n_tasks_parallel = 1
        n_procs_per_task = 1
    return n_tasks_parallel, n_procs_per_task