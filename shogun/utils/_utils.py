"""
Copyright 2015-2017 Knights Lab, Regents of the University of Minnesota.

This software is released under the GNU Affero General Public License (AGPL) v3.0 License.
"""

import hashlib
import os
import subprocess
from collections import defaultdict
from contextlib import contextmanager
from timeit import default_timer

import numpy as np
import scipy.sparse as ss

from shogun import logger


@contextmanager
def elapsed_timer():
    start = default_timer()
    elapser = lambda: default_timer() - start
    yield lambda: elapser()
    end = default_timer()
    elapser = lambda: end-start


def run_command(cmd, shell=False, stdout=subprocess.PIPE, stderr=subprocess.STDOUT):
    """
    Run prepared behave command in shell and return its output.
    :param stderr:
    :param stdout:
    :param cmd: Well-formed behave command to run.
    :param shell: Force subprocess to use shell, not recommended
    :return:
    """

    try:
        cmd = [str(i) for i in cmd]

        if not stdout:
            stdout = open(os.devnull, 'w')
        if not stderr:
            stderr = open(os.devnull, 'w')

        logger.debug(" ".join(cmd))
        with elapsed_timer() as elapsed:
            with subprocess.Popen(
                " ".join(cmd) if shell else cmd,
                stdout=stdout,
                stderr=stderr,
                shell=shell,
                universal_newlines=True,
                bufsize=1,
                cwd=os.getcwd(),
            ) as proc:
                log_subprocess_output(proc.stdout)
        logger.debug("%.2f seconds" % elapsed())
        logger.debug("Subprocess finished.")

        #if proc.returncode != 0:
            #raise AssertionError("exit code is non zero: %d\n%s" % (proc.returncode, " ".join(cmd)))
        return proc.returncode, "", ""
    except subprocess.CalledProcessError as e:
        raise AssertionError("Called Process Error: %s" % e)


def log_subprocess_output(pipe):
    for line in pipe:
        line = line.rstrip()
        if line:
            if not line.startswith('Search Progress'):
                logger.debug(line)


def hash_file(filename):
    h = hashlib.sha1()

    with open(filename, "rb") as file:
        chunk = 0
        while chunk != b'':
            chunk = file.read(1024)
            h.update(chunk)

    return h.hexdigest()


def read_checksums(filename):
    with open(filename) as inf:
        return defaultdict(str, dict([line.split() for line in inf]))


def save_csr_matrix(filename, matrix, row_names, column_names):
    """Save compressed sparse row (csr) matrix to file.

    Based on http://stackoverflow.com/a/8980156/232571

    """
    assert filename.endswith('.npz')
    attributes = {
        'data': matrix.data,
        'indices': matrix.indices,
        'indptr': matrix.indptr,
        'shape': matrix.shape,
        'rownames': row_names,
        'columnnames': column_names,
    }
    np.savez(filename, **attributes)


def load_csr_matrix(filename):
    """Load compressed sparse row (csr) matrix from file.

    Based on http://stackoverflow.com/a/8980156/232571

    """
    assert filename.endswith('.npz')
    loader = np.load(filename)
    args = (loader['data'], loader['indices'], loader['indptr'])
    matrix = ss.csr_matrix(args, shape=loader['shape'])
    return loader['rownames'], loader['columnnames'], matrix


def read_fasta(fh):
    """
    :return: tuples of (title, seq)
    """
    title = None
    data = None
    for line in fh:
        if line[0] == ">":
            if title:
                yield (title.strip(), data)
            title = line[1:]
            data = ''
        else:
            data += line.strip()
    if not title:
        yield None
    yield (title.strip(), data)
