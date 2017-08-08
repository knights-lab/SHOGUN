"""
Copyright 2015-2017 Knights Lab, Regents of the University of Minnesota.

This software is released under the GNU Affero General Public License (AGPL) v3.0 License.
"""

from .last_common_ancestor import build_lca_map
from .normalize import normalize_by_median_depth
from shogun import logger, LoggerWriter

import subprocess
import os
import hashlib
from collections import defaultdict
import numpy as np
import scipy.sparse as ss


def run_command(cmd, shell=False, stdout=subprocess.PIPE, stderr=subprocess.STDOUT):
    """
    Run prepared behave command in shell and return its output.
    :param cmd: Well-formed behave command to run.
    :param shell: Force subprocess to use shell, not recommended
    :return:
    """

    try:
        cmd = [str(i) for i in cmd]
        FNULL = open(os.devnull, 'w')
        if not stdout:
            stdout = FNULL
        if not stderr:
            stderr = FNULL

        logger.debug(' '.join(cmd))

        proc = subprocess.Popen(
            cmd,
            stdout=stdout,
            stderr=stderr,
            shell=shell,
            universal_newlines=True,
            cwd=os.getcwd(),
        )

        with proc.stdout:
            log_subprocess_output(proc.stdout)

        returncode = proc.communicate()

        if returncode != 0:
            raise AssertionError("exit code is non zero: %d\n%s\%s" % (proc.returncode, out, err))
        return proc.returncode, "", ""
    except subprocess.CalledProcessError as e:
        raise AssertionError("Called Process Error: %s" % e)

def log_subprocess_output(pipe):
    for line in iter(pipe.readline, b''):
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



__all__ = ['build_lca_map', 'run_command', 'hash_file', 'read_checksums', 'save_csr_matrix', 'load_csr_matrix',
           'normalize_by_median_depth']
