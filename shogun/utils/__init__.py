"""
Copyright 2015-2017 Knights Lab, Regents of the University of Minnesota.

This software is released under the GNU Affero General Public License (AGPL) v3.0 License.
"""

from .last_common_ancestor import build_lca_map

import subprocess
import os
import hashlib
from collections import defaultdict


def run_command(cmd, shell=False):
    """
    Run prepared behave command in shell and return its output.
    :param cmd: Well-formed behave command to run.
    :param shell: Force subprocess to use shell, not recommended
    :return:
    """

    try:
        cmd = [str(i) for i in cmd]
        proc = subprocess.Popen(
            cmd,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            shell=shell,
            universal_newlines=True,
            cwd=os.getcwd(),
        )

        out, err = proc.communicate()

        if proc.returncode != 0:
            raise AssertionError("exit code is non zero: %d\n%s\%s" % (proc.returncode, out, err))

        return proc.returncode, out, err
    except subprocess.CalledProcessError as e:
        raise AssertionError("Called Process Error: %s" % e)


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


__all__ = ['build_lca_map', 'run_command', 'hash_file', 'read_checksums']
