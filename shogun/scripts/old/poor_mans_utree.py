#!/usr/bin/env python
import sys
import subprocess
import os
import multiprocessing

fasta = sys.argv[1]
fasta = os.path.abspath(fasta)
tax = ".".join(fasta.split(".")[:-1] + ["tax"])
nthreads = str(multiprocessing.cpu_count())

def make_mapping(path, level):
    with open(path) as inf:
        outfp = '.'.join(path.split('.')[:-1] + [str(level)] + [path.split('.')[-1]])
        with open(outfp, 'w') as outf:
            for line in inf:
                tax = line.split()[-1]
                tax = ';'.join(tax.split(';')[:level])
                outf.write('%s\t%s\n' % (line.split()[0], tax))

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

for i in range(8):
    make_mapping(tax, i+1)
    tax_level = ".".join(fasta.split(".")[:-1] + [str(i+1)] + ["tax"])
    utr = ".".join(fasta.split(".")[:-1] + [str(i+1)])
    ctr = ".".join(fasta.split(".")[:-1] + [str(i+1)] + ["ctr"])
    output_time = ".".join(fasta.split(".")[:-1] + [str(i+1)] + ["txt"])

    run_command(["/usr/bin/time", "-v", "--output=" + output_time] +
                ['utree-build', fasta, tax_level, utr, nthreads])
    run_command(["/usr/bin/time", "-v", "--append", "--output=" + output_time] +
                ['utree-compress', utr, ctr])
    run_command(["/usr/bin/time", "-v", "--append", "--output=" + output_time] +
                ['rm', "-f", utr])
