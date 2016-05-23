#!/usr/bin/env python
import os
import sys

from shogun.utilities.path import verify_make_dir


class Logger(object):
    """ A convenient logging object
        Prints output to a given log file and/or stdout
    """
    def __init__(self, logfp=None, use_std_out=True):
        # note: if logfp directories don't exist, make them.

        if logfp is not None:
            outdir = os.path.abspath(os.path.dirname(os.path.realpath(logfp)))

            # Checks for output directory. Makes it if necessary.

            verify_make_dir(outdir)

        self.logfp = logfp
        logf = open(logfp, 'w')
        logf.close()
        self.use_std_out = use_std_out

    def log(self, msg):
        if not msg.endswith('\n'):
            msg += '\n'
        if self.logfp is not None:
            logf = open(self.logfp, 'a')
            logf.write(msg)
            logf.close()
        if self.use_std_out:
            sys.stdout.write(msg)
