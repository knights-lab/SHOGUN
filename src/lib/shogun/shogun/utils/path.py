#!/usr/bin/env python
import os


def verify_make_dir(outdir):
    if not os.path.exists(outdir):
        os.makedirs(outdir)
