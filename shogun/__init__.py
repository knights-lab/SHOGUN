"""
Copyright 2015-2017 Knights Lab, Regents of the University of Minnesota.

This software is released under the GNU Affero General Public License (AGPL) v3.0 License.
"""

from ._version import get_versions
__version__ = get_versions()['version']
del get_versions

from .config.settings import shogun_settings

import logging

def _logging_setup():
    # Set up the logger
    log_formatter = logging.Formatter('%(message)s', datefmt='%m/%d/%Y %I:%M:%S %p')
    root_logger = logging.getLogger()

    console_handler = logging.StreamHandler()
    console_handler.setFormatter(log_formatter)
    root_logger.addHandler(console_handler)
    return root_logger

logger = _logging_setup()

__all__ = [
    'aligners',
    'config',
    'function',
    'parsers',
    'scripts',
    'taxonomy',
    'wrappers']
