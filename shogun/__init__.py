"""
Copyright 2015-2017 Knights Lab, Regents of the University of Minnesota.

This software is released under the GNU Affero General Public License (AGPL) v3.0 License.
"""
from .config.settings import shogun_settings

# SETTINGS = Settings('shogun', shogun_settings)
# LOGGER = Logger(logfp=SETTINGS.settings['log'], log_persist=SETTINGS.settings['log_persists'])

__all__ = ['config',
           'parsers',
           'scripts',
           'wrappers']
