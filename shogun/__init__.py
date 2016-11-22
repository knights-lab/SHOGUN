from ninja_utils.config import Settings
from ninja_utils import Logger
from .config.settings import shogun_settings

SETTINGS = Settings('shogun', shogun_settings)
LOGGER = Logger(logfp=SETTINGS.settings['log'], log_persist=SETTINGS.settings['log_persists'])

__all__ = ['config',
           'parsers',
           'scripts',
           'wrappers']
