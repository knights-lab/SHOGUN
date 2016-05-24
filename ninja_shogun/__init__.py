import os

from ninja_shogun.utilities.logger import Logger
from ninja_shogun.utilities.settings import initialize_settings

SETTINGS = initialize_settings()
# Opens logger to write to log and/or stdout
# First stores original console location as a variable for error handling
shogun_log_path = os.path.join(SETTINGS.data_dir, "shogun_log.txt")
LOGGER = Logger(shogun_log_path, use_std_out=False)
