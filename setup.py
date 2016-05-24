from setuptools import setup
from glob import glob
import os

__author__ = "Knights Lab"
__copyright__ = "Copyright (c) 2016--, %s" % __author__
__credits__ = ["Benjamin Hillmann", "Dan Knights", "Gabe Al-Ghalith", "Tonya Ward", "Pajua Vangay"]
__email__ = "hillmannben@gmail.com"
__license__ = "GPL"
__maintainer__ = "Benjamin Hillmann"
__version__ = "0.0.1-dev"

long_description = ''

setup(
    name='ninja_shogun',
    version=__version__,
    packages=['ninja_shogun'],
    url='',
    license=__license__,
    author=__author__,
    author_email=__email__,
    description='',
    long_description=long_description,
    scripts=glob(os.path.join('scripts', '*py')),
    keywords='',
    install_requires=['click', 'scipy', 'numpy', 'pandas']
)
