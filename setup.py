from setuptools import setup, find_packages
from glob import glob
import os

__author__ = "Knights Lab"
__copyright__ = "Copyright (c) 2016--, %s" % __author__
__credits__ = ["Benjamin Hillmann", "Dan Knights", "Gabe Al-Ghalith", "Tonya Ward", "Pajua Vangay"]
__email__ = "hillmannben@gmail.com"
__license__ = "MIT"
__maintainer__ = "Benjamin Hillmann"
__version__ = "0.0.1-dev"

long_description = ''

setup(
    name='ninja_shogun',
    version=__version__,
    packages=find_packages(),
    url='',
    license=__license__,
    author=__author__,
    author_email=__email__,
    description='',
    long_description=long_description,
    # scripts=glob(os.path.join('ninja_shogun', 'scripts', '*py')),
    keywords='',
    install_requires=['click', 'scipy', 'numpy', 'pandas', 'cytoolz', 'pyfaidx'],
    entry_points={
        'console_scripts': [
            'shogun_bt2_align = ninja_shogun.scripts.shogun_bt2_align:shogun_bt2_align',
            'shogun_bt2_db = ninja_shogun.scripts.shogun_bt2_db:shogun_bt2_db',
            'shogun_bt2_lca = ninja_shogun.scripts.shogun_bt2_lca:shogun_bt2_lca',
            'shogun_bt2_capitalist = ninja_shogun.scripts.shogun_bt2_capitalist:shogun_bt2_capitalist',
            'shogun_utree_db = ninja_shogun.scripts.shogun_utree_db:shogun_utree_db',
            'shogun_utree_capitalist = ninja_shogun.scripts.shogun_utree_capitalist:shogun_utree_capitalist',
        ]
    },
)
