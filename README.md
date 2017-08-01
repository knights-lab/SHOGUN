Shallow shotgun sequencing
=======
Shallow seq pipeline for optimal shotgun data usage

## Installation
These installation instructions are streamlined for Linux and macOS systems. The tool SHOGUN is installable on Windows with a few minor tweaks to this tutorial. This package requires anaconda, which is a system agnostic package and virtual environment manager. Follow the installation instructions for your system at <http://conda.pydata.org/miniconda.html>.

Once anaconda is installed, get the environment file:

```
wget https://raw.githubusercontent.com/knights-lab/SHOGUN/master/environment.yml
```

Then install the requirements into the environment 'shogun':
```
conda env create -f environment.yml
```

Next, you need the tools UTree and BURST on your path.

#### BURST
TODO: Pending installation instructions...

#### UTree
TODO: Pending installation instructions...

### Development Installation

Once anaconda is installed, create an environment:
```
conda create -n shogun python=3
```

Now activate the environment.

```
# OSX, Linux
source activate shogun
```

With the shogun environment activated, install the developmental SHOGUN toolchain.

```
# If you want to use bowtie2
conda install -c bioconda bowtie2

# SHOGUN
pip install git+https://github.com/knights-lab/SHOGUN.git --no-cache-dir --upgrade
```

With the flags provided to pip, copying and pasting any of these commands will redo the installation if a failure happened.
