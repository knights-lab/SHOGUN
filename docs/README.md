SHOGUN installation tutorial
=======

### Upgrading

To upgrade, simply run the following command with your environment active.
```
TMPDIR=/dev/shm pip install git+https://github.com/knights-lab/SHOGUN.git --no-cache-dir --upgrade
```

### Installation
These installation instructions are streamlined for Linux and macOS systems. The tool SHOGUN is installable on Windows with a few minor tweaks to this tutorial. This package requires anaconda, which is a system agnostic package and virtual environment manager. Follow the installation instructions for your system at <http://conda.pydata.org/miniconda.html>.

Once anaconda is installed, get the environment file:

```
wget https://raw.githubusercontent.com/knights-lab/SHOGUN/master/environment.yml
```

Then install the requirements into the environment 'shogun':
```
conda env create --file environment.yml
```

Next, you need the tools UTree, Bowtie2 and BURST on your path. These are contained in the lab SOP.

```
export PATH="/export/scratch/ben/shogun_bin$PATH"
```

Now activate the environment.

```
# OSX, Linux
source activate shogun
```

With the shogun environment activated, we can now run the following to check the SHOGUN version number. This is linked to the Github commit string.

```
shogun --version
shogun, version v0.0.1+293.g6531389
```

To verify your installation, make sure to run the test suite. For testing, we are currently using the built in python unittests. In order to run the test suite, change directory into the root folder of the repository. Then run:

```
cd /project/flatiron2/ben/projects/SHOGUN/SHOGUN
python -m unittest discover shogun
```

After the testing is complete, we can try running our analysis. To do so, run the following command:

```
shogun --log debug pipeline --input ./shogun/tests/data/combined_seqs.fna --database /project/flatiron2/analysis_SHOGUN/data/references/rep82 --output ~/scratch_shogun --level strain
```

If everything up to the testing is working, below is the location of an example shallow shotgun file.
```
/home/grad00/hillm096/combined_seqs.fna
```

Can you run SHOGUN to answer the following questions?
1. What is the most abundant genus, species, and strain taxonomy?
2. What is the expected coverage of each taxonomy at the strain level?
3. What is the most abundant KEGG ID and the coverage of each module?
4. What is the variance between the profiles of each taxonomic profilers utree, burst, and bowtie2?
