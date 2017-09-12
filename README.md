[![DOI](https://zenodo.org/badge/51028464.svg)](https://zenodo.org/badge/latestdoi/51028464)

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

### Testing your install

For testing, we are currently using the built in python unittests. In order to run the test suite, change directory into the root folder of the repository. Then run:

```
python -m unittest discover shogun
```

# Documentation

#### SHOGUN help for Command-Line
SHOGUN is a command line application. It is meant to be run with a single command. The helpful for the command is below.

```
Usage: shogun [OPTIONS] COMMAND [ARGS]...

  SHOGUN command-line interface

  --------------------------------------

Options:
  --log [debug|info|warning|critical]
                                  The log level to record.
  --shell / --no-shell            Use the shell for Python subcommands (not
                                  recommended).
  --version                       Show the version and exit.
  -h, --help                      Show this message and exit.

Commands:
  align            Run a SHOGUN alignment algorithm.
  assign_taxonomy  Run the SHOGUN taxonomic profile algorithm on...
  coverage         Show confidence of coverage of microbes.
  functional       Run the SHOGUN functional algorithm on a...
  normalize        Normalize a taxonomic profile by median...
  pipeline         Run the SHOGUN pipeline, including taxonomic...
  redistribute     Run the SHOGUN redistribution algorithm on a...
```

#### align
  The command ```align``` runs the respective taxonomic aligner on a linearized, demultiplexed FASTA using either burst, bowtie2, or utree.

```
Usage: shogun align [OPTIONS]

  Run a SHOGUN alignment algorithm.

Options:
  -a, --aligner [all|bowtie2|burst|utree]
                                  The aligner to use.  [default: burst]
  -i, --input PATH                The file containing the combined seqs.
                                  [required]
  -d, --database PATH             The path to the database folder.
  -o, --output PATH               The output folder directory  [default: /mnt/
                                  c/Users/bhill/code/SHOGUN/results-170828]
  -t, --threads INTEGER           Number of threads to use.
  -h, --help                      Show this message and exit.
```

#### assign_taxonomy

```
Usage: shogun assign_taxonomy [OPTIONS]

  Run the SHOGUN taxonomic profile algorithm on an alignment output.

Options:
  -a, --aligner [bowtie2|burst|burst-tax|utree]
                                  The aligner to use.  [default: burst]
  -i, --input PATH                The output alignment file.
                                  [required]
  -d, --database PATH             The path to the database folder.
  -o, --output PATH               The coverage table.  [default: /mnt/c/Users/
                                  bhill/code/SHOGUN/taxatable-170828.txt]
  -h, --help                      Show this message and exit.
```


#### coverage

```
Usage: shogun coverage [OPTIONS]

  Show confidence of coverage of microbes.

Options:
  -i, --input PATH                The output BURST alignment.
                                  [required]
  -d, --database PATH             The path to the folder containing the
                                   database.  [required]
  -o, --output PATH               The coverage table.  [default: /mnt/c/Users/
                                  bhill/code/SHOGUN/coverage-170828.txt]
  -l, --level [genus|species|strain]
                                  The level to collapse to.
  -h, --help                      Show this message and exit.
```

#### functional

This command assigns function at a certain taxonomic level. Lower level KEGG IDs are assigned to higher level KEGG IDs through plurality voting. Note that plasmids are not included the KEGG ID annotation.

```
Usage: shogun functional [OPTIONS]

  Run the SHOGUN functional algorithm on a taxonomic profile.

Options:
  -i, --input PATH                The taxatable.  [required]
  -d, --database PATH             The path to the folder containing the
                                  function database.  [required]
  -o, --output PATH               The output file  [default: /mnt/c/Users/bhil
                                  l/code/SHOGUN/results-170828]
  -l, --level [genus|species|strain]
                                  The level to collapse to.
  -h, --help                      Show this message and exit.
```

### normalize

```
Usage: shogun normalize [OPTIONS]

  Normalize a taxonomic profile by median depth.

Options:
  -i, --input PATH   The output taxatable.  [required]
  -o, --output PATH  The taxatable output normalized by median depth.
                     [default: /mnt/c/Users/bhill/code/SHOGUN/taxatable.normal
                     ized-170828.txt]
  -h, --help         Show this message and exit.
```

#### pipeline

```
Usage: shogun pipeline [OPTIONS]

  Run the SHOGUN pipeline, including taxonomic and functional profiling.

Options:
  -a, --aligner [all|bowtie2|burst|utree]
                                  The aligner to use [Note: default burst is
                                  capitalist, use burst-tax if you want to
                                  redistribute].  [default: burst]
  -i, --input PATH                The file containing the combined seqs.
                                  [required]
  -d, --database PATH             The path to the database folder.
  -o, --output PATH               The output folder directory  [default: /mnt/
                                  c/Users/bhill/code/SHOGUN/results-170828]
  -l, --level [kingdom|phylum|class|order|family|genus|species|strain|all|off]
                                  The level to collapse taxatables and
                                  functions to (not required, can specify
                                  off).
  --function / --no-function      Run functional algorithms. **This will
                                  normalize the taxatable by median depth.
  --capitalist / --no-capitalist  Run capitalist with burst post-align or not.
  -t, --threads INTEGER           Number of threads to use.
  -h, --help                      Show this message and exit.
```


#### redistribute
  This command redistributes the reads at a certain taxonomic level. This assumes that you have a BIOM txt file output from SHOGUN align, or even a summarized table from redistribute at a lower level.

```
Usage: shogun redistribute [OPTIONS]

  Run the SHOGUN redistribution algorithm on a taxonomic profile.

Options:
  -i, --input PATH                The taxatable.  [required]
  -d, --database PATH             The path to the database folder.  [required]
  -l, --level [kingdom|phylum|class|order|family|genus|species|strain|all]
                                  The level to collapse to.
  -o, --output PATH               The output file  [default: /mnt/c/Users/bhil
                                  l/code/SHOGUN/taxatable-170828.txt]
  -h, --help                      Show this message and exit.
```
