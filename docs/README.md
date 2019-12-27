SHOGUN installation tutorial
=======

### Development Installation
1. Do this in a terminal:
```
conda create -n shogun -c knights-lab shogun
source activate shogun
```

2. Remove SHOGUN and install via the github master branch. This will keep all the conda dependencies installed.
```
conda uninstall shogun
pip install git+https://github.com/knights-lab/SHOGUN.git --no-cache-dir --upgrade
```

Optional: You can reinstall to the newest git version of SHOGUN at anytime via the command:
```
pip install git+https://github.com/knights-lab/SHOGUN.git --no-cache-dir --upgrade
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

If you are on MSI the database for SHOGUN is located at:
```
/home/knightsd/hillm096/globus/SHOGUN/rep82
```

Else you can access the rep82 database on AWS by running :

```
wget <path_to_folder>/shogun_db_links.txt
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
5. Can you get bash tab completion working from ```bin/shogun-complete.sh```?
6. Do any of the tests fail, and can you upload any issues to the issues page?
