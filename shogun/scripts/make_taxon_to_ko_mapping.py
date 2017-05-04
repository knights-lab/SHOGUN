# usage 
# *.py GMG.microbe.genes.tax GMG.microbe.to.KO.map.txt .8 strain2ko_output.txt species2ko_output.txt
#
# file1: (GMG.microbe.genes.tax)
# NW_002477085       k__Eukaryota; p__; c__; o__; f__Hexamitidae; g__Giardia; s__Giardia_intestinalis; t__Giardia_lamblia_ATCC_50803
#
# file 2: (GMG.microbe.to.KO.map.txt)
# NC_019701.1     K00003  K00005  K00012  K00013  K00014  K00020  K00024  K00027  K00030  K00033  K00036  K00052  K00053  K00058
#
# parameter 3 is "confidence" - only genes in this fraction of strains within a species
# will be called for that species

import sys
import re
from collections import defaultdict, Counter
from math import ceil

rs2strain = {} # {rsID:strain_name}
strain2ko = defaultdict(list) # {strain_name:[KO1, KO2, ...]}
species2strain = defaultdict(list) # {species:[strain1, strain2, ...]}
species2ko = defaultdict(list) # {species_name:[KO1, KO2, ...]}

# read in rs2strain mapping
print('Building rs2strain mapping...')
count = 0
for line in open(sys.argv[1],'r'):
    count += 1
    if count % 100000 == 0:
        sys.stdout.write(str(count) + ' ')
        sys.stdout.flush()
    words = line.strip().split('\t')
    rsID = words[0].split("|")[-1]
    strain_name = words[1].replace(' ','')
    rs2strain[rsID] = strain_name

sys.stdout.write('\n')

# read in strain2ko mapping
# also build species2strain mapping
print('Calculating strain ko content...')
sys.stdout.write('The following strains are not in the rs2strain mapping:')
sys.stdout.flush()
for line in open(sys.argv[2],'r'):
    words = line.strip().split('\t')
    rsID = '.'.join(words[0].split('.')[:-1])
    if rsID in rs2strain:
        strain_name = rs2strain[rsID]
        strain2ko[strain_name] = words[1:]
        species = ';'.join(strain_name.split(';')[:7])
        species2strain[species].append(strain_name)
    else:
        sys.stdout.write(' ' + rsID)
        sys.stdout.flush()
sys.stdout.write('\n')

# output strain KO table
print('Writing strain table...')
strain_names = sorted(strain2ko.keys())
# output new table
with open(sys.argv[4],'w') as f:
    for strain_name in strain_names:
        f.write(strain_name + '\t' + '\t'.join(sorted(strain2ko[strain_name])) + '\n')


# calculate species table - only include genes covering 80% of strains
print('Calculating species ko content...')
species_names = sorted(species2strain.keys())
for species_name in species_names:
    strains = species2strain[species_name]
    n = len(strains)
    
    ko_counter = Counter()
    for strain in strains:
        kos = set(strain2ko[strain])
        for ko in kos:
            ko_counter[ko] += 1

    # get list of kos present in at least confidence*n strains
    threshold = ceil(float(sys.argv[3]) * n)
    for ko in ko_counter:
        if ko_counter[ko] >= threshold:
            species2ko[species_name].append(ko)

# output new table
print('Writing species table...')
with open(sys.argv[5],'w') as f:
    for species_name in species_names:
        f.write(species_name + '\t' + '\t'.join(sorted(species2ko[species_name])) + '\n')
    
