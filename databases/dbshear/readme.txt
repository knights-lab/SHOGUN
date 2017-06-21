mkdir split; cd split; split -l 22000000 ../GMG_shear.fna GMG_shear

for f in split/*.fna; do echo split/GMG_shearaw.fna; mkdir GMG_shearaw; mv split/GMG_shearaw.fna GMG_shearaw; done

for f in GMG_sheara*; do echo $f; time python ../../shogun/scripts/shogun_embalmer_lca.py -i $f -r /project/flatiron/gabe/GMG.ann.450_95 -d 0.98 -p 90; done

python merge_taxon_counts.py GMG_sheara*/taxon_counts.tsv output.txt
