for f in split/*.fna; do echo split/GMG_shearaw.fna; mkdir GMG_shearaw; mv split/GMG_shearaw.fna GMG_shearaw; done
mkdir split; cd split; split -l 22000000 ../GMG_shear.fna GMG_shear
python merge_taxon_counts.py GMG_sheara*/taxon_counts.tsv output.txt
