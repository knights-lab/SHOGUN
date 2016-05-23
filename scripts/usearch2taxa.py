#!/usr/bin/env python
# usage usearch2taxa usearchout.txt (usearchout-r2.txt) taxonmap output.txt
# expects usearch blast6 output
#
# by default, if sequence header ends with /1 or /2, truncates
#
# if ref ID has an underscore, only takes the part before the underscore
#
# taxonmap is formatted like the greengenes taxonomy maps: "Refseqid\tTaxonomy"
# can also handle taxonmap that has species (e.g. Escherichia coli)
import sys

# given two tuples that are taxonomic hierarchies,
# returns the deepest shared hierarchy
# or None if no intersection


def lca(tax1, tax2):
    if tax1 is None or tax2 is None:
        return None
    i = min(len(tax1), len(tax2)) - 1
    while i >= 0:
        if tax1[i] == tax2[i]:
            return tax1[:(i + 1)]
        i -= 1
    return None

# given two lists of tuples that are taxonomic hierarchies,
# returns the set that are shared
# or None if no intersection


def commonTaxa(refIDs1, refIDs2, taxon_map):
    if refIDs1 is None or refIDs2 is None:
        return None

    # loop through readh list
    taxa1 = set([taxon_map[refID] for refID in refIDs1])
    taxa2 = set([taxon_map[refID] for refID in refIDs2])

    common = taxa1.intersection(taxa2)
    if len(common) == 0:
        return None
    return common

if __name__ == "__main__":

    paired = False
    usearch_fp = sys.argv[1]
    if len(sys.argv) == 5:
        paired = True
        usearch2_fp = sys.argv[2]
        taxon_fp = sys.argv[3]
        out_fp = sys.argv[4]
    else:
        taxon_fp = sys.argv[2]
        out_fp = sys.argv[3]

    species_only = False

    ref_ids_missing_from_taxonomy = 0

    # load the taxon map
    print "Loading taxonomy map..."
    taxa = {}
    for line in open(taxon_fp, 'r'):
        words = line.strip().split('\t')
        taxonomy = tuple(words[1].split(';'))
        if len(taxonomy) == 1 and words[1] != 'Archaea' and words[1] != 'Eukaryota':
            species_only = True
            taxonomy = tuple(taxonomy[0].split(' '))
        taxa[words[0]] = taxonomy

    # load all assignments
    refIDs1 = {}  # will hold sets of refIDs for each query
    refIDs2 = {}  # will hold sets of refIDs for each query
    taxon_assignments = {}  # will hold sets of taxon tuples for each query

    print "Loading taxonomy assignments..."
    for line in open(usearch_fp, 'r'):
        words = line.strip().split('\t')
        query_id = '/'.join(words[0].split('/')[:-1])
        ref_id = words[1].split('_')[0]
        if not refIDs1.has_key(query_id):
            refIDs1[query_id] = set()
        refIDs1[query_id].add(ref_id)

    if paired:
        for line in open(usearch2_fp, 'r'):
            words = line.strip().split('\t')
            query_id = '/'.join(words[0].split('/')[:-1])
            ref_id = words[1].split('_')[0]
            if not refIDs2.has_key(query_id):
                refIDs2[query_id] = set()
            refIDs2[query_id].add(ref_id)

        # for each query_id, keep only the intersection
        # of its R1 and R2 taxonomy assignments
        # if a query didn't have hits in both R1 and R2, throw it out
        print "Finding R2 intersection and taxonomy assignments..."
        for key in refIDs1:
            if refIDs2.has_key(key):
                common = commonTaxa(refIDs1[key], refIDs2[key], taxa)
                if common is not None:
                    taxon_assignments[key] = common
    else:
        print "Converting refIDs to taxonomy assignments..."
        for key in refIDs1:
            taxon_assignments[key] = set(
                [taxa[refID] for refID in refIDs1[key]])

    print "Determining consensus taxonomy assignments..."
    count = 1

    for key in taxon_assignments:
        if count % 1000000 == 0:
            print count
        count += 1

        taxa_set = taxon_assignments[key]
        taxon = taxa_set.pop()
        if len(taxa_set) > 0:
            for taxon2 in taxa_set:
                taxon = lca(taxon, taxon2)
        taxon_assignments[key] = taxon

    # tabulate taxon counts
    print count, "hits processed. Tabulating taxa..."
    taxon_counts = {}
    for key in taxon_assignments:
        taxon = taxon_assignments[key]
        if taxon is not None:
            if species_only:
                taxon = ' '.join(taxon)
            else:
                taxon = ';'.join(taxon)
            if not taxon_counts.has_key(taxon):
                taxon_counts[taxon] = 1
            else:
                taxon_counts[taxon] += 1

    # print results
    print "Writing results..."
    outf = open(out_fp, 'w')
    for key in sorted(taxon_counts):
        outf.write(key + '\t' + str(taxon_counts[key]) + '\n')
    outf.close()
