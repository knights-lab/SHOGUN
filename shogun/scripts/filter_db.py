import re

# drop pseudo=True

with open("./scratch/refseq_small/combined_seq.fna") as inf:
    for line in inf:
        if line.startswith(">"):
            if "pseudo=true" in line:
                print("None")
            else:
                print(re.search('protein_id=(.*)}_{location', line).group(1))
            print(line)
