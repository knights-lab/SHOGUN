import pandas as pd

def main():
    df = pd.read_csv('idxstats.csv', sep='\t', header=None, engine='c')
    df.columns = ['sequence_name', 'sequence_length', 'num_map_reads', 'num_umap_reads']

if __name__ == "__main__":
    main()
