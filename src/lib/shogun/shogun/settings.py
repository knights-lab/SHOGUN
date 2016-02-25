from collections import namedtuple
import multiprocessing
import os

Settings = namedtuple('Settings', ['data_path', 'cache_path', 'results_path', 'ncbi_taxdmp_url', 'silva_taxdmp_urls',
                                   'N_jobs'])


def initialize_settings(path=os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..', '..', '..')):

    cache_path = os.path.join(path, 'cache')
    data_path = os.path.join(path, 'data')
    results_path = os.path.join(path, 'results')
    ncbi_taxdmp_url = 'ftp://ftp.ncbi.nih.gov:/pub/taxonomy/taxdump.tar.gz'
    silva_taxdmp_urls = [
        'http://www.arb-silva.de/fileadmin/silva_databases/release_119/Exports/taxonomy/taxmap_embl_ssu_parc_119.txt',
        'http://www.arb-silva.de/fileadmin/silva_databases/release_119/Exports/taxonomy/taxmap_embl_lsu_parc_119.txt']

    for path in [path, cache_path, data_path, results_path]:
        if not os.path.exists(path):
            os.makedirs(path)

    N_jobs = multiprocessing.cpu_count()

    return Settings(silva_taxdmp_urls=silva_taxdmp_urls, ncbi_taxdmp_url=ncbi_taxdmp_url, data_path=data_path,
                    cache_path=cache_path, results_path=results_path, N_jobs=N_jobs)
