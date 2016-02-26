from collections import namedtuple
import multiprocessing
import os
import json

from shogun.utils.path import verify_make_path


def initialize_settings(lib_path=os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                                          '..', '..', '..', '..'))):
    keys = ('lib_path', 'data_path', 'results_path', 'ncbi_taxdmp_url', 'silva_taxdmp_urls', 'N_jobs')
    Settings = namedtuple('Settings', keys)

    default_values = (
        lib_path,
        os.path.join(lib_path, 'data'),
        os.path.join(lib_path, 'results'),
        'ftp://ftp.ncbi.nih.gov:/pub/taxonomy/taxdump.tar.gz',
        ['http://www.arb-silva.de/fileadmin/silva_databases/release_119/Exports/taxonomy/taxmap_embl_ssu_parc_119.txt',
         'http://www.arb-silva.de/fileadmin/silva_databases/release_119/Exports/taxonomy/taxmap_embl_lsu_parc_119.txt'],
        multiprocessing.cpu_count()
    )

    settings_dict = dict(zip(keys, default_values))

    if os.path.exists(os.path.join(lib_path, 'SETTINGS.json')):
        with open(os.path.join(lib_path, 'SETTINGS.json')) as inf_handle:
            j = json.load(inf_handle)
            for key in j.keys():
                if j[key] != 'DEFAULT' and key in settings_dict:
                    settings_dict[key] = j[key]

    for path in ['data_path', 'results_path']:
        verify_make_path(settings_dict[path])

    return Settings(**settings_dict)


def main():
    initialize_settings()

if __name__ == '__main__':
    main()
