from collections import namedtuple
import multiprocessing
import os
import json

from shogun.utils.path import verify_make_dir

Settings = namedtuple('Settings', ('default_dir', 'docs_dir', 'data_dir', 'results_dir', 'pickle_dir',
                                   'ncbi_taxdmp_url', 'ncbi_taxdmp_dir', 'silva_taxdmp_urls', 'N_jobs'))


def initialize_settings(config_dir=os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                                                '..', '..', '..', '..', '..'))):
    # I wish I could synchronize these with the Settings namedtuple but I can't
    default_keys = (
        'default_dir',
        'docs_dir',
        'data_dir',
        'results_dir',
        'pickle_dir',
        'ncbi_taxdmp_url',
        'ncbi_taxdmp_dir',
        'silva_taxdmp_urls',
        'N_jobs'
    )

    if os.path.exists(os.path.join(config_dir, 'SETTINGS.json')):
        with open(os.path.join(config_dir, 'SETTINGS.json')) as inf_handle:
            j = json.load(inf_handle)
            # result the settings dict if user changed default_dir
            if 'default_dir' in j:
                default_values = make_default_values(j['default_dir'])
                j.pop('default_dir', None)
                settings_dict = dict(zip(default_keys, default_values))

            for key in j.keys():
                if key in settings_dict:
                    settings_dict[key] = j[key]
    else:
        default_values = make_default_values(config_dir)
        settings_dict = dict(zip(default_keys, default_values))

    for outdir in [item for item in default_keys if 'dir' in item]:
        verify_make_dir(settings_dict[outdir])

    return Settings(**settings_dict)


def make_default_values(default_dir):
    return (
        default_dir,
        os.path.join(default_dir, 'docs'),
        os.path.join(default_dir, 'data'),
        os.path.join(default_dir, 'results'),
        os.path.join(default_dir, 'data', 'pickle'),
        'ftp://ftp.ncbi.nih.gov:/pub/taxonomy/taxdump.tar.gz',
        os.path.join(default_dir, 'data', 'ncbi_taxdmp'),
        ['http://www.arb-silva.de/fileadmin/silva_databases/release_119/Exports/taxonomy/taxmap_embl_ssu_parc_119.txt',
         'http://www.arb-silva.de/fileadmin/silva_databases/release_119/Exports/taxonomy/taxmap_embl_lsu_parc_119.txt'],
        multiprocessing.cpu_count()
    )


def main():
    settings = initialize_settings()
    print(settings)

if __name__ == '__main__':
    main()
