#!/usr/bin/env python
import os

from shogun.downloaders.download_txt_file import download_txt_url
from shogun import SETTINGS


def download_silva_mapping_file(_silva_taxdmp_urls=SETTINGS.silva_taxdmp_urls,
                                _silva_taxdmp_dir=SETTINGS.silva_taxdmp_dir):
        for url in _silva_taxdmp_urls:
            file_name = url.split('/')[-1]
            download_txt_url(os.path.join(_silva_taxdmp_dir, file_name), url)


def main():
    download_silva_mapping_file()

if __name__ == '__main__':
    main()
