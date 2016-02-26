#!/usr/bin/env python
import urllib.request
import tarfile

from shogun import SETTINGS


def download_ncbi_taxdmp(_ncbi_taxdmp_url=SETTINGS.ncbi_taxdmp_url, _ncbi_taxdmp_dir=SETTINGS.ncbi_taxdmp_dir):
    req = urllib.request.Request(_ncbi_taxdmp_url)
    with urllib.request.urlopen(req, 'rb') as ftp_stream:
        tfile = tarfile.open(fileobj=ftp_stream, mode='r|gz')
        tfile.extractall(_ncbi_taxdmp_dir)


def main():
    download_ncbi_taxdmp()

if __name__ == '__main__':
    main()
