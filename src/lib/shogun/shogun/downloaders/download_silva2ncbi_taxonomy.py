#!/usr/bin/env python
import os

from shogun.utilities.downloadable import Downloadable
from shogun.downloaders.utilities.stream_url import download_txt_url
from shogun import SETTINGS


class SilvaMapping(Downloadable):
    def __init__(self, _silva_taxdmp_urls=SETTINGS.silva_taxdmp_urls, _silva_taxdmp_dir=SETTINGS.silva_taxdmp_dir):
        super().__init__(_silva_taxdmp_dir)
        self.urls = _silva_taxdmp_urls

    def download(self):
        for url in self.urls:
            file_name = url.split('/')[-1]
            download_txt_url(os.path.join(self.path, file_name), url)


def main():
    SilvaMapping().run()

if __name__ == '__main__':
    main()
