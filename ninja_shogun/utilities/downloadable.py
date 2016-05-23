#!/usr/bin/env python
import json
import os.path
import time

from shogun.utilities.scroll import Scroll


class Downloadable(Scroll):
    def __init__(self, path):
        super().__init__(path)
        self.data = None
        self.urls = []

    def download(self):
        pass

    def run(self):
        if not self.verify():
            self.download()
            self.save()

    def save(self):
        scroll = os.path.join(self.path, 'scroll.json')
        json_dict = dict(zip(('urls', 'date'), (self.urls, time.strftime("%d/%m/%Y"))))
        with open(scroll, 'w') as outf:
            json.dump(json_dict, outf)

    def verify(self):
        scroll = os.path.join(self.path, 'scroll.json')
        try:
            if os.path.isfile(scroll):
                with open(scroll) as inf:
                    info = json.load(inf)
                self.urls = info['urls']
                self.date = info['date']
                return True
        except KeyError:
            return False

    def __call__(self):
        self.run()


def download(func):
    def download_wrapper(self, *args, **kwargs):
        self._downloader()
        return func(self, *args, **kwargs)
    return download_wrapper
