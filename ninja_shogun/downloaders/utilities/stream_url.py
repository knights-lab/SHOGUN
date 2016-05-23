#!/usr/bin/env python
import urllib.request


def download_txt_url(path_to_file, url):
    with urllib.request.urlopen(url) as stream:
        with open(path_to_file, 'wb') as outfile:
            for line in stream:
                outfile.write(line)
