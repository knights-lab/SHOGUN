import unittest
import sys
from nose.tools import assert_equals

import os
from tempfile import TemporaryDirectory

from shogun.scripts import kegg_parse_img_ids


class KeggParseImgIDsTest(unittest.TestCase):
    def test(self):
        with TemporaryDirectory() as outdir:
            temp_file_name = os.path.join(outdir, 'test.name')
            sys.argv = ["kegg_parse_img_ids", "-i", ".\kegg_parse_img_ids", "-o", temp_file_name]
            kegg_parse_img_ids.main()
            assert_equals(None, None)
