"""
Copyright 2015-2017 Knights Lab, Regents of the University of Minnesota.

This software is released under the GNU Affero General Public License (AGPL) v3.0 License.
"""

import os


class Aligner:
    def __init__(self, database_dir):
        check, msg = self.check_database(self.__class__, database_dir)

        if not check:
            raise Exception("Database %s is not formatted correctly: %s" % (database_dir, msg))

        self.database_dir = database_dir

    @classmethod
    def check_database(cls):
        raise NotImplementedError


class EmbalmerAligner(Aligner):
    def __init__(self, database_dir):
        super().__init__(database_dir)

    @classmethod
    def check_database(cls, dir):
        if not os.path.exists(dir):
            return False, "Directory does not exist"


class UtreeAligner(Aligner):
    def __init__(self):
        super().__init__()


class BowtieAligner(Aligner):
    def __init__(self):
        super().__init__()
