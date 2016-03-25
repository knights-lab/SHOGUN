#!/usr/bin/env python
import os
import pickle

from shogun import SETTINGS, LOGGER


class Pickleable:
    def __init__(self, _pickle_dir=SETTINGS.pickle_dir):
        self._pickle_dir = _pickle_dir
        self.dump = os.path.join(self._pickle_dir, "%s.pkl" % self.__class__.__name__)
        saved = os.path.isfile(self.dump)
        if not saved:
            self._parse()
            self.save()
        else:
            self.load()

    def _parse(self):
        pass

    def save(self):
        self_dump = os.path.join(self._pickle_dir, "%s.pkl" % self.__class__.__name__)
        with open(self_dump, 'wb') as handle:
            pickle.dump(self, handle)

    @classmethod
    def load(cls, _pickle_dir=os.path.join(SETTINGS.data_dir, 'pickle')):
        try:
            self_dump = os.path.join(_pickle_dir, "%s.pkl" % cls.__name__)
            with open(self_dump, 'rb') as handle:
                return pickle.load(handle)
        except FileNotFoundError as error:
            LOGGER.log('Failed to open the %s.pkl' % cls.__name__)
            raise error
