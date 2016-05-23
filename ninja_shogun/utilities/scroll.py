import os


class Scroll:
    def __init__(self, path):
        self.path = path

    def verify(self):
        return os.path.isfile(self.path)

    def __call__(self):
        return self.verify()


def scrolling(func):
    def scroll_wrapper(self, *args, **kwargs):
        if not self._scroll():
            raise Exception('Data file %s not found' % self._scroll.path)
        return func(self, *args, **kwargs)
    return scroll_wrapper
