class Fasta:
    def __init__(self, inf):
        """

        :param inf: str
        :return:
        """
        self.inf = inf

    def read(self):
        """

        :param f: the FASTA file
        :return: tuples of (title, seq)
        """
        title = None
        data = None
        for line in f:
            if line[0] == ">":
                if title:
                    yield (title.strip(), data)
                title = line[1:]
                data = ''
            else:
                data += line.strip()
        if not title:
            yield None
        yield (title.strip(), data)
