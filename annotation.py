"""
Module allows loading of a tab-delimited file and allows some simple handling of the data.

We refer to such a file as an `Annotation' file. Annotation files are simply a tab-delimited file containing results or
annotations for use elsewhere. A header is required and the file must be tab-delimited.
"""
__author__ = 'julianzaugg'


class Annotation:

    def __init__(self, filename=None):
        self.header = []
        self.entries = dict()
        self.number_of_annotations = 0
        if filename:
            self.load_annotations(filename)

    def load_annotations(self, filename):
        """
        Load annotations from a file
        :param filename: filename
        """
        with open(filename, 'r') as fh:
            data = [line.strip().split("\t") for line in fh.readlines()]
            self.header = data[0]
            self.entries = dict([(c, []) for c in self.header])
            for line in data[1:]:
                for h_idx in range(len(self.header)):
                    self.entries[self.header[h_idx]].append(line[h_idx])
                self.number_of_annotations += 1

    def get_column(self, column_name):
        """
        Return entries for a specific column
        :param column_name: Name of header
        :return:
        """
        assert column_name in self.entries, "No column with the header %s in annotation" % column_name
        return self.entries[column_name]

    def add_column(self, column_name, entries):
        if self.number_of_annotations != 0:
            assert len(entries) == self.number_of_annotations, "Number of entry rows is less than existing annotation"
        self.entries[column_name] = entries
        self.header.append(column_name)

    def __str__(self):
        out = "\t".join(self.header) + "\n"
        for i in range(self.number_of_annotations):
            line_str = "\t".join([self.entries[column_name][i] for column_name in self.header])
            out += line_str + '\n'
        return out


if __name__ == "__main__":
    pass
