__author__ = 'julianzaugg'


"""
Module contains methods and classes for creating, analysing and manipulating sequence objects.
"""

from collections import Counter
import numpy

from prob import *

class Sequence(object):

    def __init__(self, sequence, alphabet = None, name = '', info = '', **kwargs):
        for key, value in kwargs.items():
            setattr(self, key, value)
        self.sequence = sequence
        self.name = name
        self.alphabet = alphabet

        #Additional information
        self.length = len(sequence)
        self.info = info

    def __len__(self):
        return len(self.sequence)

    def __str__(self):
        out = "%s: %s" % (self.name, self.sequence)
        return out

    def __contains__(self, item):
        return item in self.sequence

    def __getitem__(self, ndx):
        return self.sequence[ndx]

    def __eq__(self, other):
        if isinstance(other, Sequence):
            return self.sequence == other.sequence
        return False

    def find(self, string):
        return self.sequence.find(string)

    def write_fasta(self):
        """ Write one sequence in FASTA format to a string and return it. """
        fasta = '>' + self.name + ' ' + self.info + '\n'
        data = self.sequence
        nlines = (len(self.sequence) - 1) / 60 + 1
        for i in range(nlines):
            lineofseq = ''.join(data[i*60 : (i+1)*60]) + '\n'
            fasta += lineofseq
        return fasta

class Alignment(object):

    def __init__(self, sequences):
        self.seqs = [s for s in sequences]
        self.seqs_dict = dict([(s.name, s) for s in self.seqs])
        self.alignlen = len(sequences[0])
        self.alphabet = self.seqs[0].alphabet

    def __len__(self):
        return len(self.seqs)

    def __getitem__(self, ndx):
        return self.seqs[ndx]

    def __contains__(self, item):
        return item.sequence in [s.sequence for s in self.seqs]

    def __str__(self):
        output = ""
        for seq in self.seqs:
            output += "%s\t%s\n" % (seq.name, seq.sequence)
        return output

    def add_sequence(self, sequence):
        self.seqs.append(sequence)
        self.seqs_dict = dict([(s.name, s) for s in self.seqs])

    def get_sequence(self, seq_name):
        try:
            return self.seqs_dict[seq_name]
        except KeyError:
            raise KeyError("Sequence %s was not found in the alignment" % seq_name)

    def get_probabilities(self, position, pseudo = 0, normalise = True):
        """
        Returns probabilities of each symbol in the alphabet at a position
        """
        colstr = "".join([seq[position] for seq in self.seqs])
        cnts = Counter(colstr)
        for sym in self.alphabet:
            if sym not in cnts and pseudo:
                cnts[sym] = pseudo
            else:
                cnts[sym] += pseudo
        total = sum(cnts.values())
        probs = dict()
        if not normalise: return cnts
        for sym, cnt in cnts.items():
            probs[sym] = float(cnt)/total
        return probs

    def write_clustal_file(self, filename):
        """
        Save a Alignment in CLUSTAL format.
        """
        symbolsPerLine = 60
        max_name_length = max(len(seq.name) for seq in self.seqs)
        namelen = 0
        string = ''
        for seq in self.seqs:
            namelen = max(len(seq.name), namelen)
        wholeRows = self.alignlen / symbolsPerLine
        for i in range(wholeRows):
            for j in range(len(self.seqs)):
                string += self.seqs[j].name.ljust(max_name_length) + ' '
                string += self.seqs[j][i * symbolsPerLine:(i + 1) * symbolsPerLine] + '\n'
            string += '\n'
        # Possible last row
        last_row_length = self.alignlen - wholeRows * symbolsPerLine
        if last_row_length > 0:
            for j in range(len(self.seqs)):
                if max_name_length > 0:
                    string += self.seqs[j].name.ljust(max_name_length) + ' '
                string += self.seqs[j][-last_row_length:] + '\n'
        if filename:
            fh = open(filename, 'w')
            # fake header so that clustal believes it
            fh.write('CLUSTAL O(1.2.0) multiple sequence alignment\n\n\n')
            fh.write(string)
            fh.close()
            return
        return string

    def get_profile(self, pseudo = 0.0):
        """ Determine the probability matrix from the alignment, assuming
        that each position is independent of all others. """
        p = IndepJoint([self.alphabet for _ in range(self.alignlen)], pseudo)
        for seq in self.seqs:
            p.observe(seq)
        return p

    def get_ungapped(self):
        """
        Return new alignment with gappy columns removed
        """
        gappy_columns = set()
        for seq in self:
            for i in range(self.alignlen):
                if i not in gappy_columns and seq[i] == "-":
                    gappy_columns.add(i)
        new_seqs = []
        for seq in self:
            content = "".join([seq[i] for i in range(self.alignlen) if i not in gappy_columns])
            new_seqs.append(Sequence(sequence=content, alphabet=seq.alphabet, name=seq.name))
        return Alignment(new_seqs)

    def get_ungapped_using_reference(self, seq_name):
        """
        Return a new alignment where gappy columns have been removed using in respect to
        a user specified reference sequence
        :param name: Name of template sequence
        :return:
        """
        template_seq = self.get_sequence(seq_name)
        #Find gapped column indices
        gappy_columns = [i for i in xrange(len(template_seq)) if template_seq[i] == "-"]
        new_seqs = []
        for seq in self:
            content = "".join([seq[i] for i in range(self.alignlen) if i not in gappy_columns])
            new_seqs.append(Sequence(sequence=content, alphabet=seq.alphabet, name=seq.name))
        return Alignment(new_seqs)

    def get_column(self, position):
        return [s[position] for s in self]


def read_fasta_file(filename, alphabet):
    """
    Read a Fasta file and return a set of Sequence
    {name: seq_string}
    """
    fh = open(filename, 'r')
    seqdata = dict()
    order = []
    data = [line.strip() for line in fh.readlines() if line is not None]

    for line in data:
        if not line: continue
        if line[0] == '>':
            name = line.split()[0][1:]
            order.append(name)
            seqdata[name] = ''
        else:
            seqdata[name] += line
    output = []
    for sname in order:
        sseq = seqdata[sname]
        output.append(Sequence(name=sname, sequence=sseq.strip(), alphabet=alphabet))
    return output

def read_clustal_file(filename, alpha):
    """
    Read a CLUSTAL Alignment file and return a dictionary of sequence data
    """
    fh = open(filename, 'r')
    names = []
    seqdata = dict()
    data = [line.strip('\n') for line in fh.readlines() if line is not None]
    for line in data:
        if line.startswith('CLUSTAL') or line.startswith('#'):
            continue
        if len(line) == 0:
            continue
        if line[0] == ' ' or '*' in line or ':' in line:
            continue
        sections = line.split()
        name, seqstr = sections[0], "".join(sections[1:])
        names.append(name)
        if seqdata.has_key(name):
            seqdata[name] += seqstr
        else:
            seqdata[name] = seqstr
    sequences = [Sequence(seqstr, name=seqname, alphabet=alpha) for seqname, seqstr in sorted(seqdata.items(),
                                                                            key=lambda x: names.index(x[0]))]
    return Alignment(sequences)

def write_fasta_file(filename, seqs):
    """ Write the specified sequences to a FASTA file. """
    fh = open(filename, 'w')
    for seq in seqs:
        fh.write(seq.write_fasta())
    fh.close()


class PWM(object):

    """ A position weight matrix. """

    def __init__(self, foreground, start = 0, end = None, pseudo = 0.0):
        """ Create a new PWM from the given probability matrix/ces.
        foreground: can be either an Alignment, a list of Distrib's or an instance of IndepJoint.
        background: must be a Distrib instance or None (in which case a uniform background will be used)
        Specify only a section of the matrix to use with start and end. """
        # print start,end, foreground[0].name, len(foreground), foreground.alignlen, foreground
        if isinstance(foreground, Alignment):
            foreground = foreground.get_profile(pseudo=pseudo)
        if isinstance(foreground, IndepJoint):
            foreground = foreground.store
        self.start = start
        self.end = end or len(foreground)
        self.length = self.end - self.start
        self.alphabet = foreground[self.start].alpha
        if False in [ col.alpha == self.alphabet for col in foreground[self.start + 1 : self.end] ]:
            raise RuntimeError("All positions need to be based on the same alphabet")
        self.symbols = self.alphabet.symbols
        # Set foreground probabilities from given alignment
        self.m = numpy.zeros((len(self.symbols), self.length))
        self.fg = foreground[self.start:self.end]
        self.bg = None or Distrib(self.alphabet, 1.0)
        if not self.alphabet == self.bg.alpha:
            raise RuntimeError("Background needs to use the same alphabet as the foreground")
        p = self.bg.prob()
        for i in range(self.length):
            q = self.fg[i].prob()
            for j in range(len(self.alphabet)):
                self.m[j][i] = self.logme(q[j], p[j])

    def __len__(self):
        return self.length

    MIN_VALUE = 0.00000000001

    def logme(self, fg, bg):
        if fg > self.MIN_VALUE and bg > self.MIN_VALUE:
            ratio = fg / bg
            return math.log(ratio)
        # if not, one of fg and bg is practically zero
        if fg > self.MIN_VALUE: # bg is zero
            return math.log(fg / self.MIN_VALUE)
        else: # fg is zero
            return math.log(self.MIN_VALUE)

    def get_matrix(self):
        return self.m

    def __str__(self):
        str = ''
        for j in range(len(self.alphabet)):
            str += "%s\t%s\n" % (self.alphabet[j], ' '.join("%+6.2f" % (y) for y in self.m[j]))
        return str

    def display(self, format='COLUMN'):
        if format == 'COLUMN':
            print " \t%s" % (' '.join(" %5d" % (i + 1) for i in range(self.length)))
            for j in range(len(self.alphabet)):
                print "%s\t%s" % (self.alphabet[j], ' '.join("%+6.2f" % y for y in self.m[j]))
        elif format == 'JASPAR':
            for j in range(len(self.alphabet)):
                print "%s\t[%s]" % (self.alphabet[j], ' '.join("%+6.2f" % y for y in self.m[j]))

    def search(self, sequence, lowerBound=0):
        """ Find matches to the motif in a specified sequence. Returns a list
        of  results as triples: (position, matched string, score).
        The optional argument lowerBound specifies a lower bound on reported
        scores. """
        results = []
        for i in range(len(sequence)-self.length+1):
            subseq = sequence[i:i + self.length]
            ndxseq = [self.alphabet.index(sym) for sym in subseq]
            score = 0.0
            for w in range(len(ndxseq)):
                score += self.m[ndxseq[w]][w]
            if score > lowerBound:
                results.append((i, subseq, score))
        return results

    def maxscore(self, sequence):
        """ Find matches to the motif in a specified sequence.
            Returns the maximum score found in the sequence and its index as a tuple:
            (maxscore, maxindex) """
        maxscore = 0.0
        maxindex = 0.0
        for i in range(len(sequence)-self.length+1):
            subseq = sequence[i:i + self.length]
            ndxseq = [self.alphabet.index(sym) for sym in subseq]
            score = 0.0
            for w in range(len(ndxseq)):
                score += self.m[ ndxseq[w] ][ w ]
            if maxscore is None:
                maxscore = score
                maxindex = i
            elif maxscore < score:
                maxscore = score
                maxindex = i
        return (maxscore, maxindex)
