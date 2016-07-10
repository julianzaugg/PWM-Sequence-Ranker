__author__ = 'julianzaugg'


"""
PWM SEQUENCE RANKER
Given input of:
    - A Clustal formatted sequence alignment
    - A set of sequences to use as a background (Fasta)
    - A query sequence/s to rank and score (Fasta)

This package will construct a PWM from the user supplied alignment, and will
score and save the the match scores for each background sequence. Each query
sequence will be scored against the PWM and ranked against the background distribution
of scores.
"""

import argparse
from collections import defaultdict, Counter
from itertools import groupby
from operator import itemgetter
import os
import annotation
import string
import random

from sequence import *


QUERY_SEQS = None
QUERY_SEQS_DICT = dict()  # FIXME - change QUERY_SEQS to a dict so we don't have to have this as well
# Group/Alignment name : {Sequence Name : Score}
QUERY_SEQ_SCORES_ALL = defaultdict(lambda: defaultdict(float))

# Sequence Name : (Group/Alignment name, Rank, Score)}
QUERY_SEQ_BEST = defaultdict(lambda: (None, 0, 100000000000))

# Group/Alignment name : {Sequence Name : Score}
QUERY_SEQ_RANKS = defaultdict(lambda: defaultdict(int))

BACKGROUND_SEQS = None
BACKGROUND_SEQS_DICT = dict()  # FIXME - change BACKGROUND_SEQS to a dict so we don't have to have this as well
N_BG_SEQS = None

# Group/Alignment name : {Sequence Name : Score}
BACKGROUND_SCORES = defaultdict(lambda: defaultdict(float))

# Group/Alignment name : Score
ORDERED_BACKGROUND_SCORES = defaultdict()

# Group/Alignment name : Alignment
PWM_ALIGNMENTS = defaultdict()

# Group/Alignment name : PWM
PWMS = defaultdict()

# Group/Alignment name : Start_idx, End_idx
ALIGNMENT_TRIM_INDICES = defaultdict()

LONGEST_ALN_LENGTH = 0

###############################################################################
#                                  PWM
def _create_pwms(pseudo_count, trim_threshold):
    global PWMS, ALIGNMENT_TRIM_INDICES, LONGEST_ALN_LENGTH
    for name, alignment in PWM_ALIGNMENTS.iteritems():
        aln_n_seqs = len(alignment)
        if alignment.alignlen > LONGEST_ALN_LENGTH:
            LONGEST_ALN_LENGTH = alignment.alignlen
        if aln_n_seqs == 1 or trim_threshold == 1.0:
            trim_start, trim_end = 0, None
        else:
            trim_start, trim_end = _find_trim_idxs(alignment, trim_threshold)
        ALIGNMENT_TRIM_INDICES[name] = trim_start, trim_end
        PWMS[name] = PWM(alignment, pseudo=pseudo_count, start=trim_start, end=trim_end)

def _find_trim_idxs(an_alignment, gap_threshold):
    """
    Calculate the indices/positions of columns in the pwm alignment that have
    too many gaps (determined by the gap_threshold). Specifically identifies
    the continuous sets of column indices within the full list of gappy columns and
    will return the highest of the first set and the lowest of the last set.
    Example:
    bad_columns = [0,1,2,3,4,7,8,9]
    bad_column_sets = [[0,1,2,3,4], [7,8,9]]
    return values = 4,7
    :param an_alignment: a sequence alignments
    :param gap_threshold: float between 0 and 1, determines percentage of gaps allowed in column
    :return: index_low, index_high
    """
    bad_columns = []
    for i in xrange(an_alignment.alignlen):
        column_vals = an_alignment.get_column(i)
        counts = Counter(column_vals)
        # Calculate percentage of gaps in column
        if "-" in counts and counts["-"]/float(len(an_alignment)) >= gap_threshold:
            bad_columns.append(i)
    if len(bad_columns) == 0:
        return 0, None
    bad_columns_split = []
    """
    Identify the separate bad column sets. The way this works is
    we minus the column index (identified above) from the enumeration
    index value. Since the enumeration starts from 0, each bad column will be
    a certain negative value difference from the enumeration index. Consecutive
    runs of bad columns will have the same difference to their respective enumeration indices.
    For example, if bad_columns = [1,2,3,4,7,8,9,11,12,15,19], the groups we get are:
        -1 [(0, 1), (1, 2), (2, 3), (3, 4)]
        -3 [(4, 7), (5, 8), (6, 9)]
        -4 [(7, 11), (8, 12)]
        -6 [(9, 15)]
        -9 [(10, 19)]
    itemgetter(1) simply pulls out the 1th value from the above tuples for each group (g) and
    returns a list of the values:
        -1 [1, 2, 3, 4]
        -3 [7, 8, 9]
        -4 [11, 12]
        -6 [15]
        -9 [19]
    """
    for f, g in groupby(enumerate(bad_columns), lambda (x, y): x-y):
        bad_columns_split.append(map(itemgetter(1), g))
    if len(bad_columns_split) == 1:
        # If there is only 1 consecutive run of bad columns, determine if it is in the first or
        # second half of the alignment
        if bad_columns_split[0][0] >= .5 * an_alignment.alignlen:
            return 0, bad_columns_split[0][0] + 1
        return bad_columns_split[0][-1] + 1, None
    # TODO - investigate whether we need to +1 on start AND -1 on end? I think just +1 to start.
    return bad_columns_split[0][-1] + 1, bad_columns_split[-1][0]

###############################################################################
#                                  LOADING

def _load_background_score_file(location):
    """
    Assumes the input file has columns in header of - Name, Group, Score. Where
    Name is the sequence name, Group the alignment file name it was scored against, and Score the PWM match score
    :param location: location of the input file
    """
    global BACKGROUND_SCORES, N_BG_SEQS
    bgs_annotation = annotation.Annotation(location)
    groups = bgs_annotation.get_column("Group")
    names = bgs_annotation.get_column("Name")
    scores = bgs_annotation.get_column("Score")
    N_BG_SEQS = len(set(names))
    for i in xrange(bgs_annotation.number_of_annotations):
        BACKGROUND_SCORES[groups[i]][names[i]] = scores[i]

def _load_alignments(location):
    global PWM_ALIGNMENTS
    if os.path.isdir(location):
        aln_filenames = [n for n in os.listdir(location) if n.endswith(".txt")]
        for filename in aln_filenames:
            if filename.endswith(".txt"):
                PWM_ALIGNMENTS[filename.split(".")[0]] = read_clustal_file(location + filename, Protein_Alphabet)
    elif os.path.isfile(location):
        name = location.split("/")[-1].split(".")[0]
        PWM_ALIGNMENTS[name] = read_clustal_file(location, Protein_Alphabet)
    else:
        raise StandardError("Error Loading PWM alignments")

###############################################################################
#                                  SCORING
def _score_background_sequences():
    global BACKGROUND_SCORES
    cnt = 0
    for bg_seq in BACKGROUND_SEQS:
        try:
            print "Scoring background sequence %i\t%s" % (cnt, bg_seq.name)
            for pwm_name, pwm in PWMS.iteritems():
                match_score, match_index = pwm.maxscore(bg_seq)
                BACKGROUND_SCORES[pwm_name][bg_seq.name] = match_score
            cnt += 1
        except Exception:
            print "Sequence %s cannot be processed" % bg_seq.name


def _score_query_sequences():
    global QUERY_SEQ_SCORES_ALL, QUERY_SEQ_SCORES_BEST, ORDERED_BACKGROUND_SCORES, QUERY_SEQ_RANKS
    for q_seq in QUERY_SEQS:
        print "Scoring query sequence %s" % q_seq.name
        for pwm_name, pwm in PWMS.iteritems():
            match_score, match_index = pwm.maxscore(q_seq)
            QUERY_SEQ_SCORES_ALL[pwm_name][q_seq.name] = match_score
            # Now calculate rank vs the background
            bg_scores = ORDERED_BACKGROUND_SCORES[pwm_name]
            rank_index = 0
            while rank_index < int(N_BG_SEQS):
                # FIXME? - why do these have to be converted to floats for this statement to work?
                # if float(bg_scores[rank_index]) <= float(match_score):
                if float(bg_scores[rank_index]) < float(match_score):
                    # Best rank found (since ordered we can stop at first instance)
                    break
                rank_index += 1
            QUERY_SEQ_RANKS[pwm_name][q_seq.name] = rank_index
            if rank_index < QUERY_SEQ_BEST[q_seq.name][2] and match_score != 0.0:
                QUERY_SEQ_BEST[q_seq.name] = (pwm_name, match_score, rank_index)
###############################################################################
#                                  SAVING

def _save_background_scores(location):
    with open(location, 'w') as fh:
        print >> fh, "Name\tGroup\tScore\tSequence"
        for group_name, b_seq_dict in BACKGROUND_SCORES.iteritems():
            for b_seq_name, b_seq_score in b_seq_dict.iteritems():
                print >> fh, "%s\t%s\t%f\t%s" % (b_seq_name, group_name, b_seq_score,
                                             BACKGROUND_SEQS_DICT[b_seq_name].sequence)

def _save_query_scores(location):
    global QUERY_SEQ_RANKS, QUERY_SEQ_SCORES_ALL
    with open(location, 'w') as fh:
        print >> fh, "Name\tGroup\tScore\tRank\tSequence"
        for group_name, q_seq_dict in QUERY_SEQ_SCORES_ALL.iteritems():
            for q_seq_name, score in q_seq_dict.iteritems():
                print >> fh, "%s\t%s\t%f\t%i\t%s" % (q_seq_name, group_name, score,
                                                 QUERY_SEQ_RANKS[group_name][q_seq_name],
                                                 QUERY_SEQS_DICT[q_seq_name].sequence)

def _save_best_query_scores(location):
    with open(location, 'w') as fh:
        print >> fh, "Name\tGroup\tScore\tRank\tSequence"
        for q_seq_name, best_results in QUERY_SEQ_BEST.iteritems():
            print >> fh, "%s\t%s\t%f\t%i\t%s" % (q_seq_name, best_results[0], best_results[1],
                                             QUERY_SEQ_RANKS[best_results[0]][q_seq_name],
                                             QUERY_SEQS_DICT[q_seq_name].sequence)

def _save_trim_indices(location):
    with open(location, 'w') as fh:
        print >> fh, "Group\tStart\tEnd"
        for group_name, idxs in ALIGNMENT_TRIM_INDICES.iteritems():
            print >> fh, "%s\t%s\t%s" % (group_name, str(idxs[0]), str(idxs[1]))

###############################################################################
#                                  MISC
def _is_valid_file(parser, arg):
    """Check if arg is a valid file that already exists on the file
       system.
    """
    arg = os.path.abspath(arg)
    if not os.path.exists(arg):
        parser.error("The file %s does not exist!" % arg)
    else:
        return arg

def _dir_check(directory):
    if not os.path.exists(directory):
        raise StandardError("The directory %s does not exist or the location was incorrectly specified" % directory)


def _generateRandomStringID(length=5, alpha=None, seed=None):
    """
    @param length: length of string
    @return: string
    """
    if seed:
        random.seed(seed)
    if not alpha:
        alpha = string.ascii_uppercase + string.digits
    return ''.join(random.choice(alpha) for x in range(length))

def _generate_random_sequences(amount):
    seeds = range(1, amount+1)  # seed can't equal 0
    random_sequences = [Sequence(_generateRandomStringID(LONGEST_ALN_LENGTH,Protein_Alphabet, seed=i),
                                 name=_generateRandomStringID(8, seed=i), alpha=Protein_Alphabet) for i in seeds]
    return random_sequences

def _log_parameters(location, args):
    with open(location, 'w') as fh:
        print >> fh, "#Parameters used"
        for arg, value in sorted(vars(args).items()):
            print >> fh, "Argument %s: %r" % (arg, value)

def _check_sequences(seqs):
    good = []
    for seq in seqs:
        bad = False
        for aa in seq:
            if aa not in Protein_Alphabet:
                print "Sequence %s cannot be processed, it will be ignored" % seq.name
                bad = True
                break
        if bad:
            continue
        good.append(seq)
    return good

def process_args(my_parser, my_args):
    global QUERY_SEQS, BACKGROUND_SEQS, BACKGROUND_SCORES, ORDERED_BACKGROUND_SCORES, N_BG_SEQS
    # Check files
    _dir_check(my_args.output)
    _is_valid_file(my_parser, my_args.query)
    _is_valid_file(my_parser, my_args.alignments)
    # Save the parameter values
    _log_parameters(my_args.output + "parameter_log.txt", my_args)
    # Load query sequence/s
    print "Loading Query sequences"
    QUERY_SEQS = read_fasta_file(my_args.query, Protein_Alphabet)
    for q_seq in QUERY_SEQS:
        QUERY_SEQS_DICT[q_seq.name] = q_seq
    print "Loading Alignments"
    _load_alignments(my_args.alignments)
    print "Creating PWMS"
    _create_pwms(my_args.pwm_pseudo, my_args.trim_threshold)
    print "Saving trim indices"
    _save_trim_indices(my_args.output + "pwm_trim_locations.txt")
    # Load background sequences or, if score file supplied, load the scores
    if my_args.backgroundseqs:
        BACKGROUND_SEQS = read_fasta_file(my_args.backgroundseqs, Protein_Alphabet)  # also acts a file check
        print "Checking sequences"
        BACKGROUND_SEQS = _check_sequences(BACKGROUND_SEQS)
        for bg_seq in BACKGROUND_SEQS:
            BACKGROUND_SEQS_DICT[bg_seq.name] = bg_seq
        N_BG_SEQS = len(BACKGROUND_SEQS)
        print "Scoring background sequences (this might take several hours depending on length and number of sequences)"
        _score_background_sequences()
        _save_background_scores(my_args.output + "bg_score_file.txt")
    elif my_args.backgroundfile:
        print "Loading background score file"
        _load_background_score_file(my_args.backgroundfile)
    elif my_args.backgroundnumber:
        print "Generating %i random sequences of length %i" % (my_args.backgroundnumber, LONGEST_ALN_LENGTH)
        BACKGROUND_SEQS = _generate_random_sequences(my_args.backgroundnumber)
        for bg_seq in BACKGROUND_SEQS:
            BACKGROUND_SEQS_DICT[bg_seq.name] = bg_seq
        N_BG_SEQS = my_args.backgroundnumber
        _score_background_sequences()
        _save_background_scores(my_args.output + "bg_score_file.txt")
    # We store just the ordered scores for each group to scan across
    for group_name, bg_seq_scores_dict in BACKGROUND_SCORES.iteritems():
        ORDERED_BACKGROUND_SCORES[group_name] = sorted(bg_seq_scores_dict.values())[::-1]

    # Score each query sequence
    print "Scoring and ranking query sequences"
    _score_query_sequences()
    print "Saving query results"
    _save_query_scores(my_args.output + "query_all_scores.txt")
    _save_best_query_scores(my_args.output + "query_best_scores.txt")

def _restricted_trim(x):
    x = float(x)
    if x < 0.0 or x > 1.0:
        raise argparse.ArgumentTypeError("%r not in range [0.0, 1.0]" % (x,))
    return x

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Description to be added')
    parser.add_argument('-alns', '--alignments', help='Alignments for PWM creation, currently must be clustal '
                                                            'formatted multiple sequence alignment', required=True)
    bggroup = parser.add_mutually_exclusive_group(required=True)
    bggroup.add_argument('-bgs', '--backgroundseqs', help='Background sequences, Fasta')
    bggroup.add_argument('-bgf', '--backgroundfile', help='Background sequence score file')
    bggroup.add_argument('-bgn', '--backgroundnumber', help='Background sequences to generate', type=int, default=1000)
    parser.add_argument('-q', '--query', help="Query sequences, Fasta", required=True)
    parser.add_argument('-o', '--output', help='Output location,', required=False, default="./")
    parser.add_argument('-tt', '--trim_threshold', help='Trim columns from the PWM alignments that have gaps over '
                                                        'the specified percentage, must be 0-1', required=False,
                                                        type=_restricted_trim, default=1.0)
    parser.add_argument('-pp', '--pwm_pseudo', help='Pseudo count for PWM', type=float, required=False, default=0.0)
    # parser.add_argument('-v', '--verbose', required=False, default=False)
    args = parser.parse_args()
    process_args(parser, args)


