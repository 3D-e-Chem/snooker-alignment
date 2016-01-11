#!/usr/bin/env python

import copy
import itertools
import math
import os
import socket
import string
import sys
import snooker


# THE TEXT HEADER
def txt_header():
    print ""
    print "**********************************************************"
    print "*                                                        *"
    print "*                A L I G N M E N T .py                   *"
    print "*                                                        *"
    print "**********************************************************"
    print "*          Written in 2007 by Sven van den Beld          *"
    print "*          MDI, Organon, Oss, The Netherlands            *"
    print "**********************************************************"
    print "*              sven.vandenbeld@organon.com               *"
    print "**********************************************************"
    print "**********************************************************"
    print "*   Acknowlegdements:   Jan Klomp,  Marijn Sanders       *"
    print "**********************************************************"
    print ""


#  ======================================================================
#    M A C H I N E   S P E C I F I C   S E T T I N G S
#  ======================================================================

hostname = socket.gethostname()

print hostname

# SET PATH TO MODULES ON RELEVANT MACHINES
if hostname in ["acacia.oss.intra"]:
    sys.path.insert(0, "/home/belds/zooi/")
    import logging as log

#  ======================================================================
#    E R R O R   A N D   D E B U G G I N G   F U N C T I O N   G R O U P
#  ======================================================================

# DISPLAY ERROR MESSAGE AND STOP SCRIPT
# =====================================


def err_error(err):
    log.error(err)
    print
    sys.exit(1)

# SKIP AN ERROR
# =============


def err_skiperror(err):
    text = "Skipping error - %s." % err
    log.error(text)

# LIST ERRORS
# ===========


def err_listerrors(logfile, verbose=True):
    errors, warnings = [], []
    # READ LOGFILE
    print "**********************************************************"
    print "* Logfile summary: %37s *" % logfile
    print "**********************************************************"
    content = open(logfile, "r").readlines()
    for line in content:
        if line.find("ERROR") != -1:
            errors.append(line[:-1])
        if line.find("WARNING") != -1:
            warnings.append(line[:-1])
    # GIVE OUTPUT
    print "* - Reported %3i warnings - *" % len(warnings)
    if verbose:
        for line in warnings:
            print line
    print "* - Reported %3i errors   - *" % len(errors)
    if verbose:
        for line in errors:
            print line
    print "*****************************"

#  ======================================================================
#                            F U N C T I O N S
#  ======================================================================

# retrieves the variables (like human alignment location) from file


def get_variables(variables):
    input = open(variables).readlines()
    for line in input:
        line = line.strip()
        line = line.split()
        if line:
            if line[0] == "human_alignment":
                human_alignment = line[1]
            if line[0] == "cut_off":
                cut_off = {}
                for TM in range(1, 7 + 1):
                    cut_off[TM] = float(line[TM])
            if line[0] == "allowed_low_TMs":
                allowed_low_TMs = int(line[1])
            if line[0] == "calc_treshold":
                if line[1] == "True":
                    calc_treshold = True
                else:
                    calc_treshold = False
            if line[0] == 'number_format':
                numberings = snooker.Numberings.from_filename(line[1])
    return human_alignment, cut_off, allowed_low_TMs, calc_treshold, numberings
# reads the human_alignment file and stores the information in a dictionary


def read_human_alignment(human_alignment):
    # create an empty dictionary to store human alignment information
    human_alignment_dict = {}
    # open alignment file file
    content = open(human_alignment, "r").readlines()
    # create an empty dictionary
    TMdict = {}
    for line in content:
        line = line.strip()
        line = line.split()
        if len(line) > 0:
            # get the lengths of all TM's
            if line[0] == "name":
                TMlengths = {}
                length = 0
                TM = 1
                linelength = len(line)
                for column_nr in range(1, linelength):
                    # a TM is indicated with a T and number seperated by an "-"
                    if line[column_nr][0] == "T":
                        length += 1
                    else:
                        TMlengths[TM] = length
                        length = 0
                        TM += 1
                # and when the end of the file is reached
                TMlengths[TM] = length
            # a line with all non-olfactory class A GPCR on it starts with >
            # and has an R in column 2 (R represents the class of the GPCR)
            if line[0][0] == ">" and line[1] == "R":
                # get the id
                id = line[0][1:]
                TMdict[id] = {}
                # get the length of the line
                linelength = len(line)
                # startvalues
                sequence = ""
                TMnr = 1
                for colnr in range(2, linelength):
                    # write result to dictionary and reset variables
                    if (line[colnr] == "-" and len(sequence) == TMlengths[TMnr]):
                        TMdict[id][TMnr] = {}
                        TMdict[id][TMnr]["sequence"] = sequence
                        TMdict[id][TMnr]["sequence_length"] = len(sequence)
                        TMnr += 1
                        sequence = ""
                    else:
                        sequence += line[colnr]
                # write also the last sequence
                if colnr == linelength - 1:
                    TMdict[id][TMnr] = {}
                    TMdict[id][TMnr]["sequence"] = sequence
                    TMdict[id][TMnr]["sequence_length"] = len(sequence)
    print"human alignment read succesful..."
    return TMdict, TMlengths

# Calculates a score matrix based on the human alignment


def calculate_scores(human_alignment, amino_acids, TMlengths):
    # creat an empty score matrix
    score_matrix = {}
    # number of sequences in the alignment
    total_length = len(human_alignment.keys())
    # creat an empty vector which will store the maximum scores which will be
    # used for normalization
    maxscores = {}
    for TM in TMlengths:
        score_matrix[TM] = {}
        maxscores[TM] = []
        # cycle through the positions in the sequence position
        for pos in range(len(human_alignment[human_alignment.keys()[0]][TM]["sequence"])):
            score_matrix[TM][pos] = {}
            # create a empty list where all the aminoacids on that will be
            # stored temporarily
            posx = []
            # cycle all entries
            for id in human_alignment.keys():
                posx.append(human_alignment[id][TM]["sequence"][pos])
            # calculate the score
            scorepos = []
            for aa in amino_acids:
                score = math.pow(
                    float(len(filter(lambda pos: pos == aa, posx))) / total_length, 2)
                score_matrix[TM][pos][aa] = score
                scorepos.append(score)

            maxscores[TM].append(max(scorepos))
    print"score matrix created..."
    return score_matrix, maxscores
# scale scores in a way that the maximumscore for each TM is 1


def normalize_score_matrix(score_matrix, max_scores):
    for TM in max_scores.keys():
        total = sum(max_scores[TM])
        for pos in range(len(max_scores[TM])):
            max_scores[TM][pos] = float(max_scores[TM][pos]) / total
            for aa in score_matrix[TM][pos].keys():
                score_matrix[TM][pos][aa] = score_matrix[TM][pos][aa] / total
    return score_matrix


def build_hmm(human_alignment, filedir):
    for TM in range(1, 7 + 1):
        # build a hidden markov model
        if not os.path.exists("%s/TM_hmm_%i.dat" % (filedir, TM)):
            os.system("hmmbuild %s/TM_hmm_%i.dat %s/human_TM_%i.txt" %
                      (filedir, TM, filedir, TM))
    print"hmm models build..."

# for each TM a single file is made. This enables the calculation of the
# hidden markov models


def create_TM_files(human_alignment, filedir, TMlengths):
    nr = 0
    for TM in TMlengths:
        outfile = open("%s/human_TM_%i.txt" % (filedir, TM), "w")
        for id in human_alignment.keys():
            # give a unique number
            nr += 1
            outfile.write(">%i_%s_%i\n" % (nr, id, TM))
            outfile.write("%s\n" % (human_alignment[id][TM]["sequence"]))
        outfile.close()
    print"TM-sequence file created..."

# This loads the modules needed for this script, but does not work"


def load_modules():
    os.system(". /etc/profile.d/modules.sh;module load hmmer;")
    os.system(". /etc/profile.d/modules.sh;module load python;")
    print "modules loaded..."

# This function prerates all variables which are needed to perform a
# correct alignment


def pre_alignment(human_alignment, amino_acids, filedir, numberings):
    # make a directory to write output files to.
    if not os.path.exists("files"):
        os.mkdir("files")

    load_modules()
    # Read the human alignment in and store it in a dictionary
    human_alignment_dict = numberings.tm_sequences_of_alignment(human_alignment)
    TMlengths = numberings.tm_lengths()
    preferred_gaps = numberings.gpcrdb_gaps()

    # create a different fasta file for each TM of the human alignment
    create_TM_files(human_alignment_dict, filedir, TMlengths)
    # Calculate the score matrix that will be used for scoring the alignment
    # which will be made later.
    score_matrix, max_scores = calculate_scores(
        human_alignment_dict, amino_acids, TMlengths)
    # normalize the score matrix
    score_matrix = normalize_score_matrix(score_matrix, max_scores)
    # build hidden markov models
    build_hmm(human_alignment_dict, filedir)
    # return the variables which are needed in later functions
    return TMlengths, score_matrix, preferred_gaps

# a function that tries to find the first and last character from a
# characterlist in a string


def find_first_letter(letters, sequence):
    minlist = []
    maxlist = []
#    positions = map(lambda letter: filter(lambda position: sequence[position]==letter, range(len(sequence))),letters)
#    for list in positions:
#        if len(list)>0:
#            minlist.append(min(list))
#            maxlist.append(max(list))
#    minimum = min(minlist)
#    maximum = max(maxlist)

    mins = map(lambda letter: sequence.find("%s" % letter), letters)
    maxs = map(lambda letter: sequence.rfind("%s" % letter), letters)
    mins.sort()
    ap = mins.count(-1)
    mins = mins[ap:]
    maxs.sort()
    ap = maxs.count(-1)
    maxs = maxs[ap:]
    minimum = min(mins)
    maximum = max(maxs)

    return minimum, maximum

# This function aligns all sequences


def align_sequences(sequences, filedir, alignment_file):
    for TM in range(1, 7 + 1):
        # do the alignment
        if not os.path.exists("%s_%i.dat" % (alignment_file, TM)):
            os.system("hmmalign -o %s_%i.dat %s/TM_hmm_%i.dat %s" %
                      (alignment_file, TM, filedir, TM, sequences))
    print "alignment made..."

# this alignment reads the alignment file and extracts all information
# needed and stores this in a dictionary


def read_alignment(alignment, filedir, alignment_file, TMlengths):
    d = {}
    # create dictionary for storage of the alignment
    for TM in TMlengths:
        aligned_file = "%s_%i.dat" % (alignment_file, TM)
        content = open(aligned_file, "r").readlines()
        for line in content:
            if len(line) > 0 and line[0] <> "#":
                line = line.strip()
                line = line.split()
                if len(line) > 0 and line[0] <> "//":
                    if not alignment.has_key(d[line[0]]):
                        alignment[d[line[0]]] = {}
                        alignment[d[line[0]]][TM] = {}
                        alignment[d[line[0]]][TM]['sequence'] = line[1]
                    elif not alignment[d[line[0]]].has_key(TM):
                        alignment[d[line[0]]][TM] = {}
                        alignment[d[line[0]]][TM]['sequence'] = line[1]
                    else:
                        alignment[d[line[0]]][TM]['sequence'] += line[1]
            else:
                line = line.strip().split()
                if line[0] == "#=GS" and line[2] == "DE":
                    d[line[1]] = line[1]
                    for p in range(3, len(line)):
                        d[line[1]] += " " + line[p]
    print "alignment read succesful..."
    return alignment

# The alignment is shown as uppercase characters, this part needs to be
# cut out and stored. It also extract other useful information.


def cut_out_aligned_part(alignment, TMlengths):
    lowcase = string.ascii_lowercase
    upcase = string.ascii_uppercase
    upcase += "-"
    for id in alignment.keys():
        for TM in TMlengths:
            alignment[id][TM]["sequence"] = alignment[id][TM][
                "sequence"].replace(".", "").replace("-", "")

            minimumTM, maximumTM = find_first_letter(
                upcase, alignment[id][TM]['sequence'])
            helixseq = copy.deepcopy(alignment[id][TM]['sequence'][
                                     minimumTM:maximumTM + 1])
            helixseq = string.upper(helixseq)
            helixlength = len(helixseq)
            seqlength = len(alignment[id][TM]["sequence"])
            # store information in the dictionary
            alignment[id][TM]["helixlength"] = helixlength
            alignment[id][TM]["helixseq"] = helixseq
            alignment[id][TM]["seqlength"] = seqlength
            alignment[id][TM]["firstaa"] = int(minimumTM)
            alignment[id][TM]["lastaa"] = int(maximumTM + 1)
            alignment[id][TM]["cut_off"] = True
    return alignment

# score all sequence


def score_sequence(sequence, TM, score_matrix):
    score = 0.0
    for pos in range(len(score_matrix[TM])):
        score += score_matrix[TM][pos][sequence[pos]]
    return score


# some sequence contained gaps and therefor were to short. These lengths
# need to be fixed
def fix_wrong_lengths(alignment, TMlengths, score_matrix, preferred_gaps):
    for id in alignment.keys():
        for TM in TMlengths:
            tmr = alignment[id][TM]

            sequence = tmr["helixseq"]
            # if aligned sequence is too short, when it for example contained
            # gaps. enlarge the sequence by adding aminoacids at the front and
            # back of the sequence.
            if len(sequence) < TMlengths[TM]:
                # TODO retain gaps specified in gpcrdb_gap column in numbering scheme file
                diff = TMlengths[TM] - len(sequence)
                sequence = tmr["sequence"][
                    tmr["firstaa"] - diff:tmr["lastaa"] + diff]
                sequence = string.upper(sequence)
                tmr["helixseq"] = sequence
                tmr["firstaa"] = tmr["firstaa"] - diff
                tmr["lastaa"] = tmr["lastaa"] + diff
                # now helixseq is too long, below we truncate it again

            # When sequence is too large, choose the best scoring part with the
            # correct length
            if len(sequence) > TMlengths[TM]:
                result = []
                for start in range(0, len(sequence) - TMlengths[TM]):
                    possib = sequence[start:start + TMlengths[TM]]
                    result.append([score_sequence(possib, TM, score_matrix),
                                   possib,
                                   start, 0])

                    # generate sequences with gaps added to them
                    gaps = preferred_gaps.get(TM, set())
                    for r in range(1, len(gaps) + 1):
                        for combi in itertools.combinations(gaps, r):
                            # combi will have every combination of gaps list
                            # eg gaps = [5,10] then combi will be [5],[10],[5,10]
                            possib = list(sequence[start:start + TMlengths[TM] - r])
                            for gappos in combi:
                                # gappos is sequence pos so starts at 1 while possib list starts at 0
                                possib.insert(gappos - 1, '-')
                            possib = ''.join(possib)
                            result.append([score_sequence(possib, TM, score_matrix),
                                          possib,
                                          start, r])

                # use seq with highest score
                highest = sorted(result, key=lambda d: d[0]).pop()
                tmr["helixseq"] = highest[1]
                tmr["firstaa"] = tmr["firstaa"] + highest[2]
                tmr["lastaa"] = tmr["firstaa"] + TMlengths[TM] - highest[3]

    print "gaps fixed..."
    return alignment

# each sequence needs to be scored.


def score_alignment(alignment, score_matrix, TMlengths):
    for id in alignment.keys():
        for TM in TMlengths:
            # calculate the score of each alignment
            if len(alignment[id][TM]["helixseq"]) == len(score_matrix[TM]):
                alignment[id][TM]["score"] = score_sequence(
                    alignment[id][TM]["helixseq"], TM, score_matrix)
            # for these sequences it was not possible to fill the gips, mostly
            # because of an incomplete sequence
            else:
                alignment[id][TM]["score"] = 0.0
    print "alignment scored..."
    return alignment


def calculate_cut_off(alignment, TMlengths):
    TM_cut_off = {}
    for TM in TMlengths:
        all_scores = []
        for id in alignment.keys():
            all_scores.append(float(alignment[id][TM]["score"]))
        all_scores.sort()
        TM_cut_off[TM] = all_scores[int(0.95 * (len(all_scores) + 1)) - 1]
    print TM_cut_off
    return TM_cut_off


def filter_alignment(alignment, filedir, cut_off, severe_sequences, allowed_low_TMs, TMlengths):
    #
    for TM in TMlengths:
        for id in alignment.keys():
            if alignment[id][TM]["score"] < cut_off[TM]:
                alignment[id][TM]["cut_off"] = False
    total_file = open('%s/total_alignment.dat' % (filedir), 'w')
    error_file = open('%s/error_total_alignment.dat' % (filedir), 'w')
    nr_file = open('%s/TM_numbering.dat' % (filedir), 'w')

    high = 0
    low = 0
    for id in alignment.keys():
        low_TMs = 0
        for TM in TMlengths:
            if not alignment[id][TM]["cut_off"]:
                low_TMs += 1
        # write id's + sequences to file which have less than the desired TMs
        # aligned correct

        if low_TMs <= allowed_low_TMs and \
                        id not in severe_sequences and \
                        alignment[id]["length"] and \
                        alignment[id][1]["score"] > (cut_off[1] / 2.0) and \
                        alignment[id][7]["score"] > (cut_off[7] / 2.0):
            high += 1
            total_file.write('>%s\n' % id)
            for TM in TMlengths:
                if len(alignment[id][TM]["helixseq"]) > 30:
                    print id, TM
                if TM < 7:
                    total_file.write(alignment[id][TM]["helixseq"] + '-')
                else:
                    total_file.write(alignment[id][TM]["helixseq"] + '\n')
                nr_file.write(">%s\t%s\t%s\t%s\t%s\n" % (id, TM, alignment[id][TM][
                              "firstaa"], alignment[id][TM]["lastaa"], alignment[id][TM]["score"]))
        else:
            low += 1
            eTM = ""
            for TM in TMlengths:
                eTM += "%s:%s.%s " % (alignment[id][TM]["score"], alignment[id][
                                      TM]["firstaa"], alignment[id][TM]["lastaa"])
            error_file.write('>%s (%s)\n' % (id, eTM))
            for TM in TMlengths:
                if TM < 7:
                    error_file.write(alignment[id][TM]["helixseq"] + '_')
                else:
                    error_file.write(alignment[id][TM]["helixseq"] + '\n')
    print high
    print low


def sort_matrix(matrix):
    # sort a matrix based on a column
    matrix.sort(lambda x, y: cmp(x[1], y[1]))
    return matrix


def order_check(TMlist):
    # number of incorrect involved in incorrect oridering of the helices
    incorrect = 0
    firstaa = []
    falsehelix = -1
    for nr in range(len(TMlist) - 1):
        firstaa.append([nr, TMlist[nr][0]])
        # find out if the order is correct
        if TMlist[nr][1] >= TMlist[nr + 1][0]:
            incorrect += 1
            falsehelix = nr
    # do also append the last entry
    firstaa.append([nr + 1, TMlist[nr + 1][0]])
    return incorrect, firstaa, falsehelix


def overlap(TMlist, TMlengths, IDdict, redo):
    if redo[0] == 0:
        sequencelimits = [[0, TMlist[redo[0] + 2]
                           [0] - TMlengths[(int(redo[0]) + 1)] - 1]]
        sequencelimits.append(
            [0 + TMlengths[(int(redo[0]) + 1)] + 1, TMlist[redo[0] + 2][0] - 1])
    elif redo[1] == 6:
        sequencelimits = [[TMlist[redo[0] - 1][1] + 1, IDdict[1]
                           ["seqlength"] - TMlengths[(int(redo[0]) + 1)] - 1]]
        sequencelimits.append(
            [TMlist[redo[0] - 1][1] + TMlengths[(int(redo[0]) + 1)] + 1, IDdict[1]["seqlength"]])
    else:
        sequencelimits = [[TMlist[redo[0] - 1][1] + 1,
                           TMlist[redo[0] + 2][0] - TMlengths[(redo[0] + 1)] - 1]]
        sequencelimits.append(
            [TMlist[redo[0] - 1][1] + TMlengths[(int(redo[0]) + 1)] + 1, TMlist[redo[0] + 2][0] - 1])
    return sequencelimits


def swap(TMlist, TMlengths, IDdict, redo):
    if redo[0] == 0:
        sequencelimits = [[0, TMlist[redo[0] + 2]
                           [0] - TMlengths[(int(redo[0]) + 1)] - 1]]
        sequencelimits.append(
            [0 + TMlengths[(int(redo[0]) + 1)] + 1, TMlist[redo[0] + 2][0] - 1])
    elif redo[1] == 6:
        sequencelimits = [[TMlist[redo[0] - 1][1] + 1, IDdict[1]
                           ["seqlength"] - TMlengths[(redo[0] + 1)] - 1]]
        sequencelimits.append(
            [TMlist[redo[0] - 1][1] + TMlengths[(int(redo[0]) + 1)] + 1, IDdict[1]["seqlength"]])
    else:
        sequencelimits = [[TMlist[redo[0] - 1][1] + 1,
                           TMlist[redo[0] + 2][0] - TMlengths[(redo[0] + 1)] - 1]]
        sequencelimits.append(
            [TMlist[redo[0] - 1][1] + TMlengths[(int(redo[0]) + 1)] + 1, TMlist[redo[0] + 2][0] - 1])
    return sequencelimits


def solve_minor_prob(incorrect, firstaa, redo, sequencelimits, errorcode, falsehelix, IDdict, TMlist, TMlengths):
    # There can three types of errors:
    # 1. 1 helix is placed incorrect --> replace this helix this is the helix with the highest difference. One difference is higher than 1 and the rest is 1.
    # 2. 2 helices overlap each other --> try to align with shorter sequence
    # 3. 2 helices are swapped -->try to align with shorter sequence
    checklist = []
    firstaa = sort_matrix(firstaa)
    # calculate the difference list for this list
    pos = 0
    for entry in firstaa:
        # The checklist contains the difference between the annotated TM and
        # the TM based on the sequence.
        checklist.append([entry[0], abs(entry[0] - pos)])
        pos += 1
    checklist = sort_matrix(checklist)
    checklist.reverse()
    #???
    possible = False
    if checklist[0][1] == 0:
        errorcode = 1
        # redo both helices with shortent sequences
        redo = [falsehelix, falsehelix + 1]
        # calculate both borders
        sequencelimits = overlap(TMlist, TMlengths, IDdict, redo)
    elif checklist[0][1] == 1 and checklist[1][1] == 1 and checklist[2][1] == 0:
        errorcode = 2
        redo = [checklist[0][0], checklist[1][0]]
        sequencelimits = swap(TMlist, TMlengths, IDdict, redo)
    elif checklist[0][1] > 1 and checklist[1][1] < 2:
        posredo = checklist[0][0]

        for entry in firstaa:
            if entry[0] == checklist[0][0]:
                if entry[0] == 0:
                    entry[1] = 1
                else:
                    # print "test: %s"%TMlist[checklist[0][0]-1]
                    entry[1] = TMlist[checklist[0][0] - 1][1] + 1
        firstaa = sort_matrix(firstaa)
        # calculate the difference list for this list
        pos = 0
#        print "checklist old: %s"%checklist
        checklist = []
#        print "new firstaa %s"%firstaa
        for entry in firstaa:
            checklist.append([entry[0], abs(entry[0] - pos)])
            pos += 1
        checklist = sort_matrix(checklist)
        checklist.reverse()

        if checklist[0][1] == 0:
            errorcode = 3
#            print "redoing one helix is succesfull"
            redo = [posredo]
            if redo[0] == 0:
                sequencelimits = [[0, TMlist[redo[0] + 1][0] - 1]]
            elif redo[0] == 6:
                sequencelimits = [
                    [TMlist[redo[0] - 1][1] + 1, IDdict[1]["seqlength"]]]
            else:
                sequencelimits = [
                    [TMlist[redo[0] - 1][1] + 1, TMlist[redo[0] + 1][0] - 1]]
        else:
            errorcode = 4
    else:
        errorcode = 5
    return redo, sequencelimits, errorcode


def TMorder(TMlist, IDdict, TMlengths, id, count):
    redo = []
    sequencelimits = []
    # the code reflect the severenous of the problem
    errorcode = 0

    # Count the helices that are in the wrong order
    incorrect, firstaa, falsehelix = order_check(TMlist)
    # When there is only 1 helix aligned incorrect
    if incorrect == 1:
        redo, sequencelimits, errorcode = solve_minor_prob(
            incorrect, firstaa, redo, sequencelimits, errorcode, falsehelix, IDdict, TMlist, TMlengths)
    if incorrect > 1:
        errorcode = 5

    return incorrect, redo, sequencelimits, errorcode, count


def check_TMorder(alignment, TMlengths, filedir):
    count = 0
    severe_sequences = []
    for id in alignment.keys():
        alignment[id]["TMorder"] = True
        TMlist = []
        for TM in TMlengths:
            TMlist.append([alignment[id][TM]["firstaa"],
                           alignment[id][TM]["lastaa"]])

#        print "TMlist: %s"%TMlist

        # checks the TM-order and detects which error has been made
        incorrect, redo, sequencelimits, errorcode, count = TMorder(
            TMlist, alignment[id], TMlengths, id, count)

        # errors 4 and 5 are not solvable
        if errorcode == 5 or errorcode == 4:
            severe_sequences.append(id)
            alignment[id]["TMorder"] = False
        # sequences which has been aligned incorrectly should be aligned again.
        # Now with a smaller sequence where the TM should lie.
        elif incorrect > 0:
            #            print "TMlist: %s"%TMlist
            #            print "Redo: %s"%redo
            #            print "Sequencelimits: %s"%sequencelimits
            alignment[id]["redo"] = redo
            alignment[id]["sequencelimits"] = sequencelimits

            for nr in range(len(redo)):
                if sequencelimits[nr][1] - sequencelimits[nr][0] > TMlengths[(redo[nr] + 1)]:
                    sequence = alignment[id][1]["sequence"][
                        sequencelimits[nr][0]:sequencelimits[nr][1]]
                    outfile = open("%s/redo_%i.dat" %
                                   (filedir, redo[nr] + 1), "a")
                    outfile.write(">%s\n" % (id))
                    outfile.write("%s\n" % (sequence))
                    outfile.close()
                else:
                    severe_sequences.append(id)
                    alignment[id]["TMorder"] = False

    # redo the alignment for the sequences which contained errors
    for TM in TMlengths:
        fastafile = "%s/redo_%i.dat" % (filedir, TM)
        # only do hmmalign if input file exists
        if os.access(fastafile, os.R_OK):
            os.system("hmmalign -o %s/TM_aligned_redo_%i.dat %s/TM_hmm_%s.dat %s" %
                      (filedir, TM, filedir, TM, fastafile))

    # read the re-aligned sequences
    redo_file = "%s/TM_aligned_redo" % filedir
    for TM in TMlengths:
        aligned_file = "%s_%i.dat" % (redo_file, TM)
        # only do hmmalign if output file exists
        if os.access(aligned_file, os.R_OK):
            d = {}
            content = open(aligned_file, "r").readlines()
            for line in content:
                if len(line) > 0 and line[0] != "#":
                    line = line.strip()
                    line = line.split()
                    if len(line) > 0 and line[0] != "//":
                        if alignment[d[line[0]]][TM].has_key("redo_sequence"):
                            alignment[d[line[0]]][TM][
                                'redo_sequence'] += line[1]
                        else:
                            alignment[d[line[0]]][TM][
                                'redo_sequence'] = line[1]
                else:
                    line = line.strip().split()
                    if line[0] == "#=GS" and line[2] == "DE":
                        d[line[1]] = line[1]
                        for p in range(3, len(line)):
                            d[line[1]] += " " + line[p]

    # Cut the realigned sequence and replace the old aligned sequence
    lowcase = string.ascii_lowercase
    upcase = string.ascii_uppercase
    upcase += "-"
    for id in alignment.keys():
        for TM in TMlengths:
            if alignment[id][TM].has_key('redo_sequence'):
                alignment[id][TM]['redo_sequence'] = alignment[id][TM][
                    'redo_sequence'].replace(".", "").replace("-", "")
                minimumTM, maximumTM = find_first_letter(
                    upcase, alignment[id][TM]['redo_sequence'])
                # get the limits of the redone alignment
                alignment[id][TM]["firstaa"] = int(alignment[id]["sequencelimits"][
                                                   alignment[id]["redo"].index(TM - 1)][0] + minimumTM)
                alignment[id][TM]["lastaa"] = int(alignment[id]["sequencelimits"][
                                                  alignment[id]["redo"].index(TM - 1)][0] + maximumTM) + 1
                helixseq = copy.deepcopy(alignment[id][TM]['redo_sequence'][
                                         minimumTM:maximumTM + 1])
                helixseq = string.upper(helixseq)
                helixlength = len(helixseq)
                alignment[id][TM]["helixseq"] = helixseq

        TMlist = []
        if severe_sequences.count(id) == 0:
            for TM in range(1, 7 + 1):
                TMlist.append([alignment[id][TM]["firstaa"],
                               alignment[id][TM]["lastaa"]])

        #    print "TMlist: %s"%TMlist

        if len(TMlist) > 0:
            # REcheck the TM-order and detects which error has been made
            incorrect, redo, sequencelimits, errorcode, count = TMorder(
                TMlist, alignment[id], TMlengths, id, count)
            # remove if the ordering is still incorrect !!
            if errorcode > 0:
                severe_sequences.append(id)
                alignment[id]["TMorder"] = False

    return alignment, severe_sequences


def check_lengths(alignment, TMlengths):
    for id in alignment.keys():
        alignment[id]["length"] = True
        for TM in TMlengths:
            if len(alignment[id][TM]["helixseq"]) != TMlengths[TM]:
                alignment[id]["length"] = False

    return alignment


def alignment(sequences, TMlengths, score_matrix, cut_off, calc_treshold, allowed_low_TMs, filedir, preferred_gaps):
    alignment = {}
    alignment_file = "%s/TM_aligned" % filedir
    # align the sequences
    align_sequences(sequences, filedir, alignment_file)
    # read in the alignment file produced by hmmer
    alignment = read_alignment(alignment, filedir, alignment_file, TMlengths)
    # the aligned part of the sequence need to be cut out and gaps need to be
    # fixed.
    alignment = cut_out_aligned_part(alignment, TMlengths)
    # fill gaps in the aligned part
    alignment = fix_wrong_lengths(alignment, TMlengths, score_matrix, preferred_gaps)
    # score the alignments
    alignment = score_alignment(alignment, score_matrix, TMlengths)
    # Check the TM-order and correct it if necessary
    alignment, severe_sequences = check_TMorder(alignment, TMlengths, filedir)

    if calc_treshold:
        cut_off = calculate_cut_off(alignment, TMlengths)

    ###############
    alignment = check_lengths(alignment, TMlengths)

    filter_alignment(alignment, filedir, cut_off,
                     severe_sequences, allowed_low_TMs, TMlengths)


def clean_map(filedir):
    for TM in range(1, 1 + 7):
        if os.path.exists("%s/redo_%i.dat" % (filedir, TM)):
            os.remove("%s/redo_%i.dat" % (filedir, TM))
        if os.path.exists("%s/TM_aligned_%i.dat" % (filedir, TM)):
            os.remove("%s/TM_aligned_%i.dat" % (filedir, TM))
        if os.path.exists("%s/TM_aligned_redo_%i.dat" % (filedir, TM)):
            os.remove("%s/TM_aligned_redo_%i.dat" % (filedir, TM))
        if os.path.exists("%s/redo_%i.dat" % (filedir, TM)):
            os.remove("%s/redo_%i.dat" % (filedir, TM))
#  ======================================================================
#                          S U B S C R I P T   1
#  ======================================================================cd fi


def align(variables, sequences):
    # set workdirectory to directory of script
    workdir = os.getcwd()

    filedir = "%s/files" % (workdir)

    # list of all amino-acids (one letter code)
    amino_acids = ["A", "R", "N", "D", "C", "Q", "E", "G", "H", "I", "L",
                   "K", "M", "F", "P", "S", "T", "W", "Y", "V", "X", "B", "U", "Z", "-"]

    human_alignment, cut_off, allowed_low_TMs, calc_treshold, numberings = get_variables(
        variables)

    TMlengths, score_matrix, preferred_gaps = pre_alignment(
        human_alignment, amino_acids, filedir, numberings)

    alignment(sequences, TMlengths, score_matrix, cut_off,
              calc_treshold, allowed_low_TMs, filedir, preferred_gaps)

    clean_map(filedir)

#  ======================================================================
#                          M A I N   S C R I P T
#  ======================================================================


def main():
    # PRINT WELCOME SCREEN
    txt_header()

    # get the name of the computer on which we are running
    hostname = socket.gethostname()

    # CHECK WHICH OF THE SUBSCRIPTS SHOULD BE RUN
    args = len(sys.argv)
    if (args <= 1):
        # PRINT HELP
        # ==========
        print
        print "**********************************************************"
        print "* Use the following command line options:                *"
        print "**********************************************************"
        print
        print "HELP FUNCTIONS:"
        print "==============="
        print "alignment.py"
        print "         -> Prints this help screen"
        print
        print "CLUSTER FUNCTIONS:"
        print "=================="
        print "alignment.py -run"
        print "          -> Do the alignment"
        print "=================="
        print "alignment.py -dev"
        print "         -> Just a testing function."
        print
        raise SystemExit

    else:

        # RUN FUNCTION
        # ================
        name, nops = "run", 2
        if (sys.argv[1] == "-%s" % name):
            if (args != 2 + nops):
                err_error(
                    "Number of command line parameters does not match -%s option." % name)
            # RUN THE DEVELOPMENT CODE
#      align(sys.argv[2],sys.argv[3])
            variables = sys.argv[2]
            sequences = sys.argv[3]
            align(variables, sequences)
        # END
        # ============

        # DEV FUNCTION
        # ============
        name, nops = "dev", 0
        if (sys.argv[1] == "-%s" % name):
            if (args != 2 + nops):
                err_error(
                    "Number of command line parameters does not match -%s option." % name)
            # RUN THE DEVELOPMENT CODE
#      dev(sys.argv[2],sys.argv[3])
            dev()
        # END
        # ============

# EXECUTE SCRIPT
if __name__ == "__main__":
    main()
