#!/usr/bin/env python

"""
Calculate entropies of each leaf on each branch node of a tree for each column

Usage:

entropy.py -a ali.s9t100.fa -n gpcrdb_gapped_tm_numbering.csv -t ali.s9t100.ph > ali.s9t100.entropies

"""

import argparse
import collections
import logging
import math

from Bio import Phylo, AlignIO

import snooker


def calculate_entropies(tree_file, alignment_file, numbering_file,
                        min_node_size, max_node_size, number_format):
    numberings = snooker.Numberings.from_file(numbering_file)
    ali2gpcrdb = numberings.lookup(snooker.ALIGNMENT_POSITION, number_format)

    alignment = AlignIO.read(alignment_file, 'fasta')
    id2seq = {row.id: row.seq for row in alignment}

    tree = Phylo.read(tree_file, 'newick')
    all_leafs = set([leaf.name for leaf in tree.get_terminals()])

    # for each column determine the aa distribution
    all_counters = {}
    for col in ali2gpcrdb:
        all_counters[col] = collections.Counter([seq[col - 1] for seq in id2seq.values()])

    print('{},{},{},{},{},{},{},{}'.format('node_id',
                                           'alignment_pos',
                                           number_format,
                                           'entropy_inside',
                                           'entropy_outside',
                                           'score',
                                           'variability_inside',
                                           'variability_outside',
                                           ))

    for node_id, node in enumerate(tree.get_nonterminals()):
        leafs_of_node = set([leaf.name for leaf in node.get_terminals()])
        if not (min_node_size <= len(leafs_of_node) <= max_node_size):
            msg = '{} has {} leafs, skipping'.format(node, len(leafs_of_node))
            logging.info(msg)
            continue

        leafs_outside_node = all_leafs - leafs_of_node

        seqs_inside = [id2seq[v] for v in leafs_of_node]
        nr_inside = float(len(leafs_of_node))
        nr_outside = float(len(leafs_outside_node))

        # loop over columns
        for col in ali2gpcrdb:
            aa_inside = collections.Counter([seq[col - 1] for seq in seqs_inside])
            f_i_inside = 0
            for count in aa_inside.values():
                f_i_inside += count / nr_inside * math.log(count / nr_inside)
            entropy_inside = -1 * f_i_inside
            variability_inside = len(aa_inside)

            aa_outside = all_counters[col] - aa_inside
            f_i_outside = 0
            for aa, count in aa_outside.items():
                f_i_outside += count / nr_outside * math.log(count / nr_outside)
            entropy_outside = -1 * f_i_outside
            variability_outside = len(aa_outside)

            distinct_aa = 21  # all amino acids and gap (-)
            score = math.sqrt(pow(abs(math.log(1.0 / distinct_aa)) - entropy_outside, 2)
                              + pow(entropy_inside, 2))

            print('{},{},{},{},{},{},{},{}'.format(node_id,
                                                   col,
                                                   ali2gpcrdb[col],
                                                   entropy_inside,
                                                   entropy_outside,
                                                   score,
                                                   variability_inside,
                                                   variability_outside,
                                                   ))

parser = argparse.ArgumentParser(description='Calculate entropies of each leaf on each branch node of a tree for each column')
parser.add_argument('-a', '--alignment', type=argparse.FileType('r'), required=True, help='Multiple sequence alignment (fasta format)')
parser.add_argument('-n', '--numbering', type=argparse.FileType('r'), required=True, help='Numbering file, translate sequence alignment position into generic numbering scheme')
parser.add_argument('-t', '--tree', type=argparse.FileType('r'), required=True, help='Tree of multiple sequence alignment (newick format)')
parser.add_argument('--min_node_size', type=int, default=20, help='Calculate entropies for nodes with a minimum number of leafs')
parser.add_argument('--max_node_size', type=int, default=20, help='Calculate entropies for nodes with a maximum number of leafs')
parser.add_argument('--number_format', default='gpcrdb_alignment', help='Column from numbering file to include in output')

args = parser.parse_args()

calculate_entropies(args.tree, args.alignment, args.numbering,
                    args.min_node_size, args.max_node_size, args.number_format)
