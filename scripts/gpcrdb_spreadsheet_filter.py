#!/usr/bin/env python

"""

Filter gcprdb csv alignment by removing unwanted positions.
Generates a fasta file.

Usage:

gpcrdb_spreadsheet_filter.py -n gpcrdb_gapped_tm_numbering.csv \
-f gpcr-a.members.json\
< gpcrdb_human_swissprot_tm_aligment.csv > gpcrdb_human_swissprot2.fa

The gcprdb csv alignment can be obtained by
1. Goto http://gpcrdb.org/alignment/targetselection
* Filters: human + swissprot
* Receptors: class a rhodopsin like
* Sequence segments: all tms
2. Export CSV

The family file can be fetched from http://gpcrdb.org/services/proteinfamily/proteins/001/
Where 001 equals `Class A (Rhodopsin)` family.

"""

import argparse
import collections
import csv
import json
import re
import sys
import snooker


def gpcrdbalignment2snookeralignment_mapping(gpositions,
                                             gprcdb2snookeralignment,
                                             ):
    c2l = {}
    for cpos, gpos in enumerate(gpositions):
        if gpos in gprcdb2snookeralignment:
            c2l[cpos] = gprcdb2snookeralignment[gpos]
    return c2l


def gpcrdbalignment2filteredsnookeralignment(row, gpcrdb2snooker, outlen):
    sid = row.pop(0)
    aas = row
    aa_numbered = ['-' for i in range(outlen)]
    for cpos, apos in c2l.items():
        aa_numbered[apos - 1] = aas[cpos]
    seq = ''.join(aa_numbered)
    return (sid, seq,)


def gpcrdb_name2id_mapping(family_file, tax_mapping):
    """Creates a mapping from the sequence id to its swissprot id

    Alignment file contains '[Human] bla receptor' as sequence id,
    these can be mapped using the proteins returned by the GPCRDB API.

    Species missing in tax_mapping will be skipped.
    """
    name2id = {}
    if family_file is None:
        return name2id

    # strip html tags using regexp
    TAG_RE = re.compile(r'<[^>]+>')

    family_members = json.load(family_file)
    for member in family_members:
        organism = tax_mapping.get(member['species'])
        name = TAG_RE.sub('', member['name']).replace('&', '&amp;')
        if organism is None:
            continue
        key = '[{0}] {1}'.format(organism, name)
        entry_name = member['entry_name'].upper()
        name2id[key] = entry_name
    return name2id

parser = argparse.ArgumentParser()
parser.add_argument('-n', '--numbering', type=argparse.FileType('r'))
parser.add_argument('-f', '--family', type=argparse.FileType('r'))
parser.add_argument('infile', nargs='?', type=argparse.FileType('r'), default=sys.stdin)
args = parser.parse_args()

numberings = snooker.Numberings.from_file(args.numbering)
g2a = numberings.lookup('gpcrdb_alignment', snooker.ALIGNMENT_POSITION)

outlen = max(g2a.values())

csvinreader = csv.reader(args.infile, quotechar="'")
next(csvinreader)  # skip line, not interested in transmembrane identifiers
gpositions = next(csvinreader)
gpositions.pop(0)  # gpcrdb export starts with weird utf char, remove it

c2l = gpcrdbalignment2snookeralignment_mapping(gpositions, g2a)

tax_mapping = {
    'Homo sapiens': 'Human'
}
name2id = gpcrdb_name2id_mapping(args.family, tax_mapping)

for row in csvinreader:
    if len(row) == 0:
        # Skip empty lines
        continue
    if row[1].isdigit():
        # Skip residue properties section
        continue
    if row[0] == 'CONSENSUS':
        continue
    if len(collections.Counter(row[1:])) == 1:
        # Only one kind of aa, probably all gaps
        continue
    (sid, seq) = gpcrdbalignment2filteredsnookeralignment(row, c2l, outlen)
    print('>{}\n{}'.format(name2id.get(sid, sid), seq))
