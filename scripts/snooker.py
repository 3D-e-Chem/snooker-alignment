
from collections import OrderedDict
import csv
from Bio import AlignIO

ALIGNMENT_POSITION = 'aligmnent_position'
TM = 'tm'
GPCRDB_GAP = 'gpcrdb_gap'

class Numberings(object):
    name = ''
    rows = []

    def __init__(self, name, rows):
        self.name = name
        self.rows = rows

    @classmethod
    def from_filename(cls, fn):
        return cls.from_file(open(fn))

    @classmethod
    def from_file(cls, numberformatfile):
        rows = []
        name = ''
        with numberformatfile:
            nfreader = csv.reader(numberformatfile, delimiter=",")
            header = next(nfreader)
            name = header[0]
            header[0] = ALIGNMENT_POSITION
            for line in nfreader:
                row = {header[i]: v for i, v in enumerate(line) if v != ''}
                rows.append(row)
        return Numberings(name, rows)

    def tm_lengths(self):
        lengths = OrderedDict()
        for row in self.rows:
            if ALIGNMENT_POSITION in row and TM in row:
                tm_id = int(row[TM])
                if tm_id not in lengths:
                    lengths[tm_id] = 0
                lengths[tm_id] += 1
        return lengths

    def tm_starts(self):
        starts = {}
        for row in self.rows:
            if ALIGNMENT_POSITION in row and TM in row:
                tm_id = int(row[TM])
                if tm_id not in starts:
                    starts[tm_id] = int(row[ALIGNMENT_POSITION])
        return starts

    def tm_sequences(self, sequence):
        seqs = OrderedDict()
        for row in self.rows:
            if ALIGNMENT_POSITION in row and TM in row:
                tm_id = int(row[TM])
                aln_pos = int(row[ALIGNMENT_POSITION])
                if tm_id not in seqs:
                    seqs[tm_id] = {'sequence': '', 'sequence_length': 0}
                seqs[tm_id]['sequence'] += sequence[aln_pos - 1]
                seqs[tm_id]['sequence_length'] += 1

        return seqs

    def tm_sequences_of_alignment(self, fn, format='fasta'):
        records = OrderedDict()
        alignment = AlignIO.read(open(fn), format)
        for row in alignment:
            records[row.description] = self.tm_sequences(row.seq)
        return records

    def lookup(self, key, value):
        key2value = OrderedDict()
        for row in self.rows:
            if key in row and value in row:
                keyval = row[key]
                if key is ALIGNMENT_POSITION:
                    keyval = int(keyval)
                val = row[value]
                if value is ALIGNMENT_POSITION:
                    val = int(val)
                key2value[keyval] = val
        return key2value

    def gpcrdb_gaps_alignment_positions(self):
        gaps = {}
        for row in self.rows:
            if ALIGNMENT_POSITION in row and TM in row and GPCRDB_GAP in row:
                tm_id = int(row[TM])
                ap = int(row[ALIGNMENT_POSITION])
                if tm_id not in gaps:
                    gaps[tm_id] = set()
                gaps[tm_id].add(ap)
        return gaps

    def gpcrdb_gaps(self):
        gaps = self.gpcrdb_gaps_alignment_positions()
        starts = self.tm_starts()
        ngaps = {}
        for tm in gaps:
            tmgaps = set()
            for gap in gaps[tm]:
                tmgaps.add(gap + 1 - starts[tm])
            ngaps[tm] = tmgaps
        return ngaps