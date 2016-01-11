from unittest import TestCase
from snooker import Numberings, ALIGNMENT_POSITION, GPCRDB_GAP, TM


class NumberingsTestCase(TestCase):
    def setUp(self):
        rows = [{
            TM: 1
        }, {
            ALIGNMENT_POSITION: 1,
            TM: 1
        }, {
            ALIGNMENT_POSITION: 2,
            GPCRDB_GAP: 1,
            TM: 1
        }, {
            ALIGNMENT_POSITION: 3,
            TM: 1
        }, {
            ALIGNMENT_POSITION: 4,
            TM: 1
        }, {
            ALIGNMENT_POSITION: 5,
            GPCRDB_GAP: 1,
            TM: 1
        }, {
            ALIGNMENT_POSITION: 6,
            TM: 1
        }, {
            ALIGNMENT_POSITION: 7,
            TM: 1
        }, {
            ALIGNMENT_POSITION: 8
        }, {
            ALIGNMENT_POSITION: 9,
            TM: 2
        }, {
            ALIGNMENT_POSITION: 10,
            TM: 2,
            GPCRDB_GAP: 1
        }]
        self.numbering = Numberings('gpcrdb_gapped_tm', rows)

    def test_tm_lengths(self):
        lengths = self.numbering.tm_lengths()
        expected = {1: 7, 2: 2}
        self.assertDictEqual(lengths, expected)

    def test_tm_starts(self):
        starts = self.numbering.tm_starts()
        expected = {1: 1, 2: 9}
        self.assertDictEqual(starts, expected)

    def test_gpcrdb_gaps(self):
        gaps = self.numbering.gpcrdb_gaps()
        expected = {
                    1: set([2, 5]),
                    2: set([2])
                    }
        self.assertDictEqual(gaps, expected)