"""H35's table of a left-to-right model: which reading a place takes, the minus strand turned back, the mean of the
two strands. The model is a stand-in whose answer is worked out by hand."""
import numpy as np
import pytest

from mitoage_pipelines.mitoport import evo_table as et

SEQ = "ACGTTGCAAC"          # 10 bases; its reverse complement is GTTGCAACGT


def counting(seq: str) -> np.ndarray:
    """Stand-in: at place i the log-probability of each base is log of (1 + how often it occurs among the bases
    before i), so the answer depends on what precedes the place and on nothing after it."""
    out = np.zeros((len(seq), 4))
    for i in range(len(seq)):
        c = np.array([1 + seq[:i].count(b) for b in "ACGT"], float)
        out[i] = np.log(c)
    return out


class TestOneStrand:
    def test_each_place_takes_the_reading_with_more_bases_before_it(self):
        got = et.one_strand(SEQ, counting)
        turned = SEQ[5:] + SEQ[:5]
        for p in range(10):
            if p >= 5:                                   # read from the start: p bases before
                want = counting(SEQ)[p]
            else:                                        # read from the middle: p + 5 bases before
                want = counting(turned)[p + 5]
            assert np.array_equal(got[p], want), p
        # place 0 from the start would have nothing before it; here it has the genome's second half
        assert np.array_equal(got[0], np.log([1 + SEQ[5:].count(b) for b in "ACGT"]))

    def test_a_wrong_shape_is_refused(self):
        with pytest.raises(ValueError, match="returned"):
            et.one_strand(SEQ, lambda s: np.zeros((len(s) - 1, 4)))


class TestStrands:
    def test_reverse_complement(self):
        assert et.reverse_complement(SEQ) == "GTTGCAACGT" and et.reverse_complement("AACG") == "CGTT"

    def test_the_minus_strand_is_turned_back_to_the_plus_strands_places_and_bases(self):
        minus = np.arange(12, dtype=float).reshape(3, 4)          # rows: places of the reverse complement; columns A C G T
        back = et.turned_back(minus)
        # the plus strand's place 0 is the reverse complement's place 2; its A is that strand's T
        assert back[0, 0] == minus[2, 3] and back[0, 1] == minus[2, 2] and back[2, 3] == minus[0, 0]

    def test_tables(self):
        t = et.tables(SEQ, counting)
        plus = et.normalised(et.one_strand(SEQ, counting))
        minus = et.normalised(et.one_strand("GTTGCAACGT", counting)[::-1, ::-1])
        assert np.allclose(t["plus"], plus) and np.allclose(t["minus"], minus)
        assert not np.allclose(t["plus"], t["minus"])
        assert np.allclose(np.exp(t["both"]).sum(axis=1), 1)
        raw = (plus + minus) / 2
        assert np.allclose(t["both"] - raw, (t["both"] - raw)[:, :1])          # the mean, moved by one constant per place
        assert np.allclose(np.exp(t["plus"]).sum(axis=1), 1)

    def test_an_unknown_base_is_not_scored_and_case_does_not_matter(self):
        t = et.tables("ACGTNGCAAC", counting)
        assert np.isnan(t["both"][4]).all() and not np.isnan(t["both"][3]).any() and np.isnan(t["plus"][4]).all()
        assert np.array_equal(et.tables(SEQ.lower(), counting)["both"], et.tables(SEQ, counting)["both"])
