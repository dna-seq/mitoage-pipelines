"""H31's likelihood table, with a stand-in for the model: which places are hidden together, that every place is
hidden exactly once in each reading, the padding, and from which reading a place takes its value."""
import numpy as np
import pytest

from mitoage_pipelines.mitoport import kron_table as kt


def test_encode_and_padding():
    assert kt.encode("ACGTNacgtR").tolist() == [0, 1, 2, 3, 4, 0, 1, 2, 3, 4]
    assert kt.padded(np.array([0, 1, 2]), 4).tolist() == [0, 1, 2, 4]
    assert kt.padded(np.array([0, 1, 2, 3]), 4).tolist() == [0, 1, 2, 3]


def test_groups_cover_every_place_once_and_are_a_stride_apart():
    sets = kt.groups(1000, 128)
    assert sorted(np.concatenate(sets).tolist()) == list(range(1000))
    assert all((np.diff(s) == 128).all() for s in sets)
    assert len(kt.groups(5, 128)) == 5


class Recorder:
    """Stands in for the model: at a hidden place every value is the place's index in the window, so the reading it
    came from can be told; at every other place a sentinel. It checks what it is given: the hidden places a stride
    apart, unknowns after the genome, a length the model accepts."""
    def __init__(self, n, stride):
        self.n, self.stride, self.hidden = n, stride, []

    def __call__(self, rows):
        out = np.full((*rows.shape, 4), -99.0)
        for b, r in enumerate(rows):
            at = np.flatnonzero(r == kt.MASK)
            assert len(at) > 0 and (np.diff(at) == self.stride).all()
            assert (r[self.n:] == kt.N).all() and len(r) % 64 == 0
            self.hidden.append(at)
            out[b, at, :] = at[:, None].astype(float)
        return out


def test_one_reading_hides_every_place_once_and_keeps_only_hidden_values():
    n, rec = 300, Recorder(300, 50)
    codes = np.arange(n) % 4
    got = kt.one_reading(codes, rec, stride=50, multiple=64, batch=7)
    assert sorted(np.concatenate(rec.hidden).tolist()) == list(range(n))
    assert (got[:, 0] == np.arange(n)).all() and not (got == -99.0).any()


def test_an_unknown_place_is_not_scored():
    codes = np.arange(40) % 4
    codes[7] = kt.N
    got = kt.one_reading(codes, Recorder(40, 10), stride=10, multiple=64, batch=4)
    assert np.isnan(got[7]).all() and not np.isnan(np.delete(got, 7, axis=0)).any()


def test_each_place_comes_from_the_reading_where_it_is_farther_from_the_ends():
    n = 400
    got = kt.table(np.arange(n) % 4, Recorder(n, 50), stride=50, multiple=64, batch=8)
    # the stand-in returns the place's index in its own reading: in the record's order that is the place itself,
    # in the reading from the middle it is (place - 200) mod 400
    i = np.arange(n)
    own, turned = i, (i - 200) % n
    edge = lambda k: np.minimum(k, n - 1 - k)
    want = np.where(edge(own) >= edge(turned), own, turned)
    assert (got[:, 0] == want).all()
    assert got[0, 0] == 200 and got[399, 0] == 199 and got[200, 0] == 200      # the ends come from the other reading
    assert (edge(got[:, 0].astype(int)) >= n // 4 - 1).all()                   # no place nearer an end than a quarter


def test_a_model_answer_of_the_wrong_shape_is_refused():
    with pytest.raises(ValueError):
        kt.one_reading(np.arange(40) % 4, lambda rows: np.zeros((len(rows), rows.shape[1], 5)), 10, 64, 4)


def test_write(tmp_path):
    tab = np.array([[-0.1, -2.0, -3.0, -4.0], [np.nan] * 4])
    kt.write(tmp_path / "t.tsv", "aN", tab)
    assert (tmp_path / "t.tsv").read_text() == "pos\tref\tlp_A\tlp_C\tlp_G\tlp_T\n1\tA\t-0.10000\t-2.00000\t-3.00000\t-4.00000\n2\tN\t\t\t\t\n"
    with pytest.raises(ValueError):
        kt.write(tmp_path / "t.tsv", "a", tab)
