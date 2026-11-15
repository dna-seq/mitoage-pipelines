"""The EURING longevity list reader (build2/euring.py): ages, the table's structure, and the pinned PDF itself
when `mitoage download lifespan` has put it in the data folder (the file states no licence, so it is not
copied into tests/data)."""
import pytest

from mitoage_pipelines import paths
from mitoage_pipelines.build2 import euring, pins
from mitoage_pipelines.build2.euring import Run

PDF = paths.B2_EURING / pins.EURING_PDF


class TestAge:
    def test_years_and_months(self):
        assert euring.age("37y 0m") == (37.0, None)
        assert euring.age("30y 9m")[0] == pytest.approx(30.75, abs=1e-12)
        assert euring.age("12 y 11m")[0] == 12 + 11 / 12          # the same double the site's code computes
        assert euring.age("10y 7m ")[0] == 10 + 7 / 12

    def test_an_unreadable_month_keeps_the_years_and_says_so(self):
        years, note = euring.age("11y im")
        assert years == 11 and '"im"' in note
        assert euring.age("5y 12m")[0] == 5
        assert euring.age("about 9 years") == (None, 'age "about 9 years" not readable')


HEADER = [Run(1, 24, 645, "Arial-BoldMT", "Species / Scientific name"), Run(1, 365, 645, "Arial-BoldMT", "Age")]
QUAIL = HEADER + [
    Run(1, 24, 624, "Arial-BoldMT", "Quail "), Run(1, 199, 624, "ArialMT", "Israel"), Run(1, 261, 624, "ArialMT", "ILT"),
    Run(1, 303, 624, "ArialMT", "C 7154"), Run(1, 359, 624, "ArialMT", "12y 1m"), Run(1, 449, 624, "ArialMT", "Found dead"),
    Run(1, 24, 610.8, "Arial-ItalicMT", "Coturnix coturnix"), Run(1, 199, 610.7, "ArialMT", "Spain"), Run(1, 261, 610.7, "ArialMT", "ESI"),
    Run(1, 303, 610.7, "ArialMT", "3303820"), Run(1, 359, 610.7, "ArialMT", "11y im"), Run(1, 463, 610.7, "ArialMT", "Shot"),
]


class TestStructure:
    def test_two_records_of_one_species_in_any_drawing_order(self):
        recs = euring.parse(QUAIL)
        assert euring.parse(list(reversed(QUAIL))) == recs
        assert [(r.scientific_name, r.rank, r.country, r.ring_no, r.years) for r in recs] == [
            ("Coturnix coturnix", 1, "Israel", "C 7154", 12 + 1 / 12), ("Coturnix coturnix", 2, "Spain", "3303820", 11.0)]

    def test_refusals(self):
        with pytest.raises(ValueError, match="no scientific name"):
            euring.parse(QUAIL[:8])
        with pytest.raises(ValueError, match="follows no common name"):
            euring.parse(HEADER + [Run(1, 24, 610, "Arial-ItalicMT", "Coturnix coturnix")])
        with pytest.raises(ValueError, match="no table header"):
            euring.parse(QUAIL[2:])


@pytest.fixture(scope="module")
def runs():
    return euring.text_runs(PDF.read_bytes())


@pytest.mark.skipif(not PDF.is_file(), reason="the EURING PDF is not in the data folder (mitoage download lifespan)")
class TestPinnedList:

    def test_every_age_cell_becomes_one_record(self, runs):
        import re
        cells = sorted(euring.collapse(r.text) for r in runs if r.font == "ArialMT" and re.fullmatch(r"\s*\d+\s*y\s*\S+\s*", r.text))
        recs = euring.parse(runs)
        assert sorted(r.age_text for r in recs) == cells
        assert any(r.text == "European Longevity Records - 2023-09-01" for r in runs)

    def test_one_oldest_record_per_species_and_the_second_never_older(self, runs):
        recs = euring.parse(runs)
        by = {}
        for r in recs: by.setdefault(r.scientific_name, []).append(r)
        for name, g in by.items():
            assert [r.rank for r in g] in ([1], [1, 2]), name
            if len(g) == 2: assert g[1].years <= g[0].years, name
        assert sum(r.rank == 1 for r in recs) == len(by)

    def test_observations_carry_the_pin(self):
        obs = euring.load(PDF)
        assert len({o.source_record_id for o in obs}) == len(obs)
        assert {o.source_version for o in obs} == {pins.EURING_VERSION}
        stork = [o for o in obs if o.source_name == "Ciconia ciconia"]
        assert [o.source_common_name for o in stork] == ["White Stork", "White Stork"]
