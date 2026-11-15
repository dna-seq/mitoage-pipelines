"""The hand-checked annotation patches and gene gaps of build 2 (build2/gene_annotation.py, the site's D17):
translating a patched frame, refusing a bad one, and, on the pinned RefSeq release when it is in the data folder
(`mitoage download refseq --from-snapshot`), the African grey parrot's ND6 and every recorded gap."""
import gzip

import pytest

from mitoage_pipelines import paths
from mitoage_pipelines.build2 import genome as gm
from mitoage_pipelines.build2.gene_annotation import ANNOTATION_PATCHES, GENE_GAPS, AnnotationPatch

RELEASE = paths.refseq_release("237")


class TestTranslatePatch:
    def test_a_clean_frame_reads_with_the_start_as_M_and_the_stop_dropped(self):
        # ATA is isoleucine in the standard code and methionine in code 2; as a start it reads M in both
        assert gm.translate_patch("ATAGCTTGGTAA", 2, "t") == "MAW"
        assert gm.translate_patch("ATTAAATAG", 1, "t") == "MK"

    def test_refuses_part_codons_an_internal_stop_and_a_missing_stop(self):
        with pytest.raises(ValueError, match="whole codons"): gm.translate_patch("ATGGC", 2, "t")
        with pytest.raises(ValueError, match="stop codon at codon 2 of 3"): gm.translate_patch("ATGTAAGCT", 2, "t")
        with pytest.raises(ValueError, match="no stop codon"): gm.translate_patch("ATGGCTGCT", 2, "t")
        with pytest.raises(ValueError, match="not a plain codon"): gm.translate_patch("ATGNCTTAA", 2, "t")


def _records(accessions: set[str]) -> dict[str, str]:
    """The records of these accession.versions from the pinned release, by streaming it once."""
    out, buf, version = {}, [], None
    with gzip.open(RELEASE, "rt") as f:
        for line in f:
            if line.startswith("LOCUS"):
                buf, version = [line], None
                continue
            buf.append(line)
            if line.startswith("VERSION"):
                version = line.split()[1]
            elif line.startswith("//") and version in accessions:
                out[version] = "".join(buf)
                if len(out) == len(accessions): break
    return out


@pytest.fixture(scope="module")
def records():
    return _records({p.accession for p in ANNOTATION_PATCHES} | {g.accession for g in GENE_GAPS})


@pytest.mark.skipif(not RELEASE.is_file(), reason="RefSeq release 237 is not in the data folder")
class TestOnTheRelease:
    def test_the_parrot_ND6_is_added_and_reads_like_ND6(self, records):
        g = gm.annotate(records["NC_027847.1"], 2, True)
        assert g.patched == ["ND6"]
        nd6 = [f for f in g.features if f.name == "ND6" and f.type == "CDS"]
        assert len(nd6) == 1 and (nd6[0].start, nd6[0].end, nd6[0].strand) == (17376, 17891, -1)
        cds = next(f for f in g.record.features if f.key == "CDS" and gm.feature_name(f)[0] == "ND6")
        protein = cds.qualifiers["translation"][0]
        nt = gm.extract(g.record.sequence, cds.segments)
        assert len(nt) == 17891 - 17376 + 1 == 3 * (171 + 1)           # 171 codons and the TAA stop, which is dropped
        assert nt[-3:] == "TAA" and len(protein) == 171 and protein.startswith("M")

    def test_a_patch_is_refused_when_the_record_has_the_gene_or_the_frame_is_wrong(self, records):
        rec = gm.parse_genbank(records["NC_027847.1"])
        gm.apply_gene_patches(rec, 2)
        with pytest.raises(ValueError, match="annotates ND6 itself"):
            gm.apply_gene_patches(rec, 2)                                  # the record now has the gene
        shifted = AnnotationPatch("NC_027847.1", "ND6", "complement(17377..17891)", "one base short")
        with pytest.raises(ValueError, match="whole codons"):
            gm.apply_gene_patches(gm.parse_genbank(records["NC_027847.1"]), 2, [shifted])

    def test_every_recorded_gap_is_still_a_gap_in_its_record(self, records):
        for gap in GENE_GAPS:
            g = gm.annotate(records[gap.accession], 2, True)
            assert gap.gene not in {f.name for f in g.features if f.type == "CDS"}, gap
            if gap.reason == "frameshifted_in_record":
                notes = [" ".join(f.qualifiers.get("note", []) + f.qualifiers.get("product", [])) for f in g.record.features if f.key == "misc_feature"]
                assert any("frameshift" in n.lower() or "similar to" in n.lower() for n in notes), gap
