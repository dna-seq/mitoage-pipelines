"""What build 2 read, as pinned in the site repository (dna-seq/mitoage at the commit named in
docs/BUILD2_REPRODUCTION.md): pipeline/sources/pinned.ts, docs/BUILD2_REPORT.md ("Pinned upstream versions")
and the checksums of the copies its audit kept. A download is compared with these; the build-2 code never
reads a value from here that it could read from a source.
"""

# Hugging Face dataset longevity-genie/mitoage: the revision of build 2 the manuscript cites, with the EURING birds
# (the site's D15) and the hand-checked gene annotations (D17), written by dna-seq/mitoage at 03c0949. Earlier
# uploads: f75bdc9b812d7845d2ccb24a49643d88d00d57b5 (EURING, without the gene patch; site 725da08) and
# 653afe2a3b79e18ae892e579ded9acdb17e0c334 (the first, 2026-11-13 18:02 UTC, without EURING).
HF_REPO = "longevity-genie/mitoage"
HF_BUILD2_REVISION = "0f4e970554b2fe0265e718321f422347293dde17"
SITE_COMMIT = "03c09498564c17d1445bc4c18c70d0e710d471db"   # the site repository's commit that wrote that build

BUILD_DATE = "2026-11-13T00:00:00.000Z"       # build 2's date, written into added_at / changed_at
BUILD1_DATE = "2015-07-01T00:00:00.000Z"      # build 1's release date

# NCBI RefSeq mitochondrion release. NCBI keeps only the current release on its FTP site.
REFSEQ_URL = "https://ftp.ncbi.nlm.nih.gov/refseq/release/mitochondrion/mitochondrion.1.genomic.gbff.gz"
REFSEQ_RELEASE_URL = "https://ftp.ncbi.nlm.nih.gov/refseq/release/RELEASE_NUMBER"
REFSEQ_RELEASE = "237"
REFSEQ_FILE_DATE = "2026-09-03"
REFSEQ_MD5 = "4491799651740f49c8f6a011cf9a96e5"

EFETCH_URL = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"

# NCBI Taxonomy: build 2 read the current new_taxdump of 2026-11-13 (Last-Modified 2026-11-13T00:42:20Z).
# The monthly archive keeps the dump of the first of each month; the one closest before it is 2026-11-01.
TAXDUMP_CURRENT_URL = "https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/new_taxdump/new_taxdump.tar.gz"
TAXDUMP_ARCHIVE_URL = "https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump_archive/new_taxdump_{date}.zip"
TAXDUMP_BUILD2 = "2026-11-13"
TAXDUMP_BUILD2_MD5 = "91f95787dd7c7cfd7042ab7610f2b4b2"      # new_taxdump.tar.gz as build 2 downloaded it
TAXDUMP_ARCHIVE_CLOSEST = "2026-11-01"
TAXDUMP_FILES = ("nodes.dmp", "names.dmp", "merged.dmp", "rankedlineage.dmp")

# AnAge (unversioned URL: a new AnAge build replaces the file)
ANAGE_URL = "https://genomics.senescence.info/species/dataset.zip"
ANAGE_VERSION = "AnAge build 15 (2023-07-03)"
ANAGE_MD5 = "7ec4d19f6c8f57b22595a9ffd50ea66b"

# FishBase and SeaLifeBase parquet snapshots on source.coop (cboettig/fishbase), CC BY-NC 4.0
SOURCECOOP_URL = "https://data.source.coop/cboettig/fishbase/{db}/{version}/parquet/{table}.parquet"
FISHBASE_VERSION = "v26.06"
FISHBASE_TABLES = ("species", "popchar", "popgrowth", "synonyms")   # the tables the lifespan rules read
FISHBASE_MD5 = {   # md5 of the files build 2 read (its audit's md5sums.txt)
    "fb": {"species": "0208d896c3a75868c1650c04a2eb4330", "popchar": "ad782577a255e428d3b4b57210dd0092",
           "popgrowth": "953a4110973c30b56472cc3ed4b9ce22", "synonyms": "42d92d32ee7e71d8b2dd899b7eb69f9e"},
    "slb": {"species": "f176abbbc62ee9e22451be155f75f06a", "popchar": "dcb288862abd87945f314b9a3df409f8",
            "popgrowth": "9f08f077586a1ce9417e90c3aa61a6bb", "synonyms": "31c5089e191915f881d687322251e9a1"},
}
FISHBASE_LABEL = {"fb": "FishBase v26.06 (source.coop cboettig/fishbase, files 2026-09-04)",
                  "slb": "SeaLifeBase v26.06 (source.coop cboettig/fishbase, files 2026-09-04)"}

# USGS Bird Banding Laboratory longevity records (unversioned page; the table says its version)
USGS_URL = "https://www.pwrc.usgs.gov/BBL/Bander_Portal/login/Longevity.php"
USGS_VERSION = "NABBP Longevity Records v2025.1 (current through Dec 2025)"
USGS_HTML_MD5 = "b024376be0641b560c2a9807806df68b"   # the page build 2's audit saved
USGS_CSV_MD5 = "a84743335e7f080d5891c2f719f53702"    # the table build 2 read, as CSV

# EURING longevity list of 2023-09-01. euring.org shows a maintenance page (since 2026-11-13 at least), so build 2
# reads the PDF the Internet Archive captured on 2026-10-11, the same file (one SHA-1) as every capture since
# 2023-12-16 (the site's decision D15). Version, licence and reference as the site's pipeline writes them.
EURING_URL = "https://euring.org/data-and-codes/longevity-list"
EURING_PDF_URL = ("https://web.archive.org/web/20261011111048id_/"
                  "https://euring.org/files/documents/EURING_longevity_list_20230901.pdf")
EURING_PDF = "EURING_longevity_list_20230901.pdf"
EURING_MD5 = "c82b37efe23fe95b5e2945f894c46496"
EURING_VERSION = "list of 2023-09-01 (EURING_longevity_list_20230901.pdf, Internet Archive capture of 2026-10-11)"
EURING_LICENSE = ("no licence stated on the list or its page (euring.org/data-and-codes/longevity-list-2023); "
                  "cite Fransson, Kolehmainen, Moss & Robinson (2023) EURING list of longevity records for European birds")
EURING_REFERENCE = "Fransson, Kolehmainen, Moss & Robinson (2023) EURING list of longevity records for European birds"

# The old site (July 2015), for the old species ids
LEGACY_URL = "https://www.mitoage.info/export-section/total_mtDNA/for-all-species/?page={page}"

# The dated tree
OTL_API = "https://api.opentreeoflife.org/v3"
TIMETREE_API = "https://api.timetree.org"
OTL_SYNTH_BUILD2 = "opentree16.1"

# Alignments
MAFFT_VERSION = "v7.526"

# Licences as build 2 states them per row (lifespan_license, license)
ANAGE_LICENSE = ("Free to use, subject to the Human Ageing Genomic Resources conditions of use "
                 "(https://genomics.senescence.info/legal.html); cite Tacutu et al. 2013, Nucleic Acids Res 41:D1027-D1033")
FISHBASE_LICENSE = {
    "fb": "CC BY-NC 4.0 (FishBase.org team; source.coop cboettig/fishbase dataset card): non-commercial use only; cite Froese & Pauly, FishBase",
    "slb": "CC BY-NC 4.0 (SeaLifeBase/FishBase.org team; source.coop cboettig/fishbase dataset card): non-commercial use only; cite Palomares & Pauly, SeaLifeBase",
}
USGS_LICENSE = ("U.S. Government work (USGS), presumed public domain; cite Bird Banding Laboratory, NABBP Longevity Records "
                "v2025.1, USGS EESC, Laurel MD (usgs.gov terms page returned HTTP 403, not verified)")
BUILD1_LICENSE = "Free for all purposes provided MitoAge is cited (Toren et al. 2016, Nucleic Acids Res 44:D1262-D1265)"
REFSEQ_VERSION = "release 237 (mitochondrion.1.genomic.gbff.gz, file date 2026-09-03)"
TAXONOMY_VERSION = "new_taxdump of 2026-11-13 (Last-Modified 2026-11-13T00:42:20Z)"

# The snapshot dataset of the inputs no upstream archive keeps (`--from-snapshot`)
SNAPSHOT_REPO = "longevity-genie/mitoage-build2-sources"
SNAPSHOT_REVISION = "a76215144ec33d6ed1e5ea599cb797f2651dab1a"   # 2026-11-15: the tree answers of build 2 with the EURING birds (first upload 9d55d72, 2026-11-14)
